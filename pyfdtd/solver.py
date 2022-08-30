"""Main module."""
import time

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle
from numpy import (
    add,
    array,
    diff,
    divide,
    exp,
    floor,
    linspace,
    multiply,
    ones,
    sqrt,
    subtract,
    zeros,
)


def _initialize(Nx, Ny, dx, dy, Nt, df, l, b, NPML, ur, er, freq, epssrc, musrc, src):
    start = time.time()
    e0 = 8.85e-12
    Nx2 = 2 * Nx
    Ny2 = 2 * Ny
    epssrc = 1
    musrc = 1

    URxx, URyy = _build_geometry(dx, Nx, l, b, dy, Ny, musrc, ur)

    Ezsrc, Hxsrc, dt, steps = _build_source(Nt, freq, df, musrc, epssrc, dy)

    sigx, sigy = _build_PML(Nx2, Ny2, e0, NPML, dt)

    PML_coefficients = _build_PML_coefficients(
        sigx, sigy, dt, e0, 299792458, URxx, URyy, Nx2, Ny2
    )

    Ez = _solve(*PML_coefficients, Nx, Ny, steps, dx, dy, src, Ezsrc, Hxsrc)

    plot_E(Ez, l, b, Nx, Ny)
    end = time.time()
    print(f"Time elapsed: {end - start}")


def _build_geometry(dx, Nx, l, b, dy, Ny, musrc, ur):
    nx = round(l / dx)
    nx1 = int(1 + floor((Nx - nx) / 2))
    nx2 = nx1 + nx - 1

    ny = round(b / dy)
    ny1 = int(1 + floor((Ny - ny) / 2))
    ny2 = ny1 + ny - 1

    UR2 = zeros((Nx, Ny))
    UR2[nx1 - 1 : nx2, ny1 - 1 : ny2] = 1
    UR2 = musrc * (ones(UR2.shape) - UR2) + ur * UR2
    URxx = UR2
    URyy = UR2

    return URxx, URyy


def _build_source(Nt, freq, df, musrc, epssrc, dy):
    dt = 1 / (Nt * freq)
    tau = 0.5 / freq
    t0 = 10 * tau
    steps = 1 / (dt * df)
    t = dt * linspace(0, steps, int(steps), endpoint=False)
    Ezsrc = exp(-(((t - t0) / tau) ** 2))
    nsrc = sqrt(musrc * epssrc)
    c0 = 299792458
    diract = (nsrc * dy / (2 * c0)) + dt / 2
    Hxsrc = sqrt(epssrc / musrc) * exp(-(((t + diract - t0) / tau) ** 2))

    return Ezsrc, Hxsrc, dt, steps


def _build_PML(Nx2, Ny2, e0, NPML, dt):
    sigx = zeros((Nx2, Ny2))
    for nx in range((2 * NPML[0])):
        nx1 = 2 * NPML[0] - nx + 1
        sigx[nx1][:] = (0.5 * e0 / dt) * (nx / 2 / NPML[0]) ** 3

    for nx in range((2 * NPML[1])):
        nx1 = Nx2 - 2 * NPML[1] + nx
        sigx[nx1][:] = (0.5 * e0 / dt) * (nx / 2 / NPML[1]) ** 3

    sigy = zeros((Nx2, Ny2))
    for ny in range(2 * NPML[2]):
        ny1 = 2 * NPML[2] - ny + 1
        sigy[:][ny1] = (0.5 * e0 / dt) * (ny / 2 / NPML[2]) ** 3

    for ny in range(2 * NPML[3]):
        ny1 = Ny2 - 2 * NPML[3] + ny
        sigy[:][ny1] = (0.5 * e0 / dt) * (ny / 2 / NPML[3]) ** 3

    return sigx, sigy


def _build_PML_coefficients(sigx, sigy, dt, e0, c0, URxx, URyy, Nx2, Ny2):
    sigHx = sigx[0:Nx2:2, 1:Ny2:2]
    sigHy = sigy[0:Nx2:2, 1:Ny2:2]
    mHx0 = (1 / dt) + sigHy / (2 * e0)
    mHx1 = ((1 / dt) - sigHy / (2 * e0)) / mHx0
    mHx2 = -c0 / URxx / mHx0
    mHx3 = -(c0 * dt / e0) * sigHx / URxx / mHx0
    sigHx = sigx[1:Nx2:2, 0:Ny2:2]
    sigHy = sigy[1:Nx2:2, 0:Ny2:2]
    mHy0 = (1 / dt) + sigHx / (2 * e0)
    mHy1 = ((1 / dt) - sigHx / (2 * e0)) / mHy0
    mHy2 = -c0 / URyy / mHy0
    mHy3 = -(c0 * dt / e0) * sigHy / URyy / mHy0
    sigDx = sigx[0:Nx2:2, 0:Ny2:2]
    sigDy = sigy[0:Nx2:2, 0:2:Ny2]
    mDz0 = (1 / dt) + (sigDx + sigDy) / (2 * e0) + sigDx * sigDy * (dt / 4 / e0**2)
    mDz1 = (1 / dt) - (sigDx + sigDy) / (2 * e0) - sigDx * sigDy * (dt / 4 / e0**2)
    mDz1 = mDz1 / mDz0
    mDz2 = c0 / mDz0
    mDz4 = -(dt / e0**2) * sigDx * sigDy / mDz0

    return (
        mHx1,
        mHx2,
        mHx3,
        mHy1,
        mHy2,
        mHy3,
        mDz1,
        mDz2,
        mDz4,
    )


def _solve(
    mHx1,
    mHx2,
    mHx3,
    mHy1,
    mHy2,
    mHy3,
    mDz1,
    mDz2,
    mDz4,
    Nx,
    Ny,
    steps,
    dx,
    dy,
    src,
    Ezsrc,
    Hxsrc,
):
    CEx = zeros((Nx, Ny))
    CEy = zeros((Nx, Ny))
    Ez = zeros((Nx, Ny))
    Hx = zeros((Nx, Ny))
    Hy = zeros((Nx, Ny))
    CHz = zeros((Nx, Ny))
    Dz = zeros((Nx, Ny))
    ICEx = zeros((Nx, Ny))
    ICEy = zeros((Nx, Ny))
    IDz = zeros((Nx, Ny))
    E = []

    for T in range(int(steps)):
        # print(f'Iteration {T}')

        # Find Curl of Ex and Ey
        CEx[:, : Ny - 1] = divide(diff(Ez), dy)
        CEx[:, -1] = divide(-Ez[:, -1], dy)

        # Inject source to the curl of E
        CEx[1:, src - 1] = (
            divide(subtract(Ez[1:, src], Ez[1:, src - 1]), dy) - Ezsrc[T] / dy
        )
        CEy[: Nx - 1, :] = divide(-(diff(Ez, axis=0)), dx)
        CEy[-1, :] = divide(Ez[-1, :], dx)

        ICEx = add(ICEx, CEx)
        ICEy = add(ICEy, CEy)

        Hx = add(multiply(mHx1, Hx), multiply(mHx2, CEx), multiply(mHx3, ICEx))
        Hy = add(multiply(mHy1, Hy), multiply(mHy2, CEy), multiply(mHy3, ICEy))

        CHz[0, 0] = subtract(divide(Hy[0, 0], dx), divide(Hx[0, 0], dy))
        CHz[1:, 0] = subtract(divide(diff(Hy[:, 0], axis=0), dx), divide(Hx[1:, 0], dy))
        CHz[0, 1:] = add(
            divide(-diff(Hx[0, :].reshape(1, -1), axis=1), dy), divide(Hy[0, 1:], dx)
        )
        CHz[1:, 1:] = subtract(
            divide(diff(Hy[:, 0], axis=0), dx),
            divide(diff(Hx[0, :].reshape(1, -1), axis=1), dx),
        )
        # for ny in range(2, Ny):
        #     # CHz[1,ny] = (Hy[1,ny])/dx - (Hx[1,ny] - Hx[1,ny-1])/dy
        #     for nx in range(2, Nx):
        #         CHz[nx,ny] = (Hy[nx,ny] - Hy[nx-1,ny])/dx - (Hx[nx,ny] - Hx[nx,ny-1])/dy

        CHz[1:, src] = (
            subtract(
                divide(diff(Hy[:, src], axis=0), dx),
                divide(subtract(Hx[1:, src], Hx[1:, src - 1]), dy),
            )
            + Hxsrc[T] / dy
        )

        IDz = add(IDz, Dz)

        Dz = add(multiply(mDz1, Dz), multiply(mDz2, CHz), multiply(mDz4, IDz))

        Ez = array(multiply(mDz1, Dz))

        E.append(Ez)

    return E


def plot_E(Ez, l, b, Nx, Ny):
    fig, ax = plt.subplots()

    def animate(i):
        ax.clear()
        ax.imshow(Ez[i], cmap="RdBu", vmin=-0.1, vmax=0.1)
        ax.set_title("Ez at time step {}".format(i))
        ax.set_xlabel("x")
        ax.set_ylabel("y")

    ax.add_patch(
        Rectangle(
            width=l,
            height=b,
            xy=(Nx / 2.0, Ny / 2.0),
            linewidth=1,
            edgecolor="black",
            facecolor="none",
        )
    )
    anim = FuncAnimation(fig, animate, frames=len(Ez), interval=10, repeat=False)
    # anim.save("Ez.gif")

    plt.show()


# if __name__ == "__main__":
    # _initialize(
    #     1000, 1000, 10, 10, 10, 1e6, 200, 200, [20, 20, 20, 20], 2, 6, 1e8, 1, 1, 380
    # )
