"""Console script for pyfdtd."""
import argparse
import ast
import sys

from pyfdtd import solver


def main():
    """Console script for pyfdtd."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", "-F", type=str, required=True)
    args = parser.parse_args()

    with open(args.file, "r") as f:
        lines = f.read().split("\n")
        Nx = int(lines[0])
        Ny = int(lines[1])
        dx = float(lines[2])
        dy = float(lines[3])
        Nt = int(lines[4])
        df = float(lines[5])
        l = float(lines[6])
        b = float(lines[7])
        _npml = ast.literal_eval(lines[8])
        NPML = [_npml[0], _npml[1], _npml[2], _npml[3]]
        ur = float(lines[9])
        er = float(lines[10])
        freq = float(lines[11])
        epssrc = float(lines[12])
        musrc = float(lines[13])
        src = int(lines[14])
        solver._initialize(
            Nx, Ny, dx, dy, Nt, df, l, b, NPML, ur, er, freq, epssrc, musrc, src
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
