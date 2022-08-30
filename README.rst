========
pyfdtd
========

Simple 2D FDTD Solver written in Python. This is a very simple implementation, meant entirely for the purpose of me learning the underlying theory of the Finite Difference Time Domain method, while improving my Python skills.
Currently, it supports arbitrary sized rectangles as geometries, with real valued permittivity and permeability values (EXPERIMENTAL PHASE). There is still a lot of work to be done.

To Test this code out:
1. Fork this repo and cd into it
2. run `dist/pyfdtd -F input.txt` or `python pyfdtd/solver.py`
3. For 1000 runs, it should take about a minute and a half, after which the code opens a matplotlib window with the E field movie.

* Free software: MIT license


TODO
----
Extend geometries to arbitrary shapes. Perhaps add a feature which lets the user import STL files or GDSII files.

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage

References
----------

I am following the FDTD in MATLAB course taught by Prof. Raymond Rumpf at UT El Paso:

1.  "Electromagnetic Analysis Using Finite-Difference Time-Domain," EMPossible.
