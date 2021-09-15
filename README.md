# The Mittag-Leffler function in Python

This package contains modules for calculating and approximating the Mittag-Leffler function.

The `mittag_leffler.py` module contains a Python port of
a published
[Matlab implementation](https://se.mathworks.com/matlabcentral/fileexchange/48154-the-mittag-leffler-function) of
the generalized Mittag-Leffler function, written by Konrad Hinsen. The module `ml_internal.py` contains the internal functions for `mittag_leffler.py`. The script `test_ml.py` contains tests for the functions in `mittag_leffler.py`. They cover very little of the total functionality of
the code. To use the Python port of the Matlab implementation, simply import `ml` from `mittag_leffler.py`.

The `pade_approx.py` module contains a Python implementation of a global Pade approximation for the Mittag-Leffler function as described by [Sarumi, Furati, and Khaliq](https://arxiv.org/abs/1912.10996). This approximation is valid for `z < 0`, `0 < alpha < 1`, `beta >= alpha`, `(alpha, beta) != (1, 1)`. The approximaton is highly accurate and requires 3 orders of magnitude less computation time than the Matlab algorithm.

