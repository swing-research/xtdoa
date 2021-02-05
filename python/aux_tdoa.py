"""
Copyright 2019 Robin Scheibler

Code for the computation of the sub-sample delay between two signals.
"""
import numpy as np

import pyroomacoustics as pra


def parabola(x, ref, phat=False, interp=1, fs=1):

    xcorr = pra.correlate(x, ref, interp=interp, phat=phat)

    d = np.argmax(np.abs(xcorr))

    # Fit parabola using 3 points
    y1, y2, y3 = xcorr[d - 1 : d + 2]
    a = 0.5 * (y1 - 2 * y2 + y3)
    b = 0.5 * (-3 * y1 + 4 * y2 - y3)
    # for completeness
    # c = y1

    # Now the parabola max is here:
    d_frac = -0.5 * b / a + (d - 1)

    # give proper scale
    return (d_frac / interp - (ref.shape[0] - 1)) / fs


def aux_tdoa(x, ref, n_iter=10, t0=None, phat=False, f_init=parabola, f_init_kwargs={}):
    """
    Auxilliary function based algorithm for the estimation of sub-sample precision time
    delay between two signals described in the paper:

    K. Yamaoka, R. Scheibler, N. Ono, and Y. Wakabayashi, “Sub-Sample Time
    Delay Estimation via Auxiliary-Function-Based Iterative Updates,” Proc.
    WASPAA, New Paltz, NY, USA, pp. 130–134, Oct. 2019.

    Parameters
    ----------
    x: array_like (shape: n_samples)
        The delayed signal
    ref: array_like (shape: n_samples)
        The signal with respect to which the delay is computed
    n_iter: int, optional
        The number of iterations
    t0: float, optional
        The initial estimate of the time-delay. If no estimate is provided, GCC
        with quadratic interpolation is used.
    phat: bool, optional
        If set to ``True``, the PHAT weighting is applied
    f_init: func, optional
        A function to use for initial estimation of the time delay, default is
        GCC with quadratic interpolation.  The signature of the function should
        be ``func(x, ref, **kwargs)``.  The `parabola` method gives great
        results by fitting a quadratic function on three points around the
        maximum of GCC or GCC-PHAT. Another simple choice is
        `pyroomacoustics.tdoa`, which is just GCC-PHAT.
    f_init_kwargs: dict, optional
        Keyword arguments for the initial estimation function

    Returns
    -------
    The sub-sample precision estimate of the time-delay of ``x`` with respect to ``ref``.
    """

    N1 = x.shape[0]
    N2 = ref.shape[0]

    N = N1 + N2 - 1

    X1 = np.fft.rfft(x, n=N)
    X2 = np.fft.rfft(ref, n=N)

    if phat:
        X1 /= np.abs(X1) + 1e-15
        X2 /= np.abs(X2) + 1e-15

    S = X1 * np.conj(X2)
    A = np.abs(S)
    phi = np.angle(S / A)
    A /= N
    A[1:-1] *= 2
    w = 2.0 * np.pi * np.arange(N // 2 + 1) / N

    Aw = A * w
    Aw2 = A * w ** 2

    if t0 is None:
        t = f_init(x, ref, **f_init_kwargs)
    else:
        t = t0

    for epoch in range(n_iter):
        # new phase estimate
        theta = w * t + phi

        # round to within 2 * np.pi
        n_2pi = np.round(theta / 2 / np.pi)
        theta -= n_2pi * 2 * np.pi

        # update delay estimate
        t -= np.sum(Aw * np.sin(theta)) / np.sum(Aw2 * np.sinc(theta / np.pi))

    return t
