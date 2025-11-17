import numpy as np
import scipy


def gaussian_ecf_radius(ecf, radius=1.0):
    """
    Return the radius that encloses the given count fraction.
    """
    # the enclosed fraction for a 2d Gaussian with sigma=1 is
    # ecf = 1 - exp(-r^2/2)
    return np.sqrt(-2 * np.log(1 - ecf)) * radius


def gaussian_inverse_covariance(sigma_1, sigma_2, angle):
    """
    Return the inverse covariance matrix for a 2d Gaussian with given sigmas and rotation angle.
    """
    diag = np.diag([1.0 / sigma_1**2, 1.0 / sigma_2**2])
    transform = np.array(
        [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
    )
    return transform @ diag @ transform.T


def log_normal_prob_2d(x, x0, x1, sigma_1, sigma_2, angle):
    """
    Return the log of the probability density of a 2d Gaussian with given sigmas and rotation angle.
    """
    if sigma_1 <= 0 or sigma_2 <= 0:
        return -np.inf * np.ones(x.shape[0])
    norm = 1 / (2 * np.pi * sigma_1 * sigma_2)
    loc = np.array([[x0, x1]])
    cov = gaussian_inverse_covariance(sigma_1, sigma_2, angle)
    exponent = -0.5 * np.einsum("ij,jk,ik->i", x - loc, cov, x - loc)
    return np.log(norm) + exponent


def normal_prob_2d(x, x0, x1, sigma_1, sigma_2, angle):
    """
    Return the probability density of a 2d Gaussian with given sigmas and rotation angle.
    """
    if sigma_1 <= 0 or sigma_2 <= 0:
        return np.zeros(x.shape[0])
    norm = 1 / (2 * np.pi * sigma_1 * sigma_2)
    loc = np.array([[x0, x1]])
    cov = gaussian_inverse_covariance(sigma_1, sigma_2, angle)
    exponent = -0.5 * np.einsum("ij,jk,ik->i", x - loc, cov, x - loc)
    exp = np.exp(exponent)
    return norm * exp


def normal_prob_1d(x, x0, sigma):
    """
    Return the probability density of a 1d Gaussian with given sigma.
    """
    if sigma <= 0:
        return np.zeros(x.shape[0])
    norm = 1 / np.sqrt(2 * np.pi) / sigma
    exponent = -0.5 * ((x - x0) / sigma) ** 2
    exp = np.exp(exponent)
    return norm * exp


def nan_log(x):
    """
    Safely calculate the natural logarithm, returning -inf for non-positive values.
    """
    x = np.asanyarray(x)
    res = np.full(x.shape, -np.inf)
    mask = x > 0
    res[mask] = np.log(x[mask])
    return res


def p_uniform(x, box_size=4):
    """
    Return the probability density of a uniform distribution over a box of given size.
    """
    return 1 / (4 * box_size**2) * np.ones(x.shape[0])


class Likelihood:
    def __init__(self, data, box_size=4):
        self.data = data
        self.box_size = box_size

    def __call__(self, x):
        x0, x1, sigma_1, sigma_2, rho, snr = x

        if sigma_1 <= 0 or sigma_2 <= 0 or snr < 0:
            return np.inf
        s = snr / (1 + snr)
        b = 1 / (1 + snr)
        p = s * normal_prob_2d(
            self.data, x0, x1, sigma_1, sigma_2, rho
        ) + b * p_uniform(self.data, box_size=self.box_size)
        res = -nan_log(p).sum()
        return res


def _fit(events, source, columns=("y_angle", "z_angle"), box_size=4):
    """
    Fit a 2d Gaussian to the given events around the given source position.

    Parameters
    ----------
    events : np.ndarray
        The events to fit.
    source : dict
        The source position to fit around.
    columns : tuple of str
        The columns to use for the fit. Both the events and source must have these columns.
    """
    data = np.vstack([events[columns[0]], events[columns[1]]]).T

    result = scipy.optimize.minimize(
        Likelihood(data, box_size=box_size),
        x0=[source[columns[0]], source[columns[1]], 1, 1, 0, 3],
        bounds=[
            (source[columns[0]] - 10, source[columns[0]] + 10),
            (source[columns[1]] - 10, source[columns[1]] + 10),
            (0.1, 20),
            (0.1, 20),
            (-np.pi / 2, np.pi / 2),
            (0.1, 1000),
        ],
    )
    if not result.success:
        raise RuntimeError("Fit did not converge: " + result.message)

    return result


def fit_gaussian_2d(events, source, columns=("y_angle", "z_angle"), box_size=4):
    """
    Fit a 2d Gaussian to the given events around the given source position.

    Parameters
    ----------
    events : np.ndarray
        The events to fit.
    source : dict
        The source position to fit around.
    columns : tuple of str
        The columns to use for the fit. Both the events and source must have these columns.
    """
    # enclosed count fraction (ECF) to use for determining source extent and signal to noise ratio
    ecf = 0.9
    # sigma_ecf is the Mahalanobis distance that encloses the given fraction of counts (ECF)
    sigma_ecf = gaussian_ecf_radius(ecf)

    fail_value = {
        "params": np.zeros(6),
        "hess_inv": 1e10 * np.eye(6),
        "ndof": 0,
        "fit_ok": False,
        "p_signal": 0,
        columns[0]: source[columns[0]],
        columns[1]: source[columns[1]],
        "sigma": (np.nan, np.nan),
        "rot_angle": np.nan,
        f"sigma_{columns[0]}": np.inf,
        f"sigma_{columns[1]}": np.inf,
        f"corr_{columns[0]}_{columns[1]}": 0.0,
        "psf_ratio": np.inf,
        "source_area": np.inf,
        "n": 0,
        "signal": 0,
        "background": 0,
        "snr": np.nan,
        f"ks_{columns[0]}": np.nan,
        f"ks_{columns[1]}": np.nan,
        f"ks_p_value_{columns[0]}": np.nan,
        f"ks_p_value_{columns[1]}": np.nan,
        f"ks_sign_{columns[0]}": np.nan,
        f"ks_sign_{columns[1]}": np.nan,
        "COMPONENT": source["COMPONENT"],
    }
    if len(events) < 10:
        return fail_value

    result = _fit(events, source, columns=columns, box_size=box_size)
    if not result.success:
        return fail_value

    p_signal = result.x[5] / (1 + result.x[5])
    p_bkg = 1 / (1 + result.x[5])
    inv_cov = gaussian_inverse_covariance(result.x[2], result.x[3], result.x[4])
    cov = scipy.linalg.inv(inv_cov)
    # these are the widths of the marginal distributions for x and y
    sigma_x = np.sqrt(cov[0, 0])
    sigma_y = np.sqrt(cov[1, 1])

    # these give the extent of the distribution along the principal axes for the requested ECF
    # (the result of the fit are the "1-sigma" widths, where the Mahalanobis distance is 1)
    sigma_1 = sigma_ecf * result.x[2]
    sigma_2 = sigma_ecf * result.x[3]
    source_area = np.pi * sigma_1 * sigma_2

    psf_ratio = np.sqrt(sigma_1 * sigma_2) / sigma_ecf

    # the actual number of events within the source extent (sigma_ecf)
    x = np.vstack([events[columns[0]], events[columns[1]]]).T
    x0 = np.array([result.x[0], result.x[1]])
    n = np.count_nonzero(
        np.einsum("ij,jk,ik->i", x - x0, inv_cov, x - x0) < sigma_ecf**2
    )

    # to estimate the signal to noise ratio, we calculate the fraction of signal and background
    # counts within the source extent (the integral of the signal and background distributions).
    # we then multiply that by the total number of counts, and calculate the SNR as:
    #
    #   signal / sqrt(background).

    if n == 0:
        return fail_value

    # by definition: the fraction of signal counts within sigma_ecf is given by ECF
    f_signal = p_signal * ecf
    # the relative contribution of background counts within sigma_ecf
    f_bkg = p_bkg * source_area / (4 * box_size**2)
    # the contribution of signal counts within sigma_ecf
    signal = f_signal * n / (f_signal + f_bkg)
    # the contribution of background counts within sigma_ecf
    background = f_bkg * n / (f_signal + f_bkg)
    snr = signal / np.sqrt(background)

    # An estimate of how much the data deviates from the model:
    # the 1-sample Kolmogorov-Smirnov (KS) statistic of the marginal distributions.
    # assuming the marginal distributions are Gaussian. This is clearly not perfect.
    yag_model = lambda xx: signal * normal_prob_1d(
        xx, result.x[0], sigma_x
    ) + background * 1 / (2 * box_size)
    zag_model = lambda yy: signal * normal_prob_1d(
        yy, result.x[1], sigma_y
    ) + background * 1 / (2 * box_size)

    x = np.linspace(-4, 4, 1001)
    dx = np.diff(x)[0]
    x = x[:-1] + dx / 2
    xx = x + result.x[0]
    yy = x + result.x[1]

    cdf = np.cumsum(yag_model(xx))
    cdf /= cdf[-1]

    ks_0 = scipy.stats.kstest(
        events[columns[0]],
        scipy.interpolate.interp1d(
            xx, cdf, bounds_error=False, fill_value=(cdf[0], cdf[-1])
        ),
    )

    cdf = np.cumsum(zag_model(yy))
    cdf /= cdf[-1]

    ks_1 = scipy.stats.kstest(
        events[columns[1]],
        scipy.interpolate.interp1d(
            yy, cdf, bounds_error=False, fill_value=(cdf[0], cdf[-1])
        ),
    )

    res = {
        "params": result.x,
        "hess_inv": result.hess_inv.todense(),
        "ndof": len(events) - len(result.x),
        "fit_ok": result.success,
        "p_signal": p_signal,
        columns[0]: result.x[0],
        columns[1]: result.x[1],
        "sigma": (sigma_1, sigma_2),
        "rot_angle": result.x[4],
        f"sigma_{columns[0]}": sigma_x,
        f"sigma_{columns[1]}": sigma_y,
        f"corr_{columns[0]}_{columns[1]}": cov[0, 1] / (sigma_x * sigma_y),
        "psf_ratio": psf_ratio,
        "source_area": source_area,
        "n": n,
        "signal": signal,
        "background": background,
        "snr": snr,
        f"ks_{columns[0]}": ks_0.statistic,
        f"ks_{columns[1]}": ks_1.statistic,
        f"ks_p_value_{columns[0]}": ks_0.pvalue,
        f"ks_p_value_{columns[1]}": ks_1.pvalue,
        # casting to int to avoid issues with serialization as bool
        f"ks_sign_{columns[0]}": int(ks_0.statistic_sign),
        f"ks_sign_{columns[1]}": int(ks_1.statistic_sign),
        "COMPONENT": source["COMPONENT"],
    }
    return res


def fit_gaussians(obs, sources, columns=("y_angle", "z_angle"), box_size=4):
    events = obs.periscope_drift.get_events()
    results = []
    for source in sources:
        sel = (np.abs(events[columns[0]] - source[columns[0]]) < box_size) & (
            np.abs(events[columns[1]] - source[columns[1]]) < box_size
        )
        results.append(
            fit_gaussian_2d(events[sel], source, columns=columns, box_size=box_size)
        )
    return results
