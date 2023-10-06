#!/usr/bin/env python

# import re
import argparse
import functools
import logging
import os
from pathlib import Path

import jinja2
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from cxotime import CxoTime
from cxotime import units as u
from scipy.interpolate import interp1d

# from astromon import cross_match
from Ska.Matplotlib import cxctime2plotdate, plot_cxctime

from astromon import db, utils

plt.rcParams["font.size"] = "16"
plt.rcParams["figure.max_open_warning"] = 100

JINJA2 = jinja2.Environment(
    loader=jinja2.PackageLoader("astromon.web", "templates"),
    autoescape=jinja2.select_autoescape(["html", "xml"]),
)


def plot_offsets_history(
    matches,
    title="Offsets History",
    filename="offsets-history.png",
    dy_median=None,
    dz_median=None,
):
    fig, ax = plt.subplots(1, 1, figsize=(12, 4), sharex=True)
    plt.sca(ax)
    plot_cxctime(
        matches["time"], matches["dy"], ".", label="dy", color="tab:blue", alpha=0.4
    )
    plot_cxctime(
        matches["time"], matches["dz"], ".", label="dz", color="tab:orange", alpha=0.4
    )
    if dy_median is not None:
        plot_cxctime(
            dy_median["time"], dy_median["median"], color="tab:blue", linewidth=3
        )
    if dz_median is not None:
        plot_cxctime(
            dz_median["time"], dz_median["median"], color="tab:orange", linewidth=3
        )
    plt.legend(loc="upper left")
    plt.title(title)
    plt.ylabel("offset (arcsec)")
    plt.ylim((-2, 2))
    plt.grid()
    filename = Path(filename)
    filename.parent.mkdir(exist_ok=True, parents=True)
    plt.savefig(filename)


def plot_offsets_q_history(
    matches,
    dy_median,
    dz_median,
    filename="offsets-q-history.png",
):
    assert dy_median is not None
    assert dz_median is not None

    after = matches["after_caldb"]
    # tmax = np.min(matches["time"][after])

    times = CxoTime([dy_median["start"], dy_median["stop"]], format="frac_year")
    dates = cxctime2plotdate(times)

    # dy_median = dy_median[times[0] < tmax]
    # dz_median = dz_median[times[0] < tmax]
    # times = times[:, times[0] < tmax]

    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    plt.sca(ax0)

    plot_cxctime(
        matches["time"][~after],
        matches["dy"][~after],
        ".",
        label="dy",
        color="k",
        alpha=0.2,
    )
    plot_cxctime(
        matches["time"][after],
        matches["dy"][after],
        ".",
        label="dy",
        color="r",
        alpha=0.5,
    )
    plot_cxctime(
        times, np.tile(dy_median["median"], (2, 1)), "-", linewidth=2, color="tab:blue"
    )
    for i in range(len(dy_median)):
        plt.fill_between(
            dates.T[i],
            [dy_median["sigma_minus"][i]] * 2,
            [dy_median["sigma_plus"][i]] * 2,
            color="tab:blue",
            alpha=0.4,
            linewidth=0,
        )

    # plt.legend(loc='upper left')
    plt.title(f"dY")
    plt.ylabel("offset (arcsec)")
    plt.ylim((-1.1, 1.1))
    plt.grid()

    plt.sca(ax1)
    plot_cxctime(
        matches["time"][~after],
        matches["dz"][~after],
        ".",
        label="dz",
        color="k",
        alpha=0.2,
    )
    plot_cxctime(
        matches["time"][after],
        matches["dz"][after],
        ".",
        label="dz",
        color="r",
        alpha=0.5,
    )
    plot_cxctime(
        times,
        np.tile(dz_median["median"], (2, 1)),
        "-",
        linewidth=2,
        color="tab:orange",
    )
    for i in range(len(dz_median)):
        plt.fill_between(
            dates.T[i],
            [dz_median["sigma_minus"][i]] * 2,
            [dz_median["sigma_plus"][i]] * 2,
            color="tab:orange",
            alpha=0.4,
            linewidth=0,
        )

    # plt.legend(loc='upper left')
    plt.title(f"dZ")
    plt.ylabel("offset (arcsec)")
    plt.ylim((-1.1, 1.1))
    plt.grid()
    if filename:
        filename = Path(filename)
        filename.parent.mkdir(exist_ok=True, parents=True)
        plt.savefig(filename)


def cdf_2_(matches, quantiles=(0.68, 0.90, 0.99)):
    bins, cdf, quantiles = cdf_(matches)
    cdf_err = np.zeros((40, len(bins)))
    for i in range(cdf_err.shape[0]):
        _, cdf_err[i], _ = cdf_(
            matches[np.random.choice(np.arange(len(matches)), len(matches))]
        )
    cdf_err = np.quantile(cdf_err, [0.159, 0.841], axis=0)

    return bins, cdf, quantiles, cdf_err


def cdf_(matches, quantiles=(0.68, 0.90, 0.99)):
    h, bins = np.histogram(matches["dr"], bins=np.linspace(0, 3, 301))
    cdf = np.cumsum(h) / len(matches)
    cdf = np.concatenate([[0.0], cdf])
    f = interp1d(cdf, bins)

    res = [{"q": q, "offset": f(q)} for q in quantiles]
    return bins, cdf, res


def plot_cdf_3(
    all_matches,
    col,
    quantiles=[],
    groupby="year_bin_2",
    title="",
    filename=None,
    xlims=None,
):
    _, ax = plt.subplots(1, 1, figsize=(6, 4))
    plt.sca(ax)

    bins = np.linspace(-2.0, 2.0, 1001)
    vals, bins = np.histogram(all_matches[col], bins)
    vals = vals / np.sum(vals)
    vals = np.cumsum(vals)
    edges = np.zeros((all_matches.meta[f"{groupby}_edges"].shape[0] - 1, 2))
    edges.T[0] = all_matches.meta[f"{groupby}_edges"][:-1]
    edges.T[1] = all_matches.meta[f"{groupby}_edges"][1:]

    cdf_err = np.zeros((100, len(vals)))
    for i in range(cdf_err.shape[0]):
        m = all_matches[np.random.choice(np.arange(len(all_matches)), len(all_matches))]
        cdf_err[i], _ = np.histogram(m[col], bins)
        cdf_err[i] = cdf_err[i] / np.sum(cdf_err[i])
        cdf_err[i] = np.cumsum(cdf_err[i])
    cdf_err = np.quantile(cdf_err, [0.159, 0.841], axis=0)

    plt.ylim((0, 1))
    if xlims:
        plt.xlim(xlims)
    else:
        plt.xlim((bins[0], bins[-1]))

    plt.fill_between(bins[1:], cdf_err[0], cdf_err[1], step="post", alpha=0.3)

    plt.stairs(
        vals,
        bins,
        alpha=0.5,
        linewidth=2,
    )
    # print(quantiles)
    for q in quantiles:
        r = q["offset"]
        q = q["q"]
        if r > plt.xlim()[1]:
            continue
        plt.plot([0.0, r], [q, q], "--", color="b", linewidth=1)
        plt.plot([r, r], [0.0, q], "--", color="b", linewidth=1)
        # plt.xlabel('Radial offset (arcsec)')
        # plt.ylabel('Cumulative fraction')
        plt.text(
            r + 0.02,
            0.1,
            f"{r:.2} arcsec, {q*100:.0f}%",
            rotation="vertical",
            horizontalalignment="left",
        )

    plt.title(title)
    plt.xlabel("Radial offset (arcsec)")
    plt.ylabel("Cumulative fraction")
    plt.grid()
    plt.tight_layout()
    if filename:
        filename = Path(filename)
        filename.parent.mkdir(exist_ok=True, parents=True)
        plt.savefig(filename)


def plot_cdf(
    bins,
    cdf,
    quantiles,
    title="Offset Cumulative Distribution",
    filename="offsets-cdf.png",
    cdf_err=None,
):
    _, ax = plt.subplots(1, 1, figsize=(6, 4))
    plt.sca(ax)

    if cdf_err is not None:
        plt.fill_between(bins, cdf_err[0], cdf_err[1], step="post", alpha=0.3)
    plt.step(bins, cdf, where="post")
    plt.xlim((0, 1.1))
    plt.ylim((0, 1.0))

    for q in quantiles:
        r = q["offset"]
        q = q["q"]
        if r > plt.xlim()[1]:
            continue
        plt.plot([0.0, r], [q, q], "--", color="b", linewidth=1)
        plt.plot([r, r], [0.0, q], "--", color="b", linewidth=1)
        plt.xlabel("Radial offset (arcsec)")
        plt.ylabel("Cumulative fraction")
        plt.text(
            r + 0.02,
            0.1,
            f"{r:.2} arcsec, {q*100:.0f}%",
            rotation="vertical",
            horizontalalignment="left",
        )
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    if filename:
        filename = Path(filename)
        filename.parent.mkdir(exist_ok=True, parents=True)
        plt.savefig(filename)


def year_bins(time, bins_per_year=2):
    """Compute the median value in each bin of the input data."""
    years = np.array(time.frac_year)
    ymin = np.floor(np.min(years / bins_per_year)) * bins_per_year
    ymax = np.ceil(np.max(years / bins_per_year)) * bins_per_year
    nbins = int((ymax - ymin) * bins_per_year)
    return np.linspace(ymin, ymax, nbins + 1)


def binned_median(time, vals, bins_per_year=2):
    """Compute the median value in each bin of the input data."""
    years = np.array(time.frac_year)
    bins = year_bins(time, bins_per_year)
    years_bin = np.digitize(years, bins)
    t = Table([years_bin, years, vals], names=["years_bin", "years", "vals"])
    tg = t.group_by("years_bin")
    tga = tg.groups.aggregate(np.median)
    tga["q95"] = tg[["vals"]].groups.aggregate(functools.partial(np.quantile, q=0.95))[
        "vals"
    ]
    tga["q05"] = tg[["vals"]].groups.aggregate(functools.partial(np.quantile, q=0.05))[
        "vals"
    ]
    tga["q75"] = tg[["vals"]].groups.aggregate(functools.partial(np.quantile, q=0.75))[
        "vals"
    ]
    tga["q25"] = tg[["vals"]].groups.aggregate(functools.partial(np.quantile, q=0.25))[
        "vals"
    ]
    tga["sigma_minus"] = tg[["vals"]].groups.aggregate(
        functools.partial(np.quantile, q=0.158)
    )["vals"]
    tga["sigma_plus"] = tg[["vals"]].groups.aggregate(
        functools.partial(np.quantile, q=0.842)
    )["vals"]
    tga["start"] = bins[tga["years_bin"] - 1]
    tga["stop"] = bins[tga["years_bin"]]
    tga["time"] = CxoTime(tga["years"], format="frac_year")
    tga.rename_columns(["years", "vals"], ["median_year", "median"])
    return tga


def plot_cdf_2(
    all_matches,
    col,
    groupby="year_bin_2",
    cmap="viridis",
    title="",
    filename=None,
    xlims=None,
    loc="upper left",
):
    _, ax = plt.subplots(1, 1, figsize=(6, 4))
    plt.sca(ax)
    g = all_matches.group_by(groupby)
    g.add_index(groupby)
    bmax = np.max(all_matches[groupby])
    imax = 6
    cmap = plt.get_cmap(cmap)
    bins = np.linspace(-2.0, 2.0, 1001)
    edges = np.zeros((all_matches.meta[f"{groupby}_edges"].shape[0] - 1, 2))
    edges.T[0] = all_matches.meta[f"{groupby}_edges"][:-1]
    edges.T[1] = all_matches.meta[f"{groupby}_edges"][1:]

    plt.ylim((0, 1))
    if xlims:
        plt.xlim(xlims)
    else:
        plt.xlim((bins[0], bins[-1]))

    for j, i in enumerate(range(0, imax)):
        b = bmax - i
        if b not in g[groupby]:
            continue
        vals, bins = np.histogram(g.loc[b][col], bins)
        vals = vals / np.sum(vals)
        vals = np.cumsum(vals)
        plt.stairs(
            vals,
            bins,
            color=cmap.colors[int(i * cmap.N // imax)],
            alpha=0.5,
            linewidth=2,
            label=f"{edges[b - 1][0]:.1f}-{edges[b - 1][1]:.1f}",
        )
    plt.title(title)
    plt.xlabel("Radial offset (arcsec)")
    plt.ylabel("Cumulative fraction")
    plt.legend(loc=loc, fontsize="xx-small")
    plt.tight_layout()
    if filename:
        filename = Path(filename)
        filename.parent.mkdir(exist_ok=True, parents=True)
        plt.savefig(filename)


def create_figures_mta(outdir, calalign_dir=None, use_reference_calalign=False):
    outdir = Path(outdir)

    n_years = 5
    snr = 3
    sim_z = 4  # max sim-z
    draw_median = True

    # result = {'snr': snr, 'n_years': n_years}
    end = CxoTime()
    start = end - n_years * u.year
    result = {
        "start_date": start.iso.split()[0],
        "end_date": end.iso.split()[0],
    }

    all_matches = db.get_cross_matches(
        snr=snr,
        exclude_bad_targets=True,
        sim_z=sim_z,
        exclude_categories=[
            "SN, SNR, and Isolated NS",
            "Solar System and Misc",
            "Clusters of Galaxies",
        ],
    )
    no_version = all_matches["caldb_version"] == "0.0"
    if np.any(no_version):
        logging.getLogger("celmon").warning("Some observations with no version")
        all_matches = all_matches[~no_version]

    calalign = utils.get_calalign_offsets(all_matches, calalign_dir=calalign_dir)
    all_matches["after_caldb"] = calalign["after_caldb"]
    tag = "-archive"
    if use_reference_calalign:
        all_matches["dy"] -= calalign["calalign_dy"] - calalign["ref_calalign_dy"]
        all_matches["dz"] -= calalign["calalign_dz"] - calalign["ref_calalign_dz"]
        all_matches["dr"] = np.sqrt(all_matches["dy"] ** 2 + all_matches["dz"] ** 2)
        tag = ""

    all_matches["year"] = all_matches["time"].frac_year
    year_bin_2 = year_bins(all_matches["time"], 2)
    all_matches["year_bin_2"] = np.digitize(all_matches["year"], year_bin_2)
    all_matches.meta["year_bin_2_edges"] = year_bin_2

    matches = all_matches[all_matches["time"] > CxoTime() - n_years * u.year]

    bins_per_year = 2
    dy_median = binned_median(
        all_matches["time"], all_matches["dy"], bins_per_year=bins_per_year
    )
    dz_median = binned_median(
        all_matches["time"], all_matches["dz"], bins_per_year=bins_per_year
    )

    _, _, quantiles = cdf_(matches)
    result.update({f'q{100*q["q"]:.0f}': q["offset"] for q in quantiles})
    result["max_offset"] = np.abs(np.max(matches["dr"]))
    plot_offsets_history(
        all_matches,
        title="Offsets History",
        filename=outdir / f"offsets-history{tag}.png",
        dy_median=dy_median,
        dz_median=dz_median,
    )
    plot_cdf_3(
        matches,
        "dr",
        xlims=(0, 1.1),
        quantiles=quantiles,
        title="Offset Cumulative Distribution",
        filename=outdir / f"offsets-cdf{tag}.png",
    )
    plot_cdf_2(
        matches,
        "dy",
        xlims=(-1.5, 1.5),
        filename=outdir / f"offsets-dy-cdf-multi-year{tag}.png",
        title="dY",
        loc="lower right",
    )
    plot_cdf_2(
        matches,
        "dz",
        xlims=(-1.5, 1.5),
        filename=outdir / f"offsets-dz-cdf-multi-year{tag}.png",
        title="dZ",
        loc="lower right",
    )
    plot_cdf_2(
        matches,
        "dr",
        xlims=(0, 1.5),
        filename=outdir / f"offsets-dr-cdf-multi-year{tag}.png",
        title="dR",
        loc="lower right",
    )
    plot_offsets_q_history(
        matches=all_matches,
        dy_median=dy_median,
        dz_median=dz_median,
        filename=outdir / f"offsets-q-history{tag}.png",
    )

    for det in np.unique(matches["detector"]):
        m = matches[matches["detector"] == det]
        all_m = all_matches[all_matches["detector"] == det]
        _, _, quantiles = cdf_(m)
        result.update(
            {
                det.replace("-", "_"): {
                    f'q{100*q["q"]:.0f}': q["offset"] for q in quantiles
                }
            }
        )
        plot_offsets_history(
            all_m,
            title=det,
            filename=outdir / f"offsets-{det}-history{tag}.png",
            dy_median=dy_median if draw_median else None,
            dz_median=dz_median if draw_median else None,
        )
        plot_cdf_3(
            m,
            "dr",
            xlims=(0, 1.1),
            quantiles=quantiles,
            title=f"{det} ({len(m)} points)",
            filename=outdir / f"offsets-{det}-cdf{tag}.png",
        )
    return result


def create_figures_cal(
    outdir,
    snr=5,
    n_years=5,
    draw_median=True,
    calalign_dir=None,
    use_reference_calalign=False,
):
    outdir = Path(outdir)

    sim_z = 4  # max sim-z

    end = CxoTime()
    start = end - n_years * u.year
    matches = db.get_cross_matches(
        snr=snr,
        exclude_bad_targets=True,
        sim_z=sim_z,
        exclude_categories=[
            "SN, SNR, and Isolated NS",
            "Solar System and Misc",
            "Clusters of Galaxies",
        ],
    )

    no_version = matches["caldb_version"] == "0.0"
    if np.any(no_version):
        logging.getLogger("celmon").warning("Some observations with no version")
        matches = matches[~no_version]

    calalign = utils.get_calalign_offsets(matches, calalign_dir=calalign_dir)
    matches["after_caldb"] = calalign["after_caldb"]
    tag = "-archive"
    if use_reference_calalign:
        matches["dy"] -= calalign["calalign_dy"] - calalign["ref_calalign_dy"]
        matches["dz"] -= calalign["calalign_dz"] - calalign["ref_calalign_dz"]
        matches["dr"] = np.sqrt(matches["dy"] ** 2 + matches["dz"] ** 2)
        tag = ""

    ok = matches["time"] > start
    matches = matches[ok]

    dy_median = binned_median(matches["time"], matches["dy"], bins_per_year=2)
    dz_median = binned_median(matches["time"], matches["dz"], bins_per_year=2)

    result = {
        "snr": snr,
        "n_years": n_years,
        "start_date": start.iso.split()[0],
        "end_date": end.iso.split()[0],
    }

    bins, cdf, quantiles = cdf_(matches)
    result.update({f'q{100*q["q"]:.0f}': q["offset"] for q in quantiles})
    result["max_offset"] = np.abs(np.max(matches["dr"]))
    plot_offsets_history(
        matches,
        title="Offsets History",
        filename=outdir / f"offsets-history{tag}.png",
        dy_median=dy_median,
        dz_median=dz_median,
    )
    plot_cdf(
        bins,
        cdf,
        quantiles,
        title="Offset Cumulative Distribution",
        filename=outdir / f"offsets-cdf{tag}.png",
    )
    for det in np.unique(matches["detector"]):
        result.update(
            {
                det.replace("-", "_"): {
                    f'q{100*q["q"]:.0f}': q["offset"] for q in quantiles
                }
            }
        )
        m = matches[matches["detector"] == det]
        bins, cdf, quantiles = cdf_(m)
        plot_offsets_history(
            m,
            title=det,
            filename=outdir / f"offsets-{det}-history{tag}.png",
            dy_median=dy_median if draw_median else None,
            dz_median=dz_median if draw_median else None,
        )
        for ptype in ["png", "pdf"]:
            plot_cdf(
                bins,
                cdf,
                quantiles,
                title=f"{det} ({len(m)} points)",
                filename=outdir / f"offsets-{det}-cdf{tag}.{ptype}",
            )
    return result


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", "-o", type=Path, default="./celmon")
    parser.add_argument(
        "--db-file",
        type=Path,
        default=Path(os.environ["SKA"]) / "data" / "astromon" / "astromon.h5",
    )
    parser.add_argument(
        "--calalign-dir",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"],
    )
    return parser


def main():
    from ska_helpers.logging import basic_logger

    args = get_parser().parse_args()
    logger = basic_logger(
        "celmon", format="%(levelname)s %(funcName)s %(message)s", level=args.log_level
    )

    if args.db_file:
        os.environ["ASTROMON_FILE"] = str(args.db_file)
        import importlib

        importlib.reload(db)

    (args.out / "cal").mkdir(exist_ok=True, parents=True)
    (args.out / "mta").mkdir(exist_ok=True, parents=True)

    data_cal = create_figures_cal(
        outdir=args.out / "cal", calalign_dir=args.calalign_dir
    )
    data_cal_ref = create_figures_cal(
        outdir=args.out / "cal",
        calalign_dir=args.calalign_dir,
        use_reference_calalign=True,
    )

    data_mta = create_figures_mta(
        outdir=args.out / "mta", calalign_dir=args.calalign_dir
    )
    data_mta_ref = create_figures_mta(
        outdir=args.out / "mta",
        calalign_dir=args.calalign_dir,
        use_reference_calalign=True,
    )

    tpl = JINJA2.get_template("celmon_cal.html")
    file_path = args.out / "cal" / "index.html"
    with open(file_path, "w") as out:
        out.write(tpl.render(data={"cal": data_cal, "cal_ref": data_cal_ref}))
        logger.info(f"report created at {file_path}")

    tpl = JINJA2.get_template("celmon_mta.html")
    file_path = args.out / "mta" / "index.html"

    with open(file_path, "w") as out:
        out.write(
            tpl.render(
                data={
                    "cal": data_cal,
                    "cal_ref": data_cal_ref,
                    "mta": data_mta,
                    "mta_ref": data_mta_ref,
                }
            )
        )
        logger.info(f"report created at {file_path}")


if __name__ == "__main__":
    main()
