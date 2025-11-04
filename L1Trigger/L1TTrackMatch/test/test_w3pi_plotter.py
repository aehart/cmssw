#!/usr/bin/env python3
import uproot
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator, MaxNLocator, FuncFormatter
import numpy as np
import warnings
import argparse
import os

# Suppress numpy warnings about subnormal floats
warnings.filterwarnings("ignore", category=UserWarning, module="numpy.core.getlimits")

# Constants (kept for context/extensibility)
W_PDGID = 24
PI_PDGID = 211


# ----------------------------- Plotting helpers ----------------------------- #
def is_integer_like(arr: np.ndarray) -> bool:
    """Return True if data are integers/bools or float values very close to integers."""
    if np.issubdtype(arr.dtype, np.integer) or np.issubdtype(arr.dtype, np.bool_):
        return True
    if np.issubdtype(arr.dtype, np.floating):
        if arr.size == 0:
            return False
        frac = np.abs(arr - np.round(arr))
        return np.nanmax(frac) < 1e-9
    return False


def safe_flatten(jagged):
    """
    Flatten jagged uproot arrays (lists of arrays) to 1D np.ndarray.
    Handles object-dtype arrays and empty events gracefully.
    """
    if isinstance(jagged, np.ndarray):
        if jagged.dtype == object:
            parts = [np.asarray(x) for x in jagged if x is not None]
            if not parts:
                return np.asarray([], dtype=float)
            return np.concatenate(parts) if len(parts) > 1 else parts[0]
        return jagged
    # Fallback for python lists/iterables of arrays
    try:
        return np.concatenate(jagged)
    except Exception:
        return np.asarray(jagged)


def auto_bins(data, nbins=-1, integer_like=False, max_bins=200):
    """
    Choose bins automatically:
      - if nbins >= 0: use it
      - integer-like: one bin per integer across range (capped)
      - else: Freedman–Diaconis rule (capped), fallback to sqrt(N)
    """
    data = data[~np.isnan(data) & ~np.isinf(data)]
    if data.size == 0:
        return 10
    if nbins >= 0:
        return int(nbins)

    dmin, dmax = np.min(data), np.max(data)
    if dmin == dmax:
        return 1

    if integer_like:
        rng = int(np.floor(dmax) - np.ceil(dmin) + 1)
        return min(max(rng, 1), max_bins)

    # Freedman–Diaconis
    q25, q75 = np.percentile(data, [25, 75])
    iqr = q75 - q25
    if iqr <= 0:
        return min(int(np.sqrt(data.size)), max_bins)
    bin_width = 2 * iqr / (data.size ** (1 / 3))
    nb = int(np.ceil((dmax - dmin) / bin_width))
    return min(max(nb, 5), max_bins)


def plot_branch(
    branch_name,
    branch_data,
    output_file,
    x_range=None,
    y_range=None,
    nbins=-1,
    vline=None,
    hline=None,
    xticks=None,
    yticks=None,          # explicit tick positions (lists/arrays)
    xstep=None,
    ystep=None,           # fixed step between ticks (numbers)
    x_integer=False,
    y_integer=False,      # force integer ticks
    xfmt=None,
    yfmt=None,            # formatter callables or str formats
    logy=False,
):
    """Plots the given branch data and saves it to a PDF."""
    fig, ax = plt.subplots(figsize=(10, 6))

    # Determine bins
    bins = auto_bins(branch_data, nbins=nbins, integer_like=x_integer)

    ax.hist(branch_data, bins=bins, histtype="step", label="Entries")

    ax.set_xlabel(branch_name)

    # x/y limits
    if x_range:
        ax.set_xlim(x_range)
    else:
        if branch_data.size > 0:
            xmin, xmax = float(np.min(branch_data)), float(np.max(branch_data))
            if xmin == xmax:
                xmin -= 0.5
                xmax += 0.5
            pad = 0.02 * (xmax - xmin) if xmax > xmin else 1.0
            ax.set_xlim(xmin - pad, xmax + pad)
    if y_range:
        ax.set_ylim(y_range)

    # guide lines
    if vline is not None:
        ax.axvline(vline, ls="--", lw=0.8, color="black", label=f"x = {float(vline):.3g}")
    if hline is not None:
        ax.axhline(hline, ls="--", lw=0.8, color="black", label=f"y = {float(hline):.3g}")

    # tick controls (explicit > step > integer)
    if xticks is not None:
        ax.set_xticks(xticks)
    elif xstep is not None:
        ax.xaxis.set_major_locator(MultipleLocator(xstep))
    elif x_integer:
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    if yticks is not None:
        ax.set_yticks(yticks)
    elif ystep is not None:
        ax.yaxis.set_major_locator(MultipleLocator(ystep))
    elif y_integer:
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    # optional formatters
    if xfmt is not None:
        if isinstance(xfmt, str):
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: xfmt.format(v)))
        else:
            ax.xaxis.set_major_formatter(FuncFormatter(xfmt))
    if yfmt is not None:
        if isinstance(yfmt, str):
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: yfmt.format(v)))
        else:
            ax.yaxis.set_major_formatter(FuncFormatter(yfmt))

    if logy:
        ax.set_yscale("log")

    ax.set_ylabel("Counts")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_file, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot for {branch_name} saved to {output_file}")


# ------------------------------ Script main --------------------------------- #
def main(root_file_path):
    # Define output root directory for saving plots
    out_root = "plots"

    # Ensure the plots directories exist (group subfolders created later)
    os.makedirs(out_root, exist_ok=True)

    # Open the ROOT file and navigate to the TTree
    with uproot.open(root_file_path) as file:
        directory = file["L1TrackNtuple"]
        tree = directory["w3piTree"]

        # ------------------------- Load branches ------------------------- #
        # Tracking particle branches
        tp_pt = tree["tp_pt"].array(library="np")
        tp_eta = tree["tp_eta"].array(library="np")
        tp_phi = tree["tp_phi"].array(library="np")
        tp_dxy = tree["tp_dxy"].array(library="np")
        tp_d0 = tree["tp_d0"].array(library="np")
        tp_z0 = tree["tp_z0"].array(library="np")
        tp_d0_prod = tree["tp_d0_prod"].array(library="np")
        tp_z0_prod = tree["tp_z0_prod"].array(library="np")
        tp_pdgid = tree["tp_pdgid"].array(library="np")
        tp_mother_pdgid = tree["tp_mother_pdgid"].array(library="np")
        tp_nstub = tree["tp_nstub"].array(library="np")
        tp_nstublayer = tree["tp_nstublayer"].array(library="np")
        tp_eventid = tree["tp_eventid"].array(library="np")
        tp_charge = tree["tp_charge"].array(library="np")
        tp_nmatch = tree["tp_nmatch"].array(library="np")
        tp_isHard = tree["tp_isHard"].array(library="np")
        tp_vx = tree["tp_vx"].array(library="np")
        tp_vy = tree["tp_vy"].array(library="np")
        tp_vz = tree["tp_vz"].array(library="np")
        tp_dvx = tree["tp_dvx"].array(library="np")
        tp_dvy = tree["tp_dvy"].array(library="np")
        tp_dvz = tree["tp_dvz"].array(library="np")

        # TP matched to L1 track branches
        tp_matchtrk_pt = tree["tp_matchtrk_pt"].array(library="np")
        tp_matchtrk_eta = tree["tp_matchtrk_eta"].array(library="np")
        tp_matchtrk_phi = tree["tp_matchtrk_phi"].array(library="np")
        tp_matchtrk_d0 = tree["tp_matchtrk_d0"].array(library="np")
        tp_matchtrk_chi2 = tree["tp_matchtrk_chi2"].array(library="np")
        tp_matchtrk_chi2dof = tree["tp_matchtrk_chi2dof"].array(library="np")
        tp_matchtrk_chi2rphi = tree["tp_matchtrk_chi2rphi"].array(library="np")
        tp_matchtrk_chi2rz = tree["tp_matchtrk_chi2rz"].array(library="np")
        tp_matchtrk_bendchi2 = tree["tp_matchtrk_bendchi2"].array(library="np")
        tp_matchtrk_nstub = tree["tp_matchtrk_nstub"].array(library="np")
        tp_matchtrk_nstublayer = tree["tp_matchtrk_nstublayer"].array(library="np")
        tp_matchtrk_charge = tree["tp_matchtrk_charge"].array(library="np")

        # L1 track branches
        trk_pt = tree["trk_pt"].array(library="np")
        trk_eta = tree["trk_eta"].array(library="np")
        trk_d0 = tree["trk_d0"].array(library="np")
        trk_chi2 = tree["trk_chi2"].array(library="np")
        trk_chi2dof = tree["trk_chi2dof"].array(library="np")
        trk_chi2rphi = tree["trk_chi2rphi"].array(library="np")
        trk_chi2rz = tree["trk_chi2rz"].array(library="np")
        trk_bendchi2 = tree["trk_bendchi2"].array(library="np")
        trk_nstub = tree["trk_nstub"].array(library="np")
        trk_nstublayer = tree["trk_nstublayer"].array(library="np")

        # L1 track matched to TP branches
        trk_fake = tree["trk_fake"].array(library="np")
        trk_matchtp_pt = tree["trk_matchtp_pt"].array(library="np")
        trk_matchtp_eta = tree["trk_matchtp_eta"].array(library="np")
        trk_matchtp_pdgid = tree["trk_matchtp_pdgid"].array(library="np")
        trk_matchtp_mother_pdgid = tree["trk_matchtp_mother_pdgid"].array(library="np")
        trk_matchtp_nstub = tree["trk_matchtp_nstub"].array(library="np")
        trk_matchtp_nstublayer = tree["trk_matchtp_nstublayer"].array(library="np")
        trk_matchtp_isHard = tree["trk_matchtp_isHard"].array(library="np")
        trk_matchtp_vx = tree["trk_matchtp_vx"].array(library="np")
        trk_matchtp_vy = tree["trk_matchtp_vy"].array(library="np")
        trk_matchtp_vz = tree["trk_matchtp_vz"].array(library="np")
        trk_matchtp_dvx = tree["trk_matchtp_dvx"].array(library="np")
        trk_matchtp_dvy = tree["trk_matchtp_dvy"].array(library="np")
        trk_matchtp_dvz = tree["trk_matchtp_dvz"].array(library="np")

        # Extract number of events (example; not used elsewhere)
        N_EVENTS = len(tp_pt)
        print(f"[info] Loaded {N_EVENTS} events.")

        # ------------------------- Plot many branches ------------------------- #
        # group -> [branch names]
        groups = {
            "tp": [
                "tp_pt", "tp_eta", "tp_phi", "tp_dxy", "tp_d0", "tp_z0",
                "tp_d0_prod", "tp_z0_prod", "tp_pdgid", "tp_mother_pdgid",
                "tp_nstub", "tp_nstublayer", "tp_eventid", "tp_charge",
                "tp_nmatch", "tp_isHard", "tp_vx", "tp_vy", "tp_vz",
                "tp_dvx", "tp_dvy", "tp_dvz",
            ],
            "tp_matchtrk": [
                "tp_matchtrk_pt", "tp_matchtrk_eta", "tp_matchtrk_phi",
                "tp_matchtrk_d0", "tp_matchtrk_chi2", "tp_matchtrk_chi2dof",
                "tp_matchtrk_chi2rphi", "tp_matchtrk_chi2rz", "tp_matchtrk_bendchi2",
                "tp_matchtrk_nstub", "tp_matchtrk_nstublayer", "tp_matchtrk_charge",
            ],
            "trk": [
                "trk_pt", "trk_eta", "trk_d0", "trk_chi2", "trk_chi2dof",
                "trk_chi2rphi", "trk_chi2rz", "trk_bendchi2",
                "trk_nstub", "trk_nstublayer",
            ],
            "trk_matchtp": [
                "trk_fake", "trk_matchtp_pt", "trk_matchtp_eta",
                "trk_matchtp_pdgid", "trk_matchtp_mother_pdgid",
                "trk_matchtp_nstub", "trk_matchtp_nstublayer", "trk_matchtp_isHard",
                "trk_matchtp_vx", "trk_matchtp_vy", "trk_matchtp_vz",
                "trk_matchtp_dvx", "trk_matchtp_dvy", "trk_matchtp_dvz",
            ],
        }

        # Per-branch overrides (optional): x-range, logy, integer ticks, etc.
        overrides = {
            "tp_nstub":       dict(x_range=(0, 20), x_integer=True, logy=True),
            "tp_nstublayer":  dict(x_range=(0, 10), x_integer=True, logy=True),
            "tp_charge":      dict(x_range=(-2, 2), x_integer=True, logy=False),
            "tp_isHard":      dict(x_range=(-0.5, 1.5), x_integer=True, logy=False),
            "trk_nstub":      dict(x_range=(0, 20), x_integer=True, logy=True),
            "trk_nstublayer": dict(x_range=(0, 10), x_integer=True, logy=True),
            "trk_fake":       dict(x_range=(-0.5, 1.5), x_integer=True, logy=False),
        }

        # Create group folders
        for group in groups:
            os.makedirs(os.path.join(out_root, group), exist_ok=True)

        # Use locals() map to fetch arrays by name
        loc = locals()

        for group, names in groups.items():
            for name in names:
                if name not in loc:
                    print(f"[warn] branch '{name}' not loaded, skipping.")
                    continue

                raw = loc[name]
                data = safe_flatten(raw).astype(float, copy=False)
                # Clean NaN/Inf
                data = data[~np.isnan(data) & ~np.isinf(data)]
                if data.size == 0:
                    print(f"[warn] branch '{name}' is empty after cleaning, skipping.")
                    continue

                # Decide integer-like & bins
                is_int = is_integer_like(data)
                nb = auto_bins(data, nbins=-1, integer_like=is_int)

                # Defaults
                kw = dict(
                    x_range=None, y_range=None, nbins=nb,
                    vline=None, hline=None,
                    xticks=None, yticks=None,
                    xstep=None, ystep=None,
                    x_integer=is_int, y_integer=False,
                    xfmt=None, yfmt=None,
                    logy=True,  # default to log-y for count histograms
                )

                # Apply overrides if present
                if name in overrides:
                    kw.update(overrides[name])

                # If no x_range, build from data with a small pad
                if kw["x_range"] is None:
                    xmin, xmax = np.min(data), np.max(data)
                    if xmin == xmax:
                        xmin -= 0.5
                        xmax += 0.5
                    pad = 0.02 * (xmax - xmin) if xmax > xmin else 1.0
                    kw["x_range"] = (xmin - pad, xmax + pad)

                out_file = os.path.join(out_root, group, f"{name}_dist.pdf")
                print(f"Plotting {group}/{name} -> {out_file} (bins={kw['nbins']}, int_like={kw['x_integer']})")
                plot_branch(name, data, out_file, **kw)


# -------------------------------- Entrypoint -------------------------------- #
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("root_file", help="Path to the ROOT file.")
    args = parser.parse_args()
    main(args.root_file)
