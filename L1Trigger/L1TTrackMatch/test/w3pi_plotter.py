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

# Constants
W_PDGID = 24
PI_PDGID = 211 


def w3pi_decay(pdgid_list, data_list):
    n_tps = len(data_list)
    w3pi_data_list = []
    print(f"n tps = {n_tps}")

    for i in range(n_tps):
        if (i == 10000): break
        if (abs(pdgid_list[i]) == 211):
            w3pi_data_list.append(data_list[i])

    return w3pi_data_list
        

def plot_branch(
    branch_name, branch_data, output_file,
    x_range=None, y_range=None, nbins=-1, vline=None, hline=None,
    xticks=None, yticks=None,          # explicit tick positions (lists/arrays)
    xstep=None, ystep=None,            # fixed step between ticks (numbers)
    x_integer=False, y_integer=False,  # force integer ticks
    xfmt=None, yfmt=None,              # formatter callables or str formats
    logy=False
):
    """Plots the given branch data and saves it to a PDF."""
    fig, ax = plt.subplots(figsize=(10, 6))

    bins = int(np.max(branch_data)) if nbins < 0 else nbins
    ax.hist(branch_data, bins=bins, histtype="step", label="Signal")

    ax.set_xlabel(branch_name)
    if x_range: ax.set_xlim(x_range)
    else:       ax.set_xlim(np.min(branch_data), np.max(branch_data))
    if y_range: ax.set_ylim(y_range)
    if vline is not None: ax.axvline(vline, ls='--', lw=0.5, color='black', label=f'x = {float(vline):.2f}')
    if hline is not None: ax.axhline(hline, ls='--', lw=0.5, color='black', label=f'y = {float(hline):.2f}')

    # ---- NEW tick controls ----
    # explicit ticks override everything
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

    ax.legend()
    fig.savefig(output_file, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot for {branch_name} saved to {output_file}")


def main(root_file_path):
    # Define the directories for saving plots
    out_dir = "plots/nstub"

    # Ensure the plots directories exist
    os.makedirs(out_dir, exist_ok=True)

    # Open the ROOT file and navigate to the TTree
    with uproot.open(root_file_path) as file:
        # Navigate to the TDirectoryFile and access the TTree
        directory = file["L1TrackNtuple"]
        tree = directory["w3piTree"]

        # Tracking particle pranches
        tp_pt = tree["tp_pt"].array(library="np")
        """
        tp_eta = tree["tp_eta"].array(library="np")
        tp_phi = tree["tp_phi"].array(library="np")
        tp_dxy = tree["tp_dxy"].array(library="np")
        tp_d0 = tree["tp_d0"].array(library="np")
        tp_z0 = tree["tp_z0"].array(library="np")
        tp_d0_prod = tree["tp_d0_prod"].array(library="np")
        #tp_z0_prod = tree["tp_z0_prod"].array(library="np")
        """
        tp_pdgid = tree["tp_pdgid"].array(library="np")
        tp_mother_pdgid = tree["tp_mother_pdgid"].array(library="np")
        tp_nstub = tree["tp_nstub"].array(library="np")
        """
        tp_nstublayer = tree["tp_nstublayer"].array(library="np")
        tp_evenid = tree["tp_eventid"].array(library="np")
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
        """

        # Extract N_events from data
        N_EVENTS = len(tp_pt)

        """
          - Structure:
            - branch[i][j], i -- event, j -- variable of item (either L1 track or TP) in event
            - the i in each list is the same event
            - the j in each similar group is a variable of an item (TP, L1, matched TP, matched L1) in some event i
        """
        tp_pdgid_list = np.concatenate(tp_pdgid)
        tp_mother_pdgid_list = np.concatenate(tp_mother_pdgid)
        tp_nstub_list = np.concatenate(tp_nstub)
        
        branch_data_list = w3pi_decay(tp_pdgid_list, tp_nstub_list)

        branch_name = "tp_nstub_w3pi"
        branch_data = np.asarray(branch_data_list)
        branch_data_flat = branch_data #np.concatenate(branch_data)
        output_file = os.path.join(out_dir, f"{branch_name}_dist.pdf")
        x_range = 0,20
        y_range = None
        nbins = -1
        vline = None
        hline = None
        x_ticks = None
        y_ticks = None 
        x_step = None
        y_step = None
        x_integer = True
        y_integer = False
        xfmt = None
        yfmt = None
        logy = True
        plot_branch(branch_name, branch_data_flat, output_file, x_range, y_range, nbins, vline, hline, 
                                    x_ticks, y_ticks, x_step, y_step, x_integer, y_integer, xfmt, yfmt, logy)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("root_file", help="Path to the ROOT file.")
    args = parser.parse_args()

    main(args.root_file)