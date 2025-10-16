import argparse
import matplotlib.pyplot as plt
from pychart.io_utils import load_data
from pychart.plot import PlotData, add_map_plot, add_cnt_plot
from pychart.figure import FigureBuilder  # <- new class file replacing buid_figure_layout, etc.
from pychart import cb
from pychart.cli import parse_args, sanity_checks_args, get_config
import time

# Add a print plot summary and a print of a YML to read ?
# need to think about it

def main():

    # prepare config dictionary from input arguments
    args = parse_args()
    args = sanity_checks_args(args)
    pychart_config = get_config(args)

    # --- Build the figure using the new class ---
    figure_config = pychart_config["figure"]
    map_config = pychart_config["map"]
    cnt_config = pychart_config["cnt"]
    cb_config = pychart_config["cb"]

    fb = FigureBuilder(config=figure_config, projections_yaml="projections.yml")
    fig, axes = fb.build_layout()                 # builds figure + subplots + titles

    # --- Loop through subplots ---
    for iax, ax in enumerate(axes):
        # MAP
        if iax < len(map_config["files"]):
            map_cb, pcol = add_map_plot(map_config, cb_config, figure_config, iax, ax)

        # CONTOUR
        if iax < len(cnt_config["files"]):
            add_cnt_plot(cnt_config, iax, ax)

        else:
            ax.set_visible(False)

    # --- Add colorbar using the class method ---
    fb.add_colorbar(pcol,map_cb)

    plt.show()
    
    print(f"Saving figure to {figure_config['output']}")
    fig.savefig(figure_config["output"], dpi=150)


if __name__ == "__main__":
    main()