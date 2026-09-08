#!/bin/bash
# Everything in the finding, from the two figure PNGs, in order.
set -e
python3 digitize_fig3.py       # markers -> fig3_points.json  (writes the calibration checks)
python3 fig4_coverage.py       # which (alpha,G) conditions appear in the Fig. 4 collapse
python3 structure_test.py      # product form vs additive, whole range vs the Fig-4 window
python3 time_vs_cycles.py      # cycles or seconds; and the alpha>=0.9 arithmetic
python3 fit_scaling.py         # two independent routes to the exponent nu
python3 make_figure.py         # collapse_window.png
