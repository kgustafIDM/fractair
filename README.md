# fractair
Fractional calculus and commercial air transport models used in: arxiv.org/abs/1601.07655

Note that all code is open source and any modifications or usage must also be open source (see GPL v3).

This code base, written in Matlab, is intended to (at least) reproduce the results in the aforementioned manuscript. It contains three folders for organization:
1) air_sim
2) frac_sim
3) viz_util

The first folder should reproduce the simulations of air travel, it contains four files:
runairparamscan.m, which is a wrapper for running air_paramscan23.m; airflight_input.m contains input parameters and loads input files that are constant for the entire parameter scan; airflight_diffusion.m contains the core machinery for solving the transport.

The second folder should reproduce the simulations of fractional diffusion, it contains five files:
runfracparamscan.m, which is a wrapper for frac_paramscan.m; fracflight_input.m holds parameters that are constant during the parameter scan; fracflight_diffusion.m, which has the machinery for solving the fractional diffusion equation using frac_kernel.m

The third folder contains data for the air travel networks, including some pre-processing methods for that data. Visualizations are also possible with the post-processing methods here. This folder also has the propensity.m file with disease dynamics used in the Gillespie simulations. A Google Flu spreadsheet for state-by-state historical ILI trends is also included. Of course this data is quite rough but sufficient for our purposes.
