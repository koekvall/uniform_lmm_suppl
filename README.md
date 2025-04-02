# Supplementary material

The R package used to compute test-statistics can be installed by
`devtools::install_github("koekvall/limestest@R1")`.

Please note the R package, though capable of computing the quantities needed for
the paper "Uniform inference in linear mixed models", is in development
and not yet user friendly.

The simulations for the paper were run on a computing cluster. The file
`lmm-crit-sims.sh` starts the simulations. The file `simulation_master.R` defines
settings and distributes simulations. The file `sim_funs.R` defines functions for
running one replication for the different settings considered.

The simulations and figure for the Introduction were done separately
using the file `sim_fig_intro.R`.

Figures can be created using the file `sec5_figures.R`.
