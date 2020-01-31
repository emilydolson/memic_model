# MEMIC model
Model of cancer growth under spatially heterogeneous oxygen

[![Build Status](https://travis-ci.com/emilydolson/memic_model.svg?branch=master)](https://travis-ci.com/emilydolson/memic_model) [![codecov](https://codecov.io/gh/emilydolson/memic_model/branch/master/graph/badge.svg)](https://codecov.io/gh/emilydolson/memic_model) [![CodeFactor](https://www.codefactor.io/repository/github/emilydolson/memic_model/badge)](https://www.codefactor.io/repository/github/emilydolson/memic_model)[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3633395.svg)](https://doi.org/10.5281/zenodo.3633395)



## Dependencies

-   [Empirical](https://github.com/emilydolson/Empirical/tree/memic_model): Specifically the memic_model branch of my fork (will integrate into upstream master eventually). This is the core framework the model is built on.
-   [Bootstrap](https://getbootstrap.com/): Web framework that the site is built on (takes care of responsiveness for compatability with various devices)
-   [Bootstrap-select](https://developer.snapappointments.com/bootstrap-select/): For styling the drop-down menu.

## Tentative specs

*Definitely necessary*
-   Population of cells
-   Cells can be in different states
-   Cells can reproduce (sometimes with mutations)
-   Environmental conditions vary over space
-   Environmental conditions vary over time (based on diffusion equations?)
-   Environmental conditions can affect cell reproduction, death, and behavior
-   Ability to eventually simulate applying treatments

*Possibly necessary*
-   Complex topology: 3D, free-moving cells (i.e. continuous x,y coordinates rather than lattice), something else?
-   Variation in cell proliferation ability
-   Cell movement (preferential based on environment?)
-   Non-neutral mutations (and possibly mutations with differential effects based on environment)
-   Phenotypic plasticity
-   Physics
