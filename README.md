# MEMIC model
Model of cancer growth along an oxygen gradient created by a [MEMIC plate](http://carmofon.org/3d-designs/memic/). An interactive web-based version of this model is available [here](https://emilydolson.github.io/memic_model/web/memic_model.html).

[![Build Status](https://travis-ci.com/emilydolson/memic_model.svg?branch=master)](https://travis-ci.com/emilydolson/memic_model) [![codecov](https://codecov.io/gh/emilydolson/memic_model/branch/master/graph/badge.svg)](https://codecov.io/gh/emilydolson/memic_model) [![CodeFactor](https://www.codefactor.io/repository/github/emilydolson/memic_model/badge)](https://www.codefactor.io/repository/github/emilydolson/memic_model) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3633395.svg)](https://doi.org/10.5281/zenodo.3633394)

<img src="https://raw.githubusercontent.com/emilydolson/memic_model/master/normal_memic_canvas.png" width="300">

## Dependencies

-   [Empirical](https://github.com/emilydolson/Empirical/tree/memic_model): Specifically the memic_model branch of my fork (will integrate into upstream master eventually). This is the core framework the model is built on.
-   [Bootstrap](https://getbootstrap.com/): Web framework that the site is built on (takes care of responsiveness for compatability with various devices)
-   [Bootstrap-select](https://developer.snapappointments.com/bootstrap-select/): For styling the drop-down menu.

## Running the model

This model is written in C++, using the Empirical framework to facilitate compiling it either natively (for speed and efficient use on a computing cluster) or to Javascript (to create an interactive web version). The easiest way to run it is just to go to the website for this repository and run the web version. To compile it yourself, though, you can follow these instructions:

First, clone this repository by executing the following command in a terminal:

```bash
git clone https://github.com/emilydolson/memic_model.git
```

This model has a dependency on the Empirical library. To run it, clone that as well:

```bash
git clone https://github.com/emilydolson/Empirical.git
cd Empirical
git checkout memic_model # We're using a few features that haven't been merged into master yet
cd ..
```

Then go into the directory for this repository:

```bash
cd memic_model
```

### Fast version

To run the fast version of this model, compile it and then run it on the command-line:

```bash
make  # compile the code
./memic_model  # run the code
```

This model has a few parameters:
- AGE_LIMIT:                     Age over which non-stem cells die (type=int; default=100)
- ASYMMETRIC_DIVISION_PROB:      Probability of a change in stemness (type=double; default=0)
- BASAL_OXYGEN_CONSUMPTION:      Base oxygen consumption rate (type=double; default=.00075)
- CELL_DIAMETER:                 Cell length and width in microns (type=double; default=20.0)
- DATA_RESOLUTION:               How many updates between printing data? (type=int; default=10)
- DIFFUSION_STEPS_PER_TIME_STEP: Rate at which diffusion is calculated relative to rest of model (type=int; default=10)
- DOSES:                         Number of doses of radiation to apply (type=int; default=0)
- DOSE_SIZE:                     Size of radiation dose to apply in Gy (type=double; default=2.0)
- DOSE_TIME:                     Time point at which to apply radiation (-1 means never) (type=int; default=-1)
- HYPOXIA_DEATH_PROB:            Probability of dieing, given hypoxic conditions (type=double; default=.25)
- INITIAL_OXYGEN_LEVEL:          Initial oxygen level (will be placed in all cells) (type=double; default=.5)
- INIT_POP_SIZE:                 Number of cells to seed population with (type=int; default=100)
- KM:                            Michaelis-Menten kinetic parameter (type=double; default=0.01)
- K_OER:                         Effective OER constant (type=double; default=3.28)
- MITOSIS_PROB:                  Probability of mitosis (type=double; default=.5)
- NEUTRAL_MUTATION_RATE:         Probability of a neutral mutation (only relevant for phylogenetic signature) (type=double; default=.05)
- OER_ALPHA_MAX:                 OER alpha max constant (type=double; default=1.75)
- OER_BETA_MAX:                  OER beta max constant (type=double; default=3.25)
- OER_MIN:                       OER min constant (type=double; default=1)
- OXYGEN_CONSUMPTION_DIVISION:   Amount of oxygen a cell consumes on division (type=double; default=.00075*5)
- OXYGEN_DIFFUSION_COEFFICIENT:  Oxygen diffusion coefficient (type=double; default=.1)
- OXYGEN_THRESHOLD:              How much oxygen do cells need to survive? (type=double; default=.1)
- PLATE_DEPTH:                   Depth of plate in mm (type=double; default=1.45)
- PLATE_LENGTH:                  Length of plate in mm (type=double; default=10.0)
- PLATE_WIDTH:                   Width of plate in mm (type=double; default=6.0)
- SEED:                          Random number generator seed (type=int; default=-1)
- TIME_STEPS:                    Number of time steps to run for (type=int; default=1000)

The values of these can be set with command line flags by placing a dash before the name of the parameter you would like to modify and following it with the desired parameter value:

```bash
# For example, the following runs the model with a NEUTRAL_MUTATION_RATE of .01 and TIME_STEPS set to 100
./memic_model -NEUTRAL_MUTATION_RATE .01 -TIME_STEPS 100
```

### Web version

To compile the web version, you need the [Emscripten C++ to Javascript compiler](https://emscripten.org/). Once you have it installed, you can simply run:

```bash
make web
```

This will compile the web version. To run it, open `web/memic_model.html` in your favorite web browser (it has been tested in Firefox and Chrome).

### Debug version

For additional warning messages, vector bounds checking, memory management checking, and more debugging benefits, you can compile a debug version of this model.

To compile a debugging version of the fast version, run: 

```bash
make debug
```

To compile a debugging version of the web version, run:

```bash
make debug-web
```
