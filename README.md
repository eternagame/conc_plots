# NUPACK-based Concentration Plots

## Setup

* Download and compile Nupack 3 and set the NUPACKHOME environment variable to the nupack root directory
* Install python requirements via `pip install -r requirements.txt`

## Usage

You can make a set of plots for a single sequence by running

```python plot_conc_space.py [SEQUENCE] R105_states.txt -o output.png -t [TITLE]```

The `R105_states.txt` specifies the input and reporter sequences as well as which complexes to include in the calculation.
