#!/bin/bash

# use a heuristic to generate starting parameters
test -d results_auto && rm -r results_auto
cc_fit.py -f frequencies.dat -d data.dat\
	-c 1 --plot \
	-o results_auto

# supply starting values
test -d results_user && rm -r results_user
cc_fit.py -f frequencies.dat -d data.dat\
	-c 1 --plot \
	--init_file starting_values.dat\
	-o results_user
