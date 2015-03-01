#!/bin/bash

test -d results && rm -r results
cc_fit.py -f frequencies.dat -d data.dat\
	-c 2 --plot
