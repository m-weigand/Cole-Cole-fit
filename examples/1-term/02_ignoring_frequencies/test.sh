#!/bin/bash

test -d results && rm -r results

# shoud raise an error message
cc_fit.py -f frequencies.dat -d data.dat\
	-c 1 --plot\
	--ignore "0;1;2;3;4;5"

# shoud raise an error message
cc_fit.py -f frequencies.dat -d data.dat\
	-c 1 --plot\
	--ignore "1;2;3;4;5;60"

# should run
cc_fit.py -f frequencies.dat -d data.dat\
	-c 1 --plot\
	--ignore "1;2;3;4;5"
