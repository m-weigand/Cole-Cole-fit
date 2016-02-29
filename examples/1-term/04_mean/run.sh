#!/bin/bash

test -d results && rm -r results
cc_fit.py -f fpi/extract_frequencies.dat -d fpi/rho_model_05_specs_fpi.dat -o results --mean -c 1 -p
