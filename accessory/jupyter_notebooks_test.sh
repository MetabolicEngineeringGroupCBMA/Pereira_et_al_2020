#!/usr/bin/env bash

shopt -s globstar

cd notebooks

#command to test all notebooks
pytest -vvv --current-env --capture=fd --nbval **/*.ipynb

exit $?
