#!/usr/bin/env bash
shopt -s globstar
jupyter nbconvert **/[^_^.]*.ipynb --to html
