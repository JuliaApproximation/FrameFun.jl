#!/bin/bash


jupyter nbconvert --to script examples/*.ipynb
mkdir test/scripts/
mv examples/*.jl test/scripts
cd test/
find scripts/*.jl > ../notebookscripts
cd ..
