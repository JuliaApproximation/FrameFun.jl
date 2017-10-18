#!/bin/bash


jupyter nbconvert --to script examples/*.ipynb
find examples/*.jl > notebookscripts
