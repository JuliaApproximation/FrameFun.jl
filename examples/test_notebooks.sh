#!/bin/bash


jupyter nbconvert --to script /examples/*.ipynb
ls /examples/*.jl > notebookscripts
