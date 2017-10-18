#!/bin/bash

folder=`pwd`
jupyter nbconvert --to script $folder/examples/*.ipynb
ls $folder/examples/*.jl > notebookscripts
