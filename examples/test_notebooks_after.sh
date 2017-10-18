#!/bin/bash

folder=`pwd`
echo "Remove $folder/examples/*.jl"
rm $folder/examples/*.jl
echo "Remove $folder/notebookscripts"
rm $folder/notebookscripts
