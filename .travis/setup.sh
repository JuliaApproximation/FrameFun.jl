#!/bin/bash
set -ev
if [[ $TRAVIS_BRANCH == 'development' ]]; then
  julia -e 'Pkg.checkout("BasisFunctions","development")'
fi

if [[ $TRAVIS_BRANCH == 'residual' ]]; then
  julia -e 'Pkg.checkout("BasisFunctions","residual")'
fi
