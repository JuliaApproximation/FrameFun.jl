#!/bin/bash
set -ev
if [[ $TRAVIS_BRANCH == 'development' ]]; then
  julia -e 'Pkg.checkout("BasisFunctions","development")'
fi
