#!/bin/bash
set -ev

{ # try
    julia -e 'Pkg.checkout("BasisFunctions",$TRAVIS_BRANCH)'
} || { # catch
    julia -e 'Pkg.checkout("BasisFunctions","master")'
}
# if [[ $TRAVIS_BRANCH == 'development' ]]; then
#   julia -e 'Pkg.checkout("BasisFunctions","development")'
# fi
# if [[ $TRAVIS_BRANCH == 'functionset' ]]; then
#   julia -e 'Pkg.checkout("BasisFunctions","functionset")'
# fi
