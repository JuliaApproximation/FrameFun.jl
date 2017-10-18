#!/bin/bash
set -ev

{ # try
    julia -e "Pkg.checkout(\"BasisFunctions\",\"$TRAVIS_BRANCH\")"
} || { # catch
    julia -e "Pkg.checkout(\"BasisFunctions\",\"master\")"
}
