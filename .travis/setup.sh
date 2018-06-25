#!/bin/bash
set -ev

{ # try
  julia -e "Pkg.checkout(\"BasisFunctions\",\"$TRAVIS_BRANCH\")"
} || { # catch
  julia -e "Pkg.checkout(\"BasisFunctions\",\"development\")"
}

{
  julia -e "Pkg.checkout(\"Domains\",\"$TRAVIS_BRANCH\")"
} || { # catch
  julia -e "Pkg.checkout(\"Domains\",\"development\")"
}
