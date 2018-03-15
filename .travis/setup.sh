#!/bin/bash
set -ev

{ # try
  julia -e "Pkg.checkout(\"BasisFunctions\",\"$TRAVIS_BRANCH\")"
} || { # catch
  julia -e "Pkg.checkout(\"BasisFunctions\",\"master\")"
}

{
  julia -e "Pkg.checkout(\"Domains\",\"$TRAVIS_BRANCH\")"
} || { # catch
  julia -e "Pkg.checkout(\"Domains\",\"master\")"
}
