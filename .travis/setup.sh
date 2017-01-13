#!/bin/bash
set -ev
if [[ $TRAVIS_BRANCH == 'development' ]]; then
  cd BasisFunctions/
  git checkout development
  cd ..
fi
