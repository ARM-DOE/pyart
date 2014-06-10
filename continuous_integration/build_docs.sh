#!/bin/bash

cd "$TRAVIS_BUILD_DIR"

echo "Building Docs"
conda install --yes sphinx

mv "$TRAVIS_BUILD_DIR"/doc /tmp
cd /tmp/doc
make html

exit 0
