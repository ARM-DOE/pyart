#!/bin/bash
# Adapted from the ci/build_docs.sh file from the pandas project
# https://github.com/pydata/pandas
set -e

cd "$TRAVIS_BUILD_DIR"

echo "Building Docs"
conda install -q sphinx pil

mv "$TRAVIS_BUILD_DIR"/doc /tmp
cd /tmp/doc
mv -f source/index.ci source/index.rst
make html

# upload to pyart-docs-travis repo is this is not a pull request and
# secure token is available (aka in the ARM-DOE repository.
if [ "$TRAVIS_PULL_REQUEST" == "false" ] && [ $TRAVIS_SECURE_ENV_VARS == 'true' ]; then
    cd /tmp/doc/build/html
    git config --global user.email "pyart-docs-bot@example.com"
    git config --global user.name "pyart-docs-bot"

    git init
    touch README
    git add README
    git commit -m "Initial commit" --allow-empty
    git branch gh-pages
    git checkout gh-pages
    touch .nojekyll
    git add --all .
    git commit -m "Version" --allow-empty -q
    git remote add origin https://$GH_TOKEN@github.com/ARM-DOE/pyart-docs-travis.git &> /dev/null
    git push origin gh-pages -fq &> /dev/null
fi

exit 0
