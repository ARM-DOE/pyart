#!/bin/bash
# Adapted from the ci/build_docs.sh file from the pandas project
# https://github.com/pydata/pandas

cd "$TRAVIS_BUILD_DIR"

echo "Building Docs"
conda install --yes sphinx pil

mv "$TRAVIS_BUILD_DIR"/doc /tmp
cd /tmp/doc
make html

# upload to pyart-docs-travis repo is this is not a pull request
if [ "$TRAVIS_PULL_REQUEST" == "false" ]; then
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
    git commit -m "Version" --allow-empty
    git remote add origin https://$GH_TOKEN@github.com/jjhelmus/pyart-docs-travis
    git push origin gh-pages -f
fi

exit 0
