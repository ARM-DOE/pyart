# script to rebuild Py-ART example after removing intermediates
rm -r build
rm -r source/source/auto_examples/*
BUILD_PYART_EXAMPLES=1 make html
