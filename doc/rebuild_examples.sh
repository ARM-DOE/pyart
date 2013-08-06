# script to rebuild documentation aftter removing intermediates
rm -r build
rm -r source/auto_examples/*
make html

