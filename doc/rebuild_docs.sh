# script to rebuild documentation after removing intermediates
rm -r build
rm source/API/generated/*
make html
