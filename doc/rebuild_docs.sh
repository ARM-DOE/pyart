# script to rebuild documentation aftter removing intermediates
rm -r build
rm source/user_reference/generated/*
rm source/dev_reference/generated/*
make html

