# script to rebuild complete documentation include examples after removing
# intermediates
rm -r build
rm source/API/generated/*
rm -r source/examples/*
BUILD_PYART_EXAMPLES=1 make html
