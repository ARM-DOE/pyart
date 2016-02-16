# script to rebuild complete documentation include examples after removing
# intermediates
rm -r build
rm source/user_reference/generated/*
rm source/dev_reference/generated/*
rm -r source/auto_examples/*
BUILD_PYART_EXAMPLES=1 make html
