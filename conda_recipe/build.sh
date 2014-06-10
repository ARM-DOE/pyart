#!/bin/bash

cp -r $RECIPE_DIR/.. .
export RSL_PATH=$PREFIX
$PYTHON setup.py install
