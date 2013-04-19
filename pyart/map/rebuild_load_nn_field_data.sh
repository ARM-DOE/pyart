# remove and rebuild the _load_nn_field_data module
rm  _load_nn_field_data.so
cython _load_nn_field_data.pyx
python setup.py build_ext -i
