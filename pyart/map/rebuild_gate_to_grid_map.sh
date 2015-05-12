# remove and rebuild the _load_nn_field_data module
rm  _gate_to_grid_map.so
cython _gate_to_grid_map.pyx
python setup.py build_ext -i
