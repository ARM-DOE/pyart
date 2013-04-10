# remove and rebuild the _rsl_interface.so python module
rm  _rsl_interface.so
cython _rsl_interface.pyx
python setup.py build_ext -i
