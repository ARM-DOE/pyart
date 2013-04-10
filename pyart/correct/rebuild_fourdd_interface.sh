# remove and rebuild the _fourdd_interface.so module
rm  _fourdd_interface.so
cython _fourdd_interface.pyx
python setup.py build_ext -i
