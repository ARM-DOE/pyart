# remove and rebuild _rsl_interface.so and _fourdd_interface.so modules
rm  _rsl_interface.so
cython _rsl_interface.pyx
python setup.py build_ext -i
