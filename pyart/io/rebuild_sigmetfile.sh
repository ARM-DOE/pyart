# remove and rebuild _sigmetfile.so module
rm  _sigmetfile.so
cython _sigmetfile.pyx
python setup.py build_ext -i
