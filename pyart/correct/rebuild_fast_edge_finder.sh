# remove and rebuild the fast.so module
rm  _fast_edge_finder.so
cython _fast_edge_finder.pyx
python setup.py build_ext -i
