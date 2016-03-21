# remove and rebuild the _kdp_proc.so module
rm -fv _kdp_proc.so
cython _kdp_proc.pyx
python setup.py build_ext -i
