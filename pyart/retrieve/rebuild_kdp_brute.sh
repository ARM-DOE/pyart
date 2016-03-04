# remove and rebuild the kdp_brute.so module
rm -f kdp_brute.so
cython kdp_brute.pyx
python setup.py build_ext -i
