# remove and rebuild the ball_tree.so module
rm  ball_tree.so
cython ball_tree.pyx
python setup.py build_ext -i
