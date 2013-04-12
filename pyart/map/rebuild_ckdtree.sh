# remove and rebuild the ball_tree.so module
rm  ckdtree.so
cython ckdtree.pyx
python setup.py build_ext -i
