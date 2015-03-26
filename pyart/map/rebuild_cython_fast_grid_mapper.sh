# remove and rebuild the ball_tree.so module
rm  ckdtree.so
cython cython_fast_grid_mapper.pyx
python setup.py build_ext -i
