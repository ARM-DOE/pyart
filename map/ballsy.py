"""
Ballsy wrapper

USE
---
Takes in flattened arrays of data, xyz locations of the data and maps to a
regular grid using a distance weighted scheme (currently a Cressman and
Barnes scheme are implimented)

"""

from time import time

import numpy as np

import ball_tree


class BallsyMapper:

    def __init__(self, X, z, leafsize=10, debug=False):

        assert len(X) == len(z), "len(X) %d != len(z) %d" % (len(X), len(z))

        if debug:
            print "Filling out the Ball Tree"

        # build the tree
        self.tree = ball_tree.BallTree(X, leafsize=leafsize)
        self.z = z
        self.X = X

    def __call__(self, q, r, debug=False, func='Cressman'):
        q = np.asarray(q)
        qdim = q.ndim
        if qdim == 1:
            q = np.array([q])
        t1 = time()
        if debug:
            print "Looking up nearest neighbours for points to be mapped to"
        self.locations, self.distances = self.tree.query_radius(
            q, r, return_distance=True)
        t2 = time()
        if debug:
            print "Removing voids"
        isgood = np.where(np.array([len(x) for x in self.locations]) != 0)[0]
        interpol = np.empty((len(self.locations), ) + np.shape(self.z[0]))
        interpol.fill(np.nan)
        jinterpol = 0
        t3 = time()
        if debug:
            print "Performing weighting calculations"
        my_distance = self.distances[isgood] * self.distances[isgood]
        locs = self.locations[isgood]
        my_r = r[isgood] * r[isgood]
        if func.upper() == 'CRESSMAN':
            for dist, ix, posn, roi in zip(my_distance, locs, isgood, my_r):
                interpol[isgood[jinterpol]] = np.average(
                    self.z[ix], weights=(roi - dist) / (roi + dist), axis=0)
                jinterpol += 1
        elif func.upper() == 'BARNES':
            for dist, ix, posn, roi in zip(my_distance, locs, isgood, my_r):
                w = np.exp(-dist / 2.0 * roi) + 1e-5
                w /= np.sum(w)
                interpol[isgood[jinterpol]] = np.ma.dot(w, self.z[ix])
                jinterpol += 1
        else:
            print "Sorry chap, not implimented"
        t4 = time()
        if debug:
            print "Time to query the tree:", t2 - t1
            print "time to remove voids:", t3 - t2
            print "time to calculate weights and map:", t4 - t3
        return interpol
