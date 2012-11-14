""" Utilities for reading of MDV data into numpy objects

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Argonne National Laboratory nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL Scott Collis BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

-------------------------------------------------------------------------------


Scott Collis, Argonne National Laboratory, 
Computing Environment and Life Sciences Directorate
ARM Climate Research Facility

United States Department of Energy

USE
---
Takes in flattened arrays of data, xyz locations of the data and maps to a regular grid using a distance weighted scheme (currently a Cressman and Barnes scheme are implimented


REQUIREMENTS
------------
Needs Numpy, time, ball_tree (Cython/c++ class from scikits-learn)

HISTORY
-------
2011-05-23 Start of development 
Scott Collis scollis.acrf@gmail.com
0.1: basic functionality
2011-07-19:
1.0 integrated into Py-ART
2012-10-22:
1.1 added multi-platform import of pyart
"""
__author__ = "Scott Collis <scollis@anl.gov>"
__version__ = "1.0"
import sys
import os
import numpy as np
pyart_dir=os.environ.get('PYART_DIR',os.environ['HOME']+'/python')
sys.path.append(pyart_dir+'/pyart/map/ball_tree/ball_tree')
#sys.path.append('/home/sc8/python/pyart/map/ball_tree/ball_tree')
#sys.path.append('/home/titan5/python/pyorder/cress')
import ball_tree
#import cress
from time import time


class BallsyMapper:
    def __init__( self, X, z, leafsize=10, debug=False):
        assert len(X) == len(z), "len(X) %d != len(z) %d" % (len(X), len(z))
	if debug: print "Filling out the Ball Tree"
        self.tree = ball_tree.BallTree( X, leafsize=leafsize )  # build the tree
        self.z =z
	self.X=X
    def __call__( self, q, r, debug=False, func='Cressman'):
        q = np.asarray(q)
        qdim = q.ndim
        if qdim == 1:
            q = np.array([q])
        t1=time()
	if debug: print "Looking up nearest neighbours for points to be mapped to"
        self.locations, self.distances = self.tree.query_radius( q, r, return_distance=True)
	t2=time()
	#tp=np.zeros(len(self.locations), dtype=float)
	#for i in range(len(self.locations)):
	#	tp[i]=len(self.locations[i])
	#print "Mean:", tp.mean(), " Max:", tp.max(), " Total:", tp.sum()
	if debug: print "Removing voids"
	isgood=np.where( np.array([len(x) for x in self.locations])!=0)[0]
	interpol = np.empty( (len(self.locations),) + np.shape(self.z[0]) )
	interpol.fill(np.nan)
	jinterpol=0
	t3=time()
	if debug: print "Performing weighting calculations"
	#for i in isgood:
	#interpolr=cress.cress(self.z, self.distances, r, self.locations, isgood, interpol)# cress(np.ndarray z,np.ndarray d,np.ndarray r,np.ndarray l, np.ndarray isgood, np.ndarray interpol)
	#interpolm =[ np.average(self.z[self.locations[i]], weights=(r[i]**2-self.distances[i]**2) / (r[i]**2 + self.distances[i]**2), axis=0) for i in isgood]
	#interpol[isgood]=interpolm
	my_distance=self.distances[isgood]*self.distances[isgood]
	locs=self.locations[isgood]
	my_r=r[isgood]*r[isgood]
	if func.upper()=='CRESSMAN':
		for dist, ix, posn, roi in zip(my_distance, locs, isgood, my_r):
			interpol[isgood[jinterpol]] = np.average(self.z[ix], weights=(roi-dist) / (roi + dist), axis=0)
			jinterpol += 1
	elif func.upper()=='BARNES':
		for dist, ix, posn, roi in zip(my_distance, locs, isgood, my_r):
			w=np.exp(-dist / 2.0*roi)+1e-5
			w/=np.sum(w)
			interpol[isgood[jinterpol]] = np.ma.dot(w, self.z[ix])#np.average(self.z[ix], weights=np.exp(-dist / 2.0*roi)+1e-5, axis=0)
			jinterpol += 1
	else:
		print "Sorry chap, not implimented"
	t4=time()
	if debug:
		print "Time to query the tree:", t2-t1
		print "time to remove voids:", t3-t2
		print "time to calculate weights and map:", t4-t3
	return interpol
