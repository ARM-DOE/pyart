#!/usr/bin/python
""" Code to plot a Domer RHI file

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

-------------------------------------------------------------------------------


USE
---

HISTORY
-------
5/3/2011:start 
scollis.acrf@gmail.com


"""
__author__ = "Scott Collis <scollis@anl.gov>"
__version__ = "0.1"

import sys
import os
sys.path.append('/home/titan5/python/py4dd')
import matplotlib
# chose a non-GUI backend
matplotlib.use( 'Agg' )
import py4dd
from pylab import *
import numpy as N
import getopt

def create_RHI_array(sweep):
    ppi=N.zeros([sweep.h.nrays, sweep.rays[0].h.nbins], dtype=float)+1.31072000e+05
    for raynum in range(sweep.h.nrays):
	data=sweep.rays[raynum].data
        ppi[raynum, 0:len(data)]=sweep.rays[raynum].data
    return ppi


def radar_coords_to_cart(rng, az, ele, debug=False):
    """
    Asumes standard atmosphere, ie R=4Re/3
    """
    Re=6371.0*1000.0
    #h=(r^2 + (4Re/3)^2 + 2r(4Re/3)sin(ele))^1/2 -4Re/3
    #s=4Re/3arcsin(rcos(ele)/(4Re/3+h))
    p_r=4.0*Re/3.0
    rm=rng*1000.0
    z=(rm**2 + p_r**2 + 2.0*rm*p_r*sin(ele*pi/180.0))**0.5 -p_r
    #arc length
    s=p_r*arcsin(rm*cos(ele*pi/180.)/(p_r+z))
    if debug: print "Z=", z, "s=", s
    y=s*cos(az*pi/180.0)
    x=s*sin(az*pi/180.0)
    return x,y,z

def get_optargs(argv):
  try:
    opts, args=getopt.getopt(argv, "d", ['debug'])
  except getopt.GetoptError:
    usage()
    sys.exit(2)
  return opts, args

def usage():
	print """
	plot_rhi.py -d filename outpath
	should work on UF, Sigmet/IRIS/88D (not that they do RHIs!)
	"""

def plot_rhi(xsapr, imagefilename, var, sweep, **kwargs):
	debug=kwargs.get('debug', False)
	rges=kwargs.get('rges', {'CZ':[-16.,64.], 'ZT':[-16.0, 64.0]}[var])
	var_lab=kwargs.get('var_lab', {'CZ':'Eq refl fact (dBz)', 'ZT':'Eq refl fact (dBz)'}[var])
	print rges
	#if debug: py4dd.RSL_radar_verbose_on()
	field=py4dd.fieldTypes().list.index(var)
	plain_data=create_RHI_array(xsapr.contents.volumes[field].sweeps[sweep])
	azmths=array([xsapr.contents.volumes[field].sweeps[sweep].rays[i].h.azimuth for i in range(xsapr.contents.volumes[field].sweeps[sweep].h.nrays)])
	ranges=xsapr.contents.volumes[field].sweeps[sweep].rays[0].dists/1000.0
	elevs=array([xsapr.contents.volumes[field].sweeps[sweep].rays[i].h.elev for i in range(xsapr.contents.volumes[field].sweeps[sweep].h.nrays)])
	rg,ele=meshgrid(ranges, elevs)
	azg=ones(rg.shape, dtype=float)*xsapr.contents.volumes[field].sweeps[sweep].h.azimuth
	x,y,z=radar_coords_to_cart(rg, azg, ele)
	#xsapr = RSL_anyformat_to_radar('/bm/gscratch4/scollis/cpol_radar/Gunn_pt_20060122153000_PPI_deal.uf') 
	titl=kwargs.get('titl', "RHI %(radar_name)s %(year)04d-%(month)02d-%(day)02d %(hour)02d:%(min)02d:%(sec)02d Azimuth: %(az).2f")
	rad_head=xsapr.contents.h
	titl_comp=titl %{'radar_name':rad_head.radar_name, 'year': rad_head.year,'month': rad_head.month,'day': rad_head.day, 'hour': rad_head.hour, 'min':rad_head.minute, 'sec':rad_head.sec, 'az':azmths[0]}
	f=figure(figsize=[10,4])
	pp=pcolormesh(sign(y)*sqrt(y**2+x**2)/1000.0,z/1000.0,N.ma.masked_outside(plain_data, rges[0],rges[1]) , vmin=rges[0], vmax=rges[1])
	ylim([0,17])
	ylabel('Height above radar (km)')
	xlabel('Distance from radar (km)')
	title(titl_comp)
	cbax=f.add_axes([.9, .1, 0.02, .8])
	c2=colorbar(cax=cbax, mappable=pp)
	ylabel(var_lab)
	savefig(imagefilename)
	close(f)


if __name__ == "__main__":
	pwd=os.path.abspath(os.path.curdir)+'/'
	opts, args=get_optargs(sys.argv[1:])
	opt_dict=dict(opts)
	if ('--debug' in opt_dict.keys()) | ('-d' in opt_dict.keys()):
		debug=True
	else:
		debug=False
	if debug:
		print "plot_rhi.py "+ __version__
	try:
		filename=args[0]
	except IndexError:
		usage()
		sys.exit(2)
	if debug: print "file name=", filename
	if debug: py4dd.RSL_radar_verbose_on()
	xsapr = py4dd.RSL_anyformat_to_radar(filename)
	plot_rhi(xsapr, '/home/titan5/python/test_rhi0.png', 'CZ',0, rges=[-32,64], titl="RHI " +filename.split('/')[-1][0:3]+ "-%(radar_name)sp %(year)04d-%(month)02d-%(day)02d %(hour)02d:%(min)02d:%(sec)02d Azimuth: %(az).2f")
	#plot_rhi(xsapr, '/home/titan5/python/test_rhi1.png', 'CZ',1, rges=[-32,64], titl="RHI " +filename.split('/')[-1][0:3]+ "-%(radar_name)sp %(year)04d-%(month)02d-%(day)02d %(hour)02d:%(min)02d:%(sec)02d Azimuth: %(az).2f")
	#plot_rhi(xsapr, '/home/titan5/python/test_rhi2.png', 'CZ',2, rges=[-32,64], titl="RHI " +filename.split('/')[-1][0:3]+ "-%(radar_name)sp %(year)04d-%(month)02d-%(day)02d %(hour)02d:%(min)02d:%(sec)02d Azimuth: %(az).2f")
	#plot_rhi(xsapr, '/home/titan5/python/test_rhi3.png', 'CZ',3, rges=[-32,64], titl="RHI " +filename.split('/')[-1][0:3]+ "-%(radar_name)sp %(year)04d-%(month)02d-%(day)02d %(hour)02d:%(min)02d:%(sec)02d Azimuth: %(az).2f")
	#plot_rhi(xsapr, '/home/titan5/python/test_rhi4.png', 'CZ',4, rges=[-32,64], titl="RHI " +filename.split('/')[-1][0:3]+ "-%(radar_name)sp %(year)04d-%(month)02d-%(day)02d %(hour)02d:%(min)02d:%(sec)02d Azimuth: %(az).2f")
	#plot_rhi(xsapr, '/home/titan5/python/test_rhi5.png', 'CZ',5, rges=[-32,64], titl="RHI " +filename.split('/')[-1][0:3]+ "-%(radar_name)sp %(year)04d-%(month)02d-%(day)02d %(hour)02d:%(min)02d:%(sec)02d Azimuth: %(az).2f")
	py4dd.RSL_free_radar(xsapr)