#!/usr/bin/env python
"""
This script should produce some initial comparison values with the lattice.

`Imports`: ``gg``, ``fit``, ``integrate``, ``common``, ``sys``, ``os``

Example
-------
>>> python lat.py
>>> python plot.py out/test/f
>>> python plot.py out/test/cs
"""

import gg
from fit import *
from integrate import *
from common import *
from sys import stdout
import os


def test_f(out_dir):
	"""
	plot forward amplitude for initial comparison with 9 lattice points.
	
	Parameters
	----------
	out_dir : str
		output directory for computed data files
	"""
	
	e = EtaPrime()
	p = Pi0()
	w = WholeRegion()
	h = WholeRegion(name='h', add_regge=True)
	
	h.plot_cs(var='nu', xmax=50., dx=0.01, filename=os.path.join(out_dir,'cs_0'), raw=True)
	
	gg.M_pi = 0.451
	gg.f_pi = 0.119
	gg.M_rho = 0.952
	gg.Q2 = 0.369
	
	for gg.Q1 in [0.093, 0.369, 0.828]:
		nu_max = 0.5*(gg.M_pi**2 + gg.Q1 + gg.Q2) - 0.01
		h.plot_cs(var='nu', xmax=50., dx=0.01, filename=os.path.join(out_dir,'cs_%.3f'%gg.Q1), raw=True)
		for r in [e+p, p, w, h]:
			r.plot_f(var='nu', xmax=nu_max, dx=0.01, filename=os.path.join(out_dir,'f%s_%.3f'%(r.name[0],gg.Q1)), raw=True)


if __name__ == '__main__':
	
	test_f(os.path.join('out','test'))
