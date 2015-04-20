#!/usr/bin/env python
"""
This script should produce some initial comparison values with the lattice.

`Imports`: ``gg``, ``fit``, ``integrate``, ``common``, ``sys``, ``os``, ``itertools``

Example
-------
>>> python lat.py
>>> python plot.py out/test/f
>>> python plot.py out/test/cs
>>> python plot.py out/out_m_0.451__p_0.093/f
>>> python plot.py out/out_m_0.324__p_0.041/f
"""

import gg
from fit import *
from integrate import *
from common import *
from sys import stdout
import os
from itertools import izip_longest


def make_sure_path_exists(path):
	try: 
		os.makedirs(path)
	except OSError:
		if not os.path.isdir(path):
			raise


def test_m_p(M_pi, f_pi, M_rho, Q1, out_dir='out'):
	r"""
	for a given :math:`M_\pi`, :math:`f_\pi`, :math:`M_\rho` and :math:`Q_1^2`
	plot :math:`Q_2^2`-dependence of the amplitude compared to the lattice results.
	
	Note that Regge is not included here.
	
	Parameters
	----------
	M_pi : float
		pion mass
	f_pi : float
		pion decay constant
	M_rho : float
		:math:`\rho`-meson mass
	Q1 : float
		:math:`Q_1^2`
	out_dir : str, optional
		root output directory
	
	Examples
	--------
	
	.. code-block:: python
	
		from lat.py import test_m_p
		test_m_p(0.451, 0.119, 0.952, 0.093)
	
	run and wait for a while; then plot results:
	
	>>> python plot.py out/out_m_0.451__p_0.093/f
	
	"""
	gg.M_pi, gg.f_pi, gg.M_rho, gg.Q1 = M_pi, f_pi, M_rho, Q1
	#w = TensorMesonRegion(mon_m=0.8) + EtaPrime() + Pi0()
	w = WholeRegion(add_regge=False)
	
	filename = 'm_%.3f__p_%.3f' % (M_pi, Q1)
	with open(os.path.join(out_dir, filename+'.dat')) as p:
		data = [map(float, l.split()) for l in p]
	
	d = dict()
	for nu, Q2, f, err in data:
		if nu != 0.:
			try:
				d[nu] += [[Q2, f, err]]
			except KeyError:
				d[nu] = [[Q2, f, err]]
	
	plots = []
	colors = ['b','g','r','c','m','y','k','grey','purple','orange','violet']
	maxy = {}
	i = 0
	for nu in d.keys():
		maxy[nu] = max([f for _,f,_ in d[nu]])
		
		print 'm = %.3f GeV  nu = %.3f GeV^2' % (M_pi, nu), '...'
		
		w_res = []
		Q2_min = 2*nu - gg.M_pi**2 - gg.Q1 + 0.1
		Q2_max = 6.
		for gg.Q2 in arange(Q2_min, Q2_max, 0.1):
			w_res.append([gg.Q2] + list(w.integral_f(gg.nu2k(nu))))
		
		values = []
		for a, b in izip_longest(d[nu], w_res):
			if a is not None and b is not None:
				values.append(['%.5e' % float(x) for x in a + b])
			elif a is not None:
				values.append(['%.5e' % float(x) for x in a] + ['     -     ']*3)
			elif b is not None:
				values.append(['     -     ']*3 + ['%.5e' % float(x) for x in b])
		
		i = (i + 1) % len(colors)
		make_sure_path_exists(os.path.join(out_dir, 'out_'+filename))
		with open(os.path.join(out_dir, 'out_'+filename, 'nu_%.3f' % nu), 'w') as out:
			out.write(r'xylabel: "$Q_2^2 \, [GeV^2]$" "$Re f - f_0$"' + '\n')
			out.write('color: %s\n' % colors[i])
			out.write('plot: 1 2 yerr=3 linestyle=None marker=s\n')
			out.write(r'plot: 4 5 label="$\nu : \; %.3f$"' % nu + '\n')
			out.write('\n'.join(map('    '.join, values)))
		plots.append('include: "nu_%.3f"' % nu)
	
	with open(os.path.join(out_dir, 'out_'+filename, 'f'), 'w') as out:
		out.write(r'title: $M_{\pi}:\,%.3f\,GeV \;\; Q_1^2:\,%.3f\,GeV^2$' % (M_pi, Q1) + '\n')
		out.write('maxx: 5\nmaxy: %f\n' % (max(maxy.values())+0.00002))
		out.write('\n'.join(plots))


def test_cs(M_pi, f_pi, M_rho, Q1, out_dir='out'):
	r"""
	plot cross sections for a given :math:`M_\pi`, :math:`f_\pi`, :math:`M_\rho` and :math:`Q_1^2`
	and for a set of :math:`Q_2^2` from corresponding lattice data file.
	
	Parameters
	----------
	M_pi : float
		pion mass
	f_pi : float
		pion decay constant
	M_rho : float
		:math:`\rho`-meson mass
	Q1 : float
		:math:`Q_1^2`
	out_dir : str, optional
		root output directory
	
	Examples
	--------
	
	.. code-block:: python
	
		from lat.py import test_cs
		test_cs(0.451, 0.119, 0.952, 0.093)
	
	run and wait for a while; then plot results:
	
	>>> python plot.py out/test_cs/m_0.451__p_0.093
	
	"""
	make_sure_path_exists(out_dir)
	w = WholeRegion(add_regge=True)
	gg.M_pi, gg.f_pi, gg.M_rho, gg.Q1 = M_pi, f_pi, M_rho, Q1
	
	filename = 'm_%.3f__p_%.3f' % (M_pi, Q1)
	with open(os.path.join(out_dir, filename+'.dat')) as p:
		data = [map(float, l.split()) for l in p]
	
	qq = []
	for nu, Q2, f, err in data:
		qq.append(Q2)
	qq = sorted(set(qq))
	#qq = [qq[2*i] for i in range(int(0.5*len(qq)))]
	
	data = [arange(0.001, 4., 0.001)]
	plots = []
	i = 1
	for gg.Q2 in qq:
		i += 1
		data.append([w.cs(gg.nu2k(nu)) for nu in data[0]])
		plots.append('1 %i label="$Q_2^2 :\; %.3f$"' % (i, gg.Q2))
	
	plot_to_file(data, plots, filename=os.path.join(out_dir, 'test_cs', filename), xylabel=r'"$\nu\;[GeV^2]$" "$\sigma\;[\mu b]$"')


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
		nu_max = gg.s2nu(gg.M_pi**2) - 0.01
		h.plot_cs(var='nu', xmax=50., dx=0.01, filename=os.path.join(out_dir,'cs_%.3f'%gg.Q1), raw=True)
		for r in [e+p, p, w, h]:
			r.plot_f(var='nu', xmax=nu_max, dx=0.01, filename=os.path.join(out_dir,'f%s_%.3f'%(r.name[0],gg.Q1)), raw=True)


if __name__ == '__main__':
	
	#test_cs(0.451, 0.119, 0.952, 0.093)
	
	test_f(os.path.join('out','test'))
	
	# plotting comparison for all the points currently available from lattice:
	test_m_p(0.451, 0.119, 0.952, 0.093)
	test_m_p(0.451, 0.119, 0.952, 0.369)
	test_m_p(0.451, 0.119, 0.952, 0.828)
	test_m_p(0.324, 0.109, 0.922, 0.041)
	test_m_p(0.324, 0.109, 0.922, 0.369)
