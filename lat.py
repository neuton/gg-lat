#!/usr/bin/env python
"""
This script should produce some initial comparison values with the lattice.

`Imports`: ``gg``, ``fit``, ``integrate``, ``common``, ``sys``, ``os``, ``itertools``

Example
-------
>>> python lat.py
>>> python plot.py out/test1/mpi_451_Q2_094/f
>>> python plot.py out/test2/mpi_451_Q1_094/f
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


def get_lattice_data_filename(M_pi, Q1, s='mpi_%03d_Q1_%03d'):
	return s % (int(1000*M_pi), int(1000*Q1))


def test1_mpi_q2(M_pi, f_pi, M_rho, Q2, out_dir='out'):
	r"""
	for a given :math:`M_\pi`, :math:`f_\pi`, :math:`M_\rho` and :math:`Q_2^2`
	plot :math:`\nu`-dependence of the amplitude compared to interpolated lattice results.
	
	Parameters
	----------
	M_pi : float
		pion mass
	f_pi : float
		pion decay constant
	M_rho : float
		:math:`\rho`-meson mass
	Q2 : float
		:math:`Q_2^2`
	out_dir : str, optional
		root output directory
	
	Examples
	--------
	
	.. code-block:: python
	
		from lat.py import test1_mpi_q2
		test1_mpi_q2(0.451, 0.119, 0.952, 0.094)
	
	run and wait for a while; then plot results:
	
	>>> python plot.py out/test1/mpi_451_Q2_094/f
	
	"""
	
	e = EtaPrime()
	p = Pi0()
	w = WholeRegion()
	h = WholeRegion(name='h', add_regge=True)
	
	gg.M_pi, gg.f_pi, gg.M_rho, gg.Q2 = M_pi, f_pi, M_rho, Q2
	
	filename = get_lattice_data_filename(M_pi, Q2, s='mpi_%03d_Q2_%03d')
	with open(os.path.join('lattice', filename + '__int')) as fi:
		data = []
		for l in fi:
			try:
				data .append(map(float, l.split()))
			except ValueError:
				pass
	
	d = dict()
	for nu, Q1, f, err in data:
		if Q1 != 0.:
			try:
				d[Q1] += [[nu, f, err]]
			except KeyError:
				d[Q1] = [[nu, f, err]]
	
	out_dir = os.path.join(out_dir, 'test1', filename)
	make_sure_path_exists(out_dir)
	
	h.plot_cs(var='nu', xmax=50., dx=0.01, filename=os.path.join(out_dir, 'cs_0'), raw=True)
	
	plots = []
	for gg.Q1 in d.keys():
		h.plot_cs(var='nu', xmax=50., dx=0.01, filename=os.path.join(out_dir,'cs_%.3f'%gg.Q1), raw=True)
		nu_max = gg.s2nu(gg.M_pi**2) - 0.01
		for r in [e+p, p, w, h]:
			r.plot_f(var='nu', xmax=nu_max, dx=0.01, filename=os.path.join(out_dir,'f%s_%.3f'%(r.name[0],gg.Q1)), raw=True)
		with open(os.path.join(out_dir, 'f_%.3f' % gg.Q1), 'w') as fo:
			fo.write('\n'.join([' '.join(map(str, l)) for l in d[gg.Q1]]))
	
	colors = ['r','g','b','c','m','y','k','grey','purple','orange','violet']
	colors = (colors * (1 + len(d.keys())/len(colors)))[:len(d.keys())]
	
	with open(os.path.join(out_dir, 'cs'), 'w') as fo:
		fo.write(r'xylabel: "$\nu \; [GeV ^2]$" "$\sigma \; [\mu b]$"' + '\n')
		fo.write('maxx: 3\nmaxy: 0.9\n')
		fo.write('plot: 1 2 source=cs_0     color=grey  label=physical\n')
		for Q1, color in zip(d.keys(), colors):
			fo.write('plot: 1 2 source=cs_%.3f color=%s label=' % (Q1, color) + r'"$Q_1^2 : \, ' + '%.3f'%Q1 + r' \, GeV^2$"' + '\n')
	
	with open(os.path.join(out_dir, 'f'), 'w') as fo:
		fo.write(r"xylabel: '$\nu \, [GeV^2]$' '$Re\,f - f_0$'" + '\n')
		fo.write('maxx: 0.7\nmaxy: 0.00012\nlinewidth: 1.5\n')
		for Q1, color in zip(d.keys(), colors):
			fo.write('plot: 1 2 source=fw_%.3f color=%s '%(Q1, color) + r'label="$Q_1^2 : \, ' + '%.3f'%Q1 + r' \, GeV^2$"' + '\n')
			fo.write('plot: 1 2 source=fh_%.3f color=%s linestyle=-.\n' % (Q1, color))
			fo.write('plot: 1 2 source=fp_%.3f color=%s linestyle=--\n' % (Q1, color))
			fo.write('plot: 1 2 source=fe_%.3f color=%s linestyle=:\n' % (Q1, color))
			fo.write('plot: 1 2 yerr=3 source=f_%.3f marker=o linestyle=None color=%s\n' % (Q1, color))


def test2_mpi_q1(M_pi, f_pi, M_rho, Q1, out_dir='out'):
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
	
		from lat.py import test2_mpi_q1
		test2_mpi_q1(0.451, 0.119, 0.952, 0.094)
	
	run and wait for a while; then plot results:
	
	>>> python plot.py out/test2/mpi_451_Q1_094/f
	
	"""
	#w = TensorMesonRegion(mon_m=0.8) + EtaPrime() + Pi0()
	w = WholeRegion(add_regge=False)
	
	gg.M_pi, gg.f_pi, gg.M_rho, gg.Q1 = M_pi, f_pi, M_rho, Q1
	
	filename = get_lattice_data_filename(M_pi, Q1)
	data = []
	with open(os.path.join('lattice', filename)) as p:
		for l in p:
			try:
				data .append(map(float, l.split()))
			except ValueError:
				pass
		
	d = dict()
	for nu, Q2, f, err in data:
		if nu != 0.:
			try:
				d[nu] += [[Q2, f, err]]
			except KeyError:
				d[nu] = [[Q2, f, err]]
	
	out_dir = os.path.join(out_dir, 'test2', filename)
	make_sure_path_exists(out_dir)
	
	plots = []
	colors = ['r','g','b','c','m','y','k','grey','purple','orange','violet']
	maxy = {}
	miny = {}
	i = 0
	for nu in d.keys():
		maxy[nu] = max([f for _,f,_ in d[nu]])
		miny[nu] = min([f for _,f,_ in d[nu]])
		
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
		with open(os.path.join(out_dir, 'nu_%.3f' % nu), 'w') as out:
			out.write(r'xylabel: "$Q_2^2 \, [GeV^2]$" "$Re f - f_0$"' + '\n')
			out.write('color: %s\n' % colors[i])
			out.write('plot: 1 2 yerr=3 linestyle=None marker=s\n')
			out.write(r'plot: 4 5 label="$\nu : \; %.3f$"' % nu + '\n')
			out.write('\n'.join(map('    '.join, values)))
		plots.append('include: "nu_%.3f"' % nu)
	
	with open(os.path.join(out_dir, 'f'), 'w') as out:
		out.write(r'title: $M_{\pi}:\,%.3f\,GeV \;\; Q_1^2:\,%.3f\,GeV^2$' % (M_pi, Q1) + '\n')
		out.write('maxx: 5\nminy: %f\nmaxy: %f\n' % (min(maxy.values())-0.00001, max(maxy.values())+0.00001))
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
		test_cs(0.451, 0.119, 0.952, 0.094)
	
	run and wait for a while; then plot results:
	
	>>> python plot.py out/test_cs/mpi_451_Q1_094
	
	"""
	make_sure_path_exists(out_dir)
	w = WholeRegion(add_regge=True)
	#w = PiPiRegion()
	gg.M_pi, gg.f_pi, gg.M_rho, gg.Q1 = M_pi, f_pi, M_rho, Q1
	
	filename = get_lattice_data_filename(M_pi, Q1)
	with open(os.path.join('lattice', filename)) as p:
		data = []
		for l in p:
			try:
				data .append(map(float, l.split()))
			except ValueError:
				pass
	
	qq = []
	for nu, Q2, f, err in data:
		qq.append(Q2)
	qq = sorted(set(qq))
	#qq = [qq[2*i] for i in range(int(0.5*len(qq)))]
	#qq = [qq[-1]]
	
	data = [arange(0.001, 4., 0.001)]
	plots = []
	i = 1
	for gg.Q2 in qq:
		i += 1
		data.append([w.cs(gg.nu2k(nu)) for nu in data[0]])
		plots.append('1 %i label="$Q_2^2 :\; %.3f$"' % (i, gg.Q2))
	
	make_sure_path_exists(os.path.join(out_dir, 'test_cs'))
	plot_to_file(data, plots, filename=os.path.join(out_dir, 'test_cs', filename), xylabel=r'"$\nu\;[GeV^2]$" "$\sigma\;[\mu b]$"')


if __name__ == '__main__':
	
	# see cross sections extrapolation:
	test_cs(0.451, 0.119, 0.952, 0.094)
	
	# M_pi = 277 MeV:
	m = 0.277, 0.104, 0.878
	test1_mpi_q2(*m, Q2=0.377)
	test2_mpi_q1(*m, Q1=0.042)
	test2_mpi_q1(*m, Q1=0.377)
	
	# M_pi = 324 MeV:
	m = 0.324, 0.109, 0.922
	test1_mpi_q2(*m, Q2=0.377)
	test2_mpi_q1(*m, Q1=0.042)
	test2_mpi_q1(*m, Q1=0.377)
	
	# M_pi = 451 MeV:
	m = 0.451, 0.119, 0.952
	test1_mpi_q2(*m, Q2=0.377)
	test2_mpi_q1(*m, Q1=0.094)
	test2_mpi_q1(*m, Q1=0.377)
	test2_mpi_q1(*m, Q1=0.845)
