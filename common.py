"""
some common helper functions definitions

`Imports`: ``math``, ``itertools``

------------------------------------------------------------------------
"""


def read_xy(filename, xmin = None, xmax = None):
	"""
	read first 2 columns of numbers (assumed to be X-values and Y-values, respectively) from file.
	
	Parameters
	----------
	filename : string
		a file to read numbers from
	xmin : float, optional
		minimum X-value to read from
	xmax : float, optional
		maximum X-value to read to
	
	Returns
	-------
	list
		a list of 2 tuples of floats
	"""
	data = []
	with open(filename) as fi:
		for l in fi:
			if l.strip()!='' and l.strip()[0]!='#':
				try:
					data.append(map(float, l.split()[:2]))
				except:
					pass
	compare = lambda a, b: 1 if a[0]>b[0] else -1
	data = sorted(data, compare)
	r = zip(*[(d[0], d[1]) for d in data if (xmin or d[0]) <= d[0] <= (xmax or d[0])])
	return [(),()] if r==[] else r


def read_xye(filename, xmin = None, xmax = None):
	"""
	read first 3 columns of numbers (assumed to be X-values, Y-values and error-values, respectively) from file.
	
	Parameters
	----------
	filename : string
		a file to read numbers from
	xmin : float, optional
		minimum X-value to read from
	xmax : float, optional
		maximum X-value to read to
	
	Returns
	-------
	list
		a list of 3 tuples of floats
	"""
	data = []
	with open(filename) as fi:
		for l in fi:
			if l.strip()!='' and l.strip()[0]!='#':
				try:
					data.append(map(float, l.split()[:3]))
				except:
					pass
	compare = lambda a, b: 1 if a[0]>b[0] else -1
	data = sorted(data, compare)
	r = zip(*[(d[0], d[1], 0. if len(d) < 3 else d[2]) for d in data if (xmin or d[0]) <= d[0] <= (xmax or d[0])])
	return [(),(),()] if r==[] else r


import math

def arange(a, b, dx, logx=False):
	"""
	generate a list of float values between `a` and `b` with a step-size `dx`.
	
	Parameters
	----------
	a : float
		minimum value
	b : float
		maximum value
	dx : float
		step-size
	logx : bool, optional
		a flag for logarithmic stepping
	
	Returns
	-------
	list
		a list of floats
	"""
	if a > b:
		return []
	elif logx:
		return [math.exp(math.log(a) + dx*i) for i in range(int((math.log(b)-math.log(a))/dx))] + [b]
	else:
		return [a + dx*i for i in range(int((b-a)/dx))] + [b]


from itertools import izip_longest

def plot_to_file(data, plots=[], filename='out/out', title=None, xylabel=None, xlogscale=False):
	"""
	write columns of floats from given data to a file for further plotting with ``plot.py``
	
	Parameters
	----------
	data : iterable
		a list of lists (columns) of floats
	plots : list, optional
		a list of strings of plotting directives parsed by ``plot.py``
	filename : string, optional
		a file to write to
	title : string, optional
		a figure title parsed by ``plot.py``
	xylabel : string, optional
		contains labels for X-axis and Y-axis parsed by ``plot.py``
	xlogscale : bool, optional
		a flag to inform `plot.py` whether to use logarithmic X-scale
	"""
	with open(filename, 'w') as fout:
		if xlogscale:
			fout.write('xlogscale:\n')
		if title is not None:
			fout.write('title: ' + title + '\n')
		if xylabel is not None:
			fout.write('xylabel: ' + xylabel + '\n')
		for p in plots:
			fout.write('plot: ' + p + '\n')
		fout.write('\n'.join([' '.join(map(str, l)) for l in izip_longest(*data)]))
