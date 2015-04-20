r"""
Integration module for :math:`\gamma\gamma` cross section.

Defines classes for forward :math:`\gamma\gamma` scattering amplitude computation.
By "amplitude" I mean :math:`Re(f-f_0)`.

`Imports`: ``scipy``, ``sys``, ``fit``, ``gg``, ``common``

Examples
--------
- see lat.py__
- plot some amplitude and cross-section:

	``test.py``:

	.. code-block:: python
	
		w = WholeRegion(add_regge=False) # create the whole region instance without Regge
		w.plot_cs(dx=0.01, xmax=10., filename='out/cs_0_0_0') # plot cross-section
		gg.f_pi, gg.M_pi, gg.M_rho = 0.119, 0.451, 0.952 # set new pion mass
		gg.Q2 = gg.Q1 = 0.369 # set new virtualies
		w.plot_cs(dx=0.01, xmax=10., filename='out/cs_0.451_0.37_0.37') # plot cross-section again
		w.plot_f(filename='out/f_0.451_0.37_0.37', xmin=1.5, xmax=2.2, dx=0.001) # plot some amplitude (slow)

__ lat.html

------------------------------------------------------------------------
"""

import gg
from common import *
from fit import *

from scipy import pi, sqrt, inf, exp
from scipy.integrate import simps, trapz, quad
from scipy.misc import derivative as diff
from sys import stdout



class Region:
	r"""
	base class for integration.
	
	This class represents an energy region for cross section integration.
	It holds fit and data from the region and is intended to simplify amplitude computation procedure.
	
	Note that :math:`\pi` and :math:`\eta'` :math:`\delta`-contributions are initially *not* included here.
	
	Parameters
	----------
	function_fit : FunctionFit, optional
		a fit for the region
	k_min : float, optional
		initial *physical* value of lower energy bound of the region
	k_max : float, optional
		initial *physical* value of upper energy bound of the region
	name : str, optional
		name of the region to be displayed on plots
	init : bool, optional
		choose whether to fit the data on instance creation
	"""
	def __init__(self, function_fit=None, k_min=0., k_max=inf, name='region', init=True):
		if function_fit is None:
			self.fit = FunctionFit(p=[0., 0.])
		else:
			self.fit = function_fit
		self.k_min, self.k_max = k_min, k_max
		self.name = name
		self.update_data()
		if init:
			self.update_fit()
	
	def update_data(self, filename='data/data.dat'):
		"""
		update data for fitting from file.
		
		Parameters
		----------
		filename : str, optional
			name of the file to read data from
		
		Returns
		-------
		Region
			self
		"""
		self.x, self.y, self.err = map(list, read_xye(filename, *self.k_bounds()))
		return self
	
	def update_fit(self, quiet=False):
		"""
		update parameters of fit for the data in the region.
		
		Parameters
		----------
		quiet : bool, optional
			whether to print info
		
		Returns
		-------
		Region
			self
		"""
		if len(self.x) > 4:
			self.fit.fit(self.x, self.y, self.err)
		if not quiet:
			self.chisqr = chisqr(self.fit, self.x, self.y, self.err)
			print 'total chisqr/point =', self.chisqr/len(self.x)
		return self
	
	def get_fit_err(self):
		"""
		get fitting errors estimation.
		
		Returns
		-------
		list
			list of errors for each X-value of data in the region.
		"""
		if len(self.x) > 4:
			return self.fit.get_errors(self.x, self.y, self.err)
		else:
			return [0. for e in self.x]
	
	def __add__(self, region):
		"""
		create new region with sum of cross sections of 2 regions.
		
		Warning: this is broken by design.
		
		Parameters
		----------
		region : Region
			second region to add self to
		
		Returns
		-------
		Region
			sum of self + region
		"""
		sum_region = Region(name=self.name+' + '+region.name, init=False)
		sum_region.cs = lambda k: map(sum, zip(self.cs(k), region.cs(k)))
		sum_region.integral_f = lambda k: map(sum, zip(self.integral_f(k), region.integral_f(k)))
		return sum_region
	
	def k_bounds(self):
		"""
		get current bounds of the region.
		
		This method should be called when parameters from ``gg`` are changed.
		This should most probably be overriden when inherited.
		
		Returns
		-------
		tuple
			(k_min, k_max)
		"""
		return 0., inf
	
	def cs(self, k):
		r"""
		get cross section value.
		
		Parameters
		----------
		k : float
			:math:`\sqrt{2 q_1 \cdot q_2}` in GeV
		
		Returns
		-------
		float
			fit value at specified `k`
		"""
		return self.fit(k)
	
	def integral_f(self, k):
		r"""
		integrate cross section to obtain forward amplitude:
		
		.. math:: Re(f(\nu)-f(0)) = \frac{8\nu^2}{\pi} \int_{\nu_{min}}^{\nu_{max}}
			\frac{\sqrt{X(\nu')}}{\nu'} \frac{d\nu'}{(\nu'^2-\nu^2)} \sigma_{TT}(\nu'),
		
		where :math:`\nu_{min}` and :math:`\nu_{max}` are bounds of the integration region.
		
		Note that :math:`\pi` and :math:`\eta'` :math:`\delta`-contributions are *not* included here.
		
		Parameters
		----------
		k : float
			:math:`\sqrt{2 q_1 \cdot q_2}, \; [GeV]`
		
		Returns
		-------
		tuple
			value of Re(f(k)-f(0)), some error
		
		Notes
		-----
		Divergency extraction is done here (when :math:`\nu` is close enough to threshold):
		
		.. math:: I = \frac{8\nu^2}{\pi} \left(\int_{\nu_{min}}^{\nu_{max}} \frac{(a(\nu') - a(\nu))}{(\nu'^2-\nu^2)} d\nu' \; + \;
			\frac{a(\nu)}{2\nu} \ln\left|\frac{\nu'-\nu}{\nu'+\nu}\right| \Big|_{\nu_{min}}^{\nu_{max}} \right),
		
		where
		
		.. math:: a(\nu) = \frac{\sqrt{X(\nu)}}{\nu}\sigma(\nu).
		"""
		nu = gg.k2nu(k)
		nu0, nu1 = map(gg.k2nu, self.k_bounds())
		a = lambda x: self.cs(x) * gg.nu2X(x)**0.5 / x
		if gg.nu2X(nu) <= 0:
			# temporar hack: separating integration region
			q1 = quad(lambda p: a(p)/(p**2 - nu**2), nu0, 4*nu0)[:2]
			q2 = quad(lambda p: a(p)/(p**2 - nu**2), 4*nu0, nu1)[:2]
			r, err = map(sum, zip(q1, q2))
		else:
			a_nu = a(nu)
			singular = lambda p: 0. if p == inf else log(abs((p-nu)/(p+nu)))/(2*nu)
			r, err = quad(lambda p: (a(p) - a_nu)/(p**2 - nu**2), nu0, nu1)[:2]
			r +=  a_nu * (singular(nu1) - singular(nu0))
		prefactor = 8 * nu**2 / (pi * 10**4 * hc**2) # don't forget to convert units
		return prefactor * r, prefactor * err
	
	def plot_cs(self, xmin=None, xmax=None, filename=None, dx=None, var='k', color=None, raw=False, logx=False):
		"""
		put cross section values to file for further plotting.
		
		Parameters
		----------
		xmin : float, optional
			minimum X value for output points
		xmax : float, optional
			maximum X value for output points
		filename : str, optional
			output filename
		dx : float, optional
			step-size for X-values
		var : str, optional
			what X variable should stand for;
			available values: 'k'(default), 's', 'nu'
		color : str, optional
			a color to plot cross section with;
			ignored when `raw` is `True`
		raw : bool, optional
			whether to put raw columns of values to the file
			or to add plotting directives for ``plot.py`` there as well
		logx : bool, optional
			whether to produce values with logarithmic X-scale
		
		Returns
		-------
		Region
			self
		"""
		k_min, k_max = self.k_bounds()
		if var == 's':
			convert2k = gg.s2k
		elif var == 'nu':
			convert2k = gg.nu2k
		else:
			convert2k = gg.x2x
		if filename is None:
			filename = 'out/cs_'+self.name+'.out'
		if xmin is None:
			xmin = gg.reverse_convert(convert2k)(k_min)
		if xmax is None:
			if k_max is inf:
				xmax = 10**8
			else:
				xmax = gg.reverse_convert(convert2k)(k_max)
		#data = zip(*[(x, y, e) for x, y, e in zip(self.x, self.y, self.err) if xmin <= x <= xmax])
		if dx is None:
			xf = [x for x in self.x if xmin <= x <= xmax]
		else:
			xf = arange(xmin, xmax, dx, logx)
		fx = []
		stdout.write('\r[0%  ] plotting to "'+filename+'"')
		for i in range(len(xf)):
			stdout.write('\r[%d%%' % (100*(i+1)/len(xf)))
			stdout.flush()
			fx.append(self.cs(convert2k(xf[i])))
		stdout.write('\n')
		if not raw:
			plots = ['1 2 label="'+self.name+' fit" linewidth=2']
			if color is not None:
				plots[0] += 'color="'+color+'"'
			title = 'Total cross-section'
			xylabel = var+' "total cross-section [mcb]"'
		else:
			plots = []
			title = xylabel = None
		plot_to_file([xf, fx], plots, filename, title=title, xylabel=xylabel, xlogscale=logx)
		return self
	
	def plot_f(self, xmin=0., xmax=5., filename=None, dx=0.005, var='k', raw=False, logx=False):
		"""
		put amplitude values to file for further plotting.
		
		Note that this is usually slow, so don't set too small `dx`.
		
		Parameters
		----------
		xmin : float, optional
			minimum X value for output points
		xmax : float, optional
			maximum X value for output points
		filename : str, optional
			output filename
		dx : float, optional
			step-size for X-values
		var : str, optional
			what X variable should stand for;
			available values: 'k'(default), 's', 'nu'
		raw : bool, optional
			whether to put raw columns of values to the file
			or to add plotting directives for ``plot.py`` there as well
		logx : bool, optional
			whether to produce values with logarithmic X-scale
		
		Returns
		-------
		Region
			self
		"""
		if var == 's':
			convert2k = gg.s2k
		elif var == 'nu':
			convert2k = gg.nu2k
		else:
			convert2k = gg.x2x
		if filename is None:
			filename = 'out/f_'+self.name+'.out'
		xf = arange(xmin + dx, xmax, dx, logx)
		fx = []
		err = []
		stdout.write('\r[0%  ] plotting to "'+filename+'"')
		for i in range(len(xf)):
			stdout.write('\r[%d%%' % (100*(i+1)/len(xf)))
			stdout.flush()
			int_f, int_err = self.integral_f(convert2k(xf[i]))
			fx.append(int_f)
			err.append(int_err)
		stdout.write('\n')
		if not raw:
			plots = ['1 2 label="'+self.name+'"']
			title = self.name
			xylabel = var + ' "Re f"'
		else:
			plots = []
			title = xylabel = None
		plot_to_file([xf, fx, err], plots=plots, filename=filename, title=title, xylabel=xylabel, xlogscale=logx)
		return self
	
	def print_f(self, k):
		"""
		print value of amplitude.
		
		Parameters
		----------
		k : float
			:math:`\sqrt{2 q_1 \cdot q_2}, \; [GeV]`
		
		Returns
		-------
		Region
			self
		"""
		a, e = self.integral_f(k)
		print 'f('+str(k)+') =', a, '+-', e
		return self



class ResonanceRegion(Region):
	r"""
	region with resonances from :math:`2M_\pi` to specified upper bound.
	
	:math:`\pi` and :math:`\eta'` :math:`\delta`-contributions are *not* included here.
	
	Parameters
	----------
	k_max : float, optional
		upper bound of the region
	name : str, optional
		name of the region for plots
	init : bool, optional
		whether to update fits on instance creation
	"""
	def __init__(self, k_max=inf, name='resonance region', init=False):
		Region.__init__(self, ResonanceRegionFit(), k_min=2*gg.M_pi, k_max=k_max, name=name, init=init)
	
	def k_bounds(self):
		r"""
		get current bounds of the region.
		
		Here the lower bound is shifted as a pion mass and
		the upper bound is shifted according to :math:`\rho`-meson mass shift model for mesons.
		
		Returns
		-------
		tuple
			(k_min, k_max)
		"""
		return map(gg.s2k, (4*gg.M_pi**2, gg.shift_mass(self.k_max)**2))



class ReggeRegion(Region):
	r"""
	high-energy Regge behavior region from specified lower bound up to :math:`\infty`.
	
	:math:`\pi` and :math:`\eta'` :math:`\delta`-contributions are *not* included here.
	
	Parameters
	----------
	k_min : float, optional
		lower bound of the region
	name : str, optional
		name of the region for plots
	init : bool, optional
		whether to update fits on instance creation
	"""
	def __init__(self, k_min=1.5, name='high-energy region', init=True):
		Region.__init__(self, ReggeFit(), k_min=k_min, k_max=inf, name=name, init=init)
	
	def k_bounds(self):
		r"""
		get current bounds of the region.
		
		Here the lower bound is shifted according to :math:`\rho`-meson mass shift model for mesons.
		
		Returns
		-------
		tuple
			(k_min, k_max)
		"""
		return gg.s2k(gg.shift_mass(self.k_min)**2), inf



class WholeRegion(Region):
	r"""
	whole energies region.
	
	Note that :math:`\pi` and :math:`\eta'` :math:`\delta`-contributions *are included* here.
	
	Parameters
	----------
	name : str, optional
		name of the region for plots
	filename1 : str
		data file with points for resonances fitting
	filename2 : str
		data file with points for Regge fitting
	add_regge : bool, optional
		whether to include Regge contribution
	init : bool, optional
		whether to update fits on instance creation
	"""
	def __init__(self, name='whole region', filename1='data/pennington.dat', filename2='data/full_data.dat', add_regge=False, init=False):
		self.name = name
		m = 1.5
		self.rr = ResonanceRegion(init=False)
		self.add_regge = add_regge
		if add_regge:
			self.hr = ReggeRegion(init=False)
			self.region = Region(Smooth(self.rr.fit, self.hr.fit, m, a0=0.1), *self.k_bounds(), init=False)
		else:
			self.region = Region(self.rr.fit, *self.k_bounds(), init=False)
		self.fit = self.region.fit
		self.update_data(filename1=filename1, filename2=filename2)
		if init:
			self.update_fit()
		self.pi_0 = Pi0()
		self.eta_p = EtaPrime()
	
	def update_data(self, filename1='data/pennington.dat', filename2='data/full_data.dat'):
		self.rr.update_data(filename1)
		if self.add_regge:
			self.hr.update_data(filename2)
		self.region.update_data(filename2)
		return Region.update_data(self, filename2)
	
	def update_fit(self, quiet=False):
		self.rr.update_fit()
		if self.add_regge:
			self.hr.update_fit()
		return self
		#return Region.update_fit(self, quiet)
	
	def k_bounds(self):
		return gg.s2k(4*gg.M_pi**2), inf
	
	def integral_f(self, k):
		return map(sum, zip(Region.integral_f(self, k), self.pi_0.integral_f(k), self.eta_p.integral_f(k)))



class Pi0(Region):
	r"""
	single :math:`\pi^0` :math:`\delta`-contribution.
	
	Pseudo-region class to mimic part of normal region instance usage.
	Only `self.plot_f()` and `self.integral_f()` should be used.
	
	Parameters
	----------
	name : str, optional
		name to display on graphs
	"""
	def __init__(self, name='pion delta-contribution'):
		self.name=name
	
	def integral_f(self, k):
		r"""
		.. math:: \frac{\alpha^2}{\pi^2f_\pi^2}\; \frac{\nu^2 X_0}{\nu_0(\nu_0^2-\nu^2)}\; \left| \frac{F_{Q_1Q_2}}{F_{00}} \right|^2
		
		Parameters
		----------
		k : float
			:math:`\sqrt{2 q_1 \cdot q_2}, \; [GeV]`
		
		Returns
		-------
		tuple
			value of Re(f(k)-f(0)), 0.
		"""
		mon_m = 0.776 # monopole mass parameter
		nu = gg.k2nu(k)
		FF = 1./((1 + gg.Q1/mon_m**2)*(1 + gg.Q2/mon_m**2))
		nu0 = gg.s2nu(gg.M_pi**2) # note that we use same mass for pi+ and pi0
		X0 = gg.nu2X(nu0)
		return (alpha*nu*FF / (pi*gg.f_pi))**2 * X0 / (nu0*(nu0**2 - nu**2)), 0.



class EtaPrime(Region):
	r"""
	single :math:`\eta'` :math:`\delta`-contribution.
	
	Pseudo-region class to mimic part of normal region instance usage.
	Only `self.plot_f()` and `self.integral_f()` should be used.
	
	Parameters
	----------
	name : str, optional
		name to display on graphs
	"""
	def __init__(self, name='eta-prime delta-contribution'):
		self.name=name
	
	def integral_f(self, k):
		r"""
		.. math::
			64\pi\frac{\nu^2\Gamma_{\gamma\gamma} X_0}{\nu_0(\nu_0^2-\nu^2)m^3} \left| \frac{F_{Q_1Q_2}}{F_{00}} \right|^2
		
		Note that the mass is shifted and :math:`\Gamma_{\gamma\gamma}` is assumed to be proportional to :math:`m^3`,
		while :math:`F_{00}` assumed to remain constant.
		
		Parameters
		----------
		k : float
			:math:`\sqrt{2 q_1 \cdot q_2}, \; [GeV]`
		
		Returns
		-------
		tuple
			value of Re(f(k)-f(0)), 0.
		"""
		m = 0.95766 # physical mass
		mon_m = 0.859 # monopole mass parameter
		W_gg = 4.3*10**-6 # hardcoded physical value of gg-width
		nu0 = gg.s2nu(gg.shift_mass(m)**2)
		# note that there's no mistake in shifting the mass and width here
		X0 = gg.nu2X(nu0)
		nu = gg.k2nu(k)
		FF = 1./((1 + gg.Q1/mon_m**2)*(1 + gg.Q2/mon_m**2))
		return 64*pi * nu**2 * W_gg * X0 * FF**2 / (m**3 * nu0 * (nu0**2 - nu**2)), 0.


class TensorMesonRegion(Region):
	r"""
	single tensor region :math:`\delta`-contribution.
	
	Pseudo-region class to mimic part of normal region instance usage.
	Only `self.plot_f()` and `self.integral_f()` should be used.
	
	Parameters
	----------
	name : str, optional
		name to display on graphs
	"""
	def __init__(self, m=1.27, W_gg=0.5*(3.49+2.93)*10**-6, mon_m=inf, name='tensor meson delta-contribution'):
		self.name = name
		self.m = m
		self.W_gg = W_gg
		self.mon_m = mon_m
	
	def integral_f(self, k):
		r"""
		.. math::
			64\pi\frac{\nu^2\Gamma_{\gamma\gamma} X_0}{\nu_0(\nu_0^2-\nu^2)m^3} \left( \frac{20X_0}{m^4} + \frac{5\nu_0^2}{X_0} \right)
		
		Parameters
		----------
		k : float
			:math:`\sqrt{2 q_1 \cdot q_2}, \; [GeV]`
		
		Returns
		-------
		tuple
			value of Re(f(k)-f(0)), 0.
		"""
		m = gg.shift_mass(self.m)
		W_gg = self.W_gg * m/self.m
		nu0 = gg.s2nu(m**2)
		X0 = gg.nu2X(nu0)
		nu = gg.k2nu(k)
		FF = 1./((1 + gg.Q1/self.mon_m**2)*(1 + gg.Q2/self.mon_m**2))
		return 64*pi * nu**2 * W_gg * X0 * (20*X0/m**4 + 5*nu0**2/X0) / (m**3 * nu0 * (nu0**2 - nu**2)) * FF**2, 0.
