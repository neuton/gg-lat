r"""
contains classes for :math:`\gamma\gamma` total absorptive cross-section fitting procedure.

`Imports`: ``scipy``, ``numpy``, ``gg``, ``common``

Examples
--------
- See fits for real photons:

	``test.py``:
	
	.. code-block:: python

		import gg
		from common import *
		from fit import *
		xdata, ydata, errdata = read_xye('data/pennington.dat')
		r = ResonanceRegionFit().fit(xdata, ydata, errdata) # fitting to the data
		p = PiPiFit() # to plot with resonances
		rx = arange(0, 2., 0.01)
		ry = [r(x) for x in rx]
		py = [p(x) for x in rx]
		data = [rx, ry, py] + [[res(x) for x in rx] for res in r.resonances[:3]]
		plots = ['1 2 color=orange linewidth=2 label=fit', '1 3 label="$\pi^+ + \pi^-$ scalar QED"']
		plots += ['1 '+str(i+4) for i in range(3)]
		plot_to_file(data, plots, 'out/fit', xlogscale=False) # putting values to "out/fit" file

	Then invoke it with

	>>> python test.py

	And plot results with ``plot.py`` (where "out/fit" is a file with values just created):

	>>> python plot.py out/fit

- See how resonances fit extrapolates to higher pion masses:

	``test.py``:
	
	.. code-block:: python

		import gg
		from common import *
		from fit import *
		xdata, ydata, errdata = read_xye('data/pennington.dat')
		r = ResonanceRegionFit().fit(xdata, ydata, errdata)
		rx = arange(0, 2., 0.0002)
		for gg.f_pi, gg.M_pi, gg.M_rho in zip([0.0924, 0.109, 0.119], [0.1396, 0.324, 0.451], [0.7755, 0.922, 0.952]):
		# ^ note that you have to explicitely set gg.f_pi and gg.M_rho as well
			ry = [r(x) for x in rx]
			plot = r'1 2 label="$M_{\pi}: \, ' + str(gg.M_pi) + '$"'
			plot_to_file([rx, ry], [plot], 'out/fit_m_%.2f' % gg.M_pi)

	Then invoke it and plot results:

	>>> python test.py
	>>> python plot.py out/fit_m_0.14 out/fit_m_0.32 out/fit_m_0.45

- See how resonances fit extrapolates to higher Q-square:

	``test.py``:

	.. code-block:: python

		import gg
		from common import *
		from fit import *
		xdata, ydata, errdata = read_xye('data/pennington.dat')
		r = ResonanceRegionFit().fit(xdata, ydata, errdata)
		rx = arange(0, 2., 0.0002)
		gg.f_pi, gg.M_pi, gg.M_rho = 0.119, 0.451, 0.952
		for Q in [0., 0.2, 0.4, 0.8]:
			gg.Q1 = Q
			ry = [r(x) for x in rx]
			plot = r'1 2 label="$Q_1^2: \, %.1f$"' % Q
			plot_to_file([rx, ry], [plot], 'out/fit_q_%.1f' % Q)

	Then invoke it and plot results:

	>>> python test.py
	>>> python plot.py out/fit_q_0.0 out/fit_q_0.2 out/fit_q_0.4 out/fit_q_0.8

------------------------------------------------------------------------
"""


import numpy
from scipy.optimize import curve_fit, leastsq
from scipy import log, exp, pi, inf

from gg import alpha, hc
import gg
from common import *



def chisqr(f, xdata, ydata, errdata):
	r"""
	compute a total :math:`\chi`-square value for a given X-data, Y-data and error-data:
	
	.. math:: \chi^2 = \sum_x \frac{(y_x - f(x))^2}{(\Delta y_x)^2}
	"""
	return sum([((y-f(x))/err)**2 for x, y, err in zip(xdata, ydata, errdata)])



class FunctionFit:
	"""
	a base class for fits.
	
	When inherited, `f(self, k, *p)` method should be overriden with a desired fitting function.
	
	`__call__(self, k, *p)` method is defined here, which simply invokes `f(self, k, *p)`,
	so that the class instance may be thought of as a simple function.
	
	Parameters
	----------
	f : function, optional
		the fitting function
	p : list, optional
		list of initial values for fitting function parameters
	"""
	def __init__(self, f=None, p=[]):
		if f is not None:
			self.f = f
		self.set_p(*p)
	
	def __call__(self, k, *p):
		p = self.full_p(p)
		if type(k) == numpy.ndarray:
			return numpy.array([self.f(x, *p) for x in k], dtype=float)
		else:
			return self.f(k, *p)
	
	def f(k, *p):
		"""
		the fitting function (to be overriden when inherited)
		
		Parameters
		----------
		k : float
			:math:`\sqrt{2 q_1 \cdot q_2}` in GeV
		p : args, optional
			explicitely specified list of fitting parameters
		
		Returns
		-------
		float
			fit function value at specified `k` with (optionally) specified parameters
		"""
		return 0.
	
	def get_p(self):
		"""
		Returns
		-------
		list
			current fit parameters list
		"""
		return self.p
	
	def set_p(self, *p):
		"""
		explicitely set fit parameters.
		
		Parameters
		----------
		p : args
			list of fit parameters to set
		"""
		#self.p = self.full_p(p)
		self.p = p
	
	def full_p(self, p):
		"""
		Parameters
		----------
		p : list
			list of parameters (len(p) may be < len(self.p))
		
		Returns
		-------
		list
			full parameters list with first len(p) parameters changed to `p`
		"""
		return list(p) + list(self.get_p())[len(p):]
	
	def get_errors(self, xdata, ydata, errdata):
		"""
		Parameters
		----------
		xdata : list
			list of X-values of data
		ydata : list
			list of Y-values of data
		errdata : list
			list of error-values of data
		
		Returns
		-------
		list
			list of errors estimation (based on weightening with :math:`\chi`-square) for each point in data
		"""
		errs = [(y-self(x))**2/abs(err) for x,y,err in zip(xdata,ydata,errdata)]
		#errs = [abs(y-self(x)) for x,y,err in zip(xdata,ydata,errdata)]
		#print 'max fit error:', max(errs), ';  min fit error:', min(errs)
		return errs
	
	def fit(self, xdata, ydata, errdata, quiet=False):
		"""
		fit the fitting function to the data.
		
		uses `curve_fit` and then `leastsq` from ``scipy.optimize``
		
		Parameters
		----------
		xdata : list
			list of X-values of data
		ydata : list
			list of Y-values of data
		errdata : list
			list of error-values of data
		quiet : bool, optional
			stay quiet
		
		Returns
		-------
		FunctionFit
			self
		"""
		if not quiet: print 'fitting ...'
		p, pcov = curve_fit(self, xdata, ydata, p0=self.get_p())
		if not quiet:
			print 'chisqr =', chisqr(self, xdata, ydata, errdata)
			print 'curve_fit: success'
		r = lambda p: numpy.array([abs(y - self(x, *p))/err for x,y,err in zip(xdata,ydata,errdata)], dtype=float)
		p, pcov = leastsq(func=r, x0=p, full_output=True)[:2]
		self.set_p(*p)
		self.chisqr = chisqr(self, xdata, ydata, errdata)
		if not quiet:
			print 'chisqr =', self.chisqr
			print 'leastsq: success'
		return self



class Smooth(FunctionFit):
	"""
	smooth transition function between 2 given functions.
	
	Parameters
	----------
	f1 : function
		left function
	f2 : function
		right function
	m : float
		the transition point
	a0 : float
		initial width of the transition
	"""
	def __init__(self, f1, f2, m, a0=0.1):
		self._f1, self._f2, self._m = f1, f2, m
		FunctionFit.__init__(self, p=[a0])
	
	def f(self, k, a):
		r"""
		function of a form:
		
		.. math:: f_1(k) \, \frac{1}{1 + \mathrm{exp}(\frac{s - m^2}{a})} \; + \;
			f_2(k) \, \frac{\mathrm{exp}(\frac{s - m^2}{a})}{1 + \mathrm{exp}(\frac{s - m^2}{a})}
		
		Note that `m` is shifted according to meson mass shift.
		
		Parameters
		----------
		k : float
			:math:`\sqrt{2 q_1 \cdot q_2}, \; [GeV]`
		a : float
			width of the transition (a single parameter of the fit)
		
		Returns
		-------
		float
			function value at `k` and with specified `a`
		"""
		numpy.seterr(all='ignore') # ignore the warning about overflow in exp and underflow in double_scalars
		e = 1./(1. + exp((gg.k2s(k) - gg.shift_mass(self._m)**2)/a))
		e1 = 1. - e
		numpy.seterr(all='warn')
		if e == 0 or e1 == 1:
			return self._f2(k)
		if e == 1 or e1 == 0:
			return self._f1(k)
		r = e * self._f1(k) + e1 * self._f2(k)
		return r



class ResonanceFit(FunctionFit):
	r"""
	a base fit for a resonance.
	
	Parameters
	----------
	J : int
		total spin :math:`J`
	m : float
		*physical* mass of the resonance (:math:`m`), [GeV] --- initial fitting parameter
	W_tot : float
		*physical* width of the resonance (:math:`\Gamma_{tot}`), [GeV] --- initial fitting parameter
	W_gg : float
		*physical* :math:`\gamma\gamma` decay width of the resonance (:math:`\Gamma_{\gamma\gamma}`), [GeV] --- initial fitting parameter
	"""
	def __init__(self, J, m, W_tot, W_gg):
		FunctionFit.__init__(self, p=[m, W_tot, W_gg])
		self.J = J
	
	def fall_off(self, k, *p):
		r"""
		non-model empirical fall-off to 0 function of a form:
		
		.. math:: f(s) = 0 \;\;\;\;\;\; \mathrm{for} \;\; s \le 4M_{\pi}^2 \; ,
		.. math:: f(s) = 1 - 2 \cdot exp\left(-\frac{\ln(2)}{\left(1 - \frac{2(s - 4M_{\pi}^2)}{m^2 - 4M_{\pi}^2}\right)}\right) \;\;\;
			\mathrm{for} \;\; 4M_{\pi}^2 < s < (m^2 + 4M_{\pi}^2)/2 \; ,
		.. math:: f(s) = 1 \;\;\;\;\;\; \mathrm{for} \;\; s \ge (m^2 + 4M_{\pi}^2)/2 \; ,
		
		To be multiplied by a Breit-Wigner resonance, so that it will be zero below 2-pion production threshold.
		
		Parameters
		----------
		k : float
			:math:`\sqrt{2 q_1 \cdot q_2}, \; [GeV]`
		p : args, optional
			explicitely specified list of fitting parameters: `m` (physical :math:`m`) and `W_tot` (physical :math:`\Gamma_{tot}`)
		
		Returns
		-------
		float
			function value at specified `k`, in :math:`[\mu b / GeV]`
		"""
		m = gg.shift_mass(self.full_p(p)[0])
		s = gg.k2s(k)
		s_0 = 4*gg.M_pi**2
		s_1 = (m**2 + s_0)*0.5
		if s <= s_0:
			return 0.
		elif s >= s_1:
			return 1.
		else:
			numpy.seterr(all='ignore') # ignore the warning about underflow in exp
			r = 1 - 2*exp(- log(2) / abs(1 - (s - s_0)/(s_1 - s_0)))
			numpy.seterr(all='warn')
			return r
	
	def breit_wigner(self, k, *p):
		r"""
		a function of a relativistic Breit-Wigner form (multiplied by some convinience factor):
		
		.. math:: \frac{1}{2} 16\pi \frac{\Gamma_{tot}}{(s - m^2)^2 + m^2\Gamma_{tot}^2} \;\; [\mu b / GeV]
		
		Parameters
		----------
		k : float
			:math:`\sqrt{2 q_1 \cdot q_2}, \; [GeV]`
		p : args, optional
			explicitely specified list of fitting parameters: `m` (physical :math:`m`) and `W_tot` (physical :math:`\Gamma_{tot}`)
		
		Returns
		-------
		float
			function value at specified `k`, in :math:`[\mu b / GeV]`
		"""
		m, w = self.full_p(p)[:2]
		m = gg.shift_mass(m)
		s = gg.k2s(k)
		return 10**4 * hc**2 * 8 * pi * w / ((s - m*m)**2 + w*w * m*m)
	
	def form_factor_fraction(self):
		r"""
		default form-factor fraction:
		
		.. math:: \left| \frac{F(Q_1, Q_2)}{F(0, 0)} \right| = 1
		
		Returns
		-------
		float
			form-factor fraction
		"""
		return 1.
	
	def shift_W_gg(self, *p):
		r"""
		shift :math:`\Gamma_{\gamma\gamma}` *linearly* with the pion mass shift.
		
		Parameters
		----------
		p : args, optional
			`m` (physical :math:`m`), `W_tot` (physical :math:`\Gamma_{tot}`) and `W_gg` (physical :math:`\Gamma_{gg}`)
		
		Returns
		-------
		float
			shifted :math:`\Gamma_{\gamma\gamma}`
		"""
		m, _, W_gg = self.full_p(p)
		return W_gg * gg.shift_mass(m)/m
	
	def k_factor(self, *p):
		r"""
		default kinematical factor:
		
		.. math:: \frac{k_0^2}{m^2}
		
		Parameters
		----------
		p : args, optional
			`m` (physical :math:`m`)
		
		Returns
		-------
		float
			kinematical factor
		"""
		m = gg.shift_mass(self.full_p(p)[0])
		k0 = gg.s2k(m**2)
		return k0**2 / m**2
	
	def f(self, k, *p):
		r"""
		a fit function for a meson resonance of a general form:
		
		.. math::
			\mathrm{fall\_off\_function}(k) \cdot \;\;
			\frac{1}{2} 16\pi \;
			\frac{\Gamma_{tot}}{(s - m^2)^2 + m^2\Gamma_{tot}^2} \;
			(2J + 1) \Gamma_{\gamma\gamma} \; \cdot \; \mathrm{k\_factor} \; \cdot \;
			\left| \frac{F(Q_1, Q_2)}{F(0, 0)} \right|^2 \;\;\;\;\; [\mu b]
		
		Parameters
		----------
		k : float
			:math:`\sqrt{2 q_1 \cdot q_2}, \; [GeV]`
		p : args, optional
			explicitely specified list of fitting parameters: `m` (physical :math:`m`), `W_tot` (physical :math:`\Gamma_{tot}`), `W_gg` (physical :math:`\Gamma_{gg}`) and `J` (spin :math:`J`)
		
		Returns
		-------
		float
			fit function value at specified `k` with (optionally) specified parameters, in :math:`[\mu b]`
		"""
		p = self.full_p(p)
		return self.fall_off(k, *p) * (2*self.J + 1) * self.breit_wigner(k, *p) * self.shift_W_gg(*p) * self.k_factor(*p) * self.form_factor_fraction()**2
	
	def get_p(self):
		return self.m, self.W_tot, self.W_gg
	
	def set_p(self, m, W_tot, W_gg):
		self.m, self.W_tot, self.W_gg = m, W_tot, W_gg



class ScalarResonanceFit(ResonanceFit):
	r"""
	a base fit for scalar meson resonance.
	
	Parameters
	----------
	m : float
		*physical* mass of the resonance (:math:`m`), [GeV] --- initial fitting parameter
	W_tot : float
		*physical* width of the resonance (:math:`\Gamma_{tot}`), [GeV] --- initial fitting parameter
	W_gg : float
		*physical* :math:`\gamma\gamma` decay width of the resonance (:math:`\Gamma_{\gamma\gamma}`), [GeV] --- initial fitting parameter
	mon_m : float, optional
		monopole mass for the form factor fraction, [GeV]
	"""
	def __init__(self, m, W_tot, W_gg, mon_m=inf):
		ResonanceFit.__init__(self, 0, m, W_tot, W_gg)
		self.mon_m = mon_m
	
	def k_factor(self, *p):
		r"""
		kinematical factor.
		
		.. math:: \frac{2 \nu_0^2}{m^2 \sqrt{X_0}}
		
		Parameters
		----------
		p : args, optional
			`m` (physical :math:`m`)
		
		Returns
		-------
		float
			kinematical factor
		"""
		m = gg.shift_mass(self.full_p(p)[0])
		nu0 = gg.s2nu(m**2)
		X0 = gg.nu2X(nu0)
		return 2 * nu0**2 / (m**2 * X0**0.5)
	
	def form_factor_fraction(self):
		r"""
		monopole-form form-factor fraction:
		
		.. math:: \left| \frac{F(Q_1, Q_2)}{F(0, 0)} \right| = \frac{1}{1 + Q_1^2/\Lambda^2} \cdot \frac{1}{1 + Q_2^2/\Lambda^2}
		
		Returns
		-------
		float
			form-factor fraction
		"""
		return 1. / (1 + gg.Q1/self.mon_m**2) / (1 + gg.Q2/self.mon_m**2)



class AxialVectorResonanceFit(ResonanceFit):
	r"""
	a base fit for axial vector meson resonance.
	
	Parameters
	----------
	m : float
		*physical* mass of the resonance (:math:`m`), [GeV] --- initial fitting parameter
	W_tot : float
		*physical* width of the resonance (:math:`\Gamma_{tot}`), [GeV] --- initial fitting parameter
	W_gg : float
		*physical* :math:`\gamma\gamma` decay width of the resonance (:math:`\Gamma_{\gamma\gamma}`), [GeV] --- initial fitting parameter
	dip_m : float, optional
		dipole mass for the form factor fraction, [GeV]
	"""
	def __init__(self, m, W_tot, W_gg, dip_m=inf):
		ResonanceFit.__init__(self, 0, m, W_tot, W_gg)
		self.dip_m = dip_m
	
	def k_factor(self, *p):
		r"""
		kinematical factor.
		
		.. math:: \frac{(Q_1^2 - Q_2^2)^2}{m^4} \; \frac{2 \nu_0^2}{m^2 \sqrt{X_0}}
		
		Parameters
		----------
		p : args, optional
			`m` (physical :math:`m`)
		
		Returns
		-------
		float
			kinematical factor
		"""
		m = gg.shift_mass(self.full_p(p)[0])
		nu0 = gg.s2nu(m**2)
		X0 = gg.nu2X(nu0)
		return (gg.Q1 - gg.Q2)**2 * 2 * nu0**2 / (m**6 * X0**0.5)
	
	def form_factor_fraction(self):
		r"""
		dipole-form form-factor fraction:
		
		.. math:: \left| \frac{F(Q_1, Q_2)}{F(0, 0)} \right| = \left( \frac{1}{1 + Q_1^2/\Lambda^2} \cdot \frac{1}{1 + Q_2^2/\Lambda^2} \right)^2
		
		Returns
		-------
		float
			form-factor fraction
		"""
		l2 = self.dip_m**2
		return ( 1. / (1 + gg.Q1/l2) / (1 + gg.Q2/l2) )**2



class TensorResonanceFit(ResonanceFit):
	r"""
	a base fit for a tensor meson resonance.
	
	Parameters
	----------
	m : float
		*physical* mass of the resonance (:math:`m`), [GeV] --- initial fitting parameter
	W_tot : float
		*physical* width of the resonance (:math:`\Gamma_{tot}`), [GeV] --- initial fitting parameter
	W_gg : float
		*physical* :math:`\gamma\gamma` decay width of the resonance (:math:`\Gamma_{\gamma\gamma}`), [GeV] --- initial fitting parameter
	"""
	def __init__(self, m, W_tot, W_gg):
		ResonanceFit.__init__(self, 0, m, W_tot, W_gg)
	
	def k_factor(self, *p):
		r"""
		kinematical factor.
		
		.. math:: \frac{2 \nu_0^2}{m^2 \sqrt{X_0}} \; + \; \frac{8X_0\sqrt{X_0}}{m^6}
		
		Note that 2 factors coming from :math:`\sigma_0` and :math:`\sigma_2` are added here.
		:math:`\Gamma_{\gamma\gamma}` widths are assumed to be equal for both :math:`\sigma_0` and :math:`\sigma_2`.
		
		Parameters
		----------
		p : args, optional
			`m` (physical :math:`m`)
		
		Returns
		-------
		float
			kinematical factor
		"""
		m = gg.shift_mass(self.full_p(p)[0])
		nu0 = gg.s2nu(m**2)
		X0 = gg.nu2X(nu0)
		return 2 * nu0**2 / (m**2 * X0**0.5) + 8*X0**1.5/m**6



# scalar resonances:

class f0_500(ScalarResonanceFit):
	def __init__(self, m=0.5, W_tot=0.5, W_gg=2.05*10**-6):
		ScalarResonanceFit.__init__(self, m, W_tot, W_gg)

class f0_980(ScalarResonanceFit):
	def __init__(self, m=0.998, W_tot=0.042, W_gg=0.32*10**-6):
		ScalarResonanceFit.__init__(self, m, W_tot, W_gg)

class f0_1370(ScalarResonanceFit):
	def __init__(self, m=1.44, W_tot=0.4, W_gg=4.*10**-6):
		ScalarResonanceFit.__init__(self, m, W_tot, W_gg)


# axial vector resonances:

class f1_1285(AxialVectorResonanceFit):
	def __init__(self, m=1.2819, W_tot=24.2*10**-3, W_gg=3.5*10**-6):
		AxialVectorResonanceFit.__init__(self, m, W_tot, W_gg, dip_m=1.040)

class f1_1420(AxialVectorResonanceFit):
	def __init__(self, m=1.4264, W_tot=54.9*10**-3, W_gg=3.2*10**-6):
		AxialVectorResonanceFit.__init__(self, m, W_tot, W_gg, dip_m=0.926)


# tensor resonances:

class a2_1320(TensorResonanceFit):
	def __init__(self, m=1.313, W_tot=0.107, W_gg=1.04*10**-6):
		TensorResonanceFit.__init__(self, m, W_tot, W_gg)

class f2_1270(TensorResonanceFit):
	def __init__(self, m=1.27, W_tot=0.1851, W_gg=0.5 * (3.49 + 2.93)*10**-6):
		TensorResonanceFit.__init__(self, m, W_tot, W_gg)



class ResonanceRegionFit(FunctionFit):
	r"""
	a compound fit for all included resonances + background for the resonance region.
	
	Currently includes:
	
	.. math:: f_0(980), \; a_2(1320), \; f_2(1270), \; f_1(1285), \; f_1(1420)
	"""
	def __init__(self):
		self.resonances = [f0_980(), a2_1320(), f2_1270()] + [f1_1285(), f1_1420()]
		# setting already fitted parameters here, so that we don't need to refit each time:
		p0 = [0.98245702190502038, 0.034460567775547173, 2.1173056323462041e-07,
			1.2625605391874659, 0.13163851473973218, 3.8906514125532965e-06,
			1.1782303951657658, 0.15573290936433748, 1.787955427585119e-06]
		self.background = PiPiFit()
		FunctionFit.__init__(self, p=p0)#self.get_p())
	
	def f(self, k, *p):
		r"""
		a sum of all included resonances + background of a form of scalar QED :math:`\pi^+\pi^-` decay cross-section.
		
		Parameters
		----------
		k : float
			:math:`\sqrt{2 q_1 \cdot q_2}, \; [GeV]`
		p : args, optional
			explicitely specified list of fitting parameters --- sequence of fitting parameters of included resonances
		
		Returns
		-------
		float
			fit function value at specified `k` with (optionally) specified parameters
		"""
		self.set_p(*self.full_p(p))
		return sum([r(k) for r in self.resonances]) + self.background(k)
	
	def get_p(self):
		p = ()
		for r in self.resonances[:3]:
			p += r.get_p()
		return p
	
	def set_p(self, *p):
		for i in range(3):
			self.resonances[i].set_p(*p[3*i:3*i+3])



class ReggeFit(FunctionFit):
	r"""
	a fit of a Regge form for the high energy cross-section behaviour.
	"""
	def __init__(self):
		p0 = [0.0046331113535695153, 0.42660207416760498, 0.45493098451754294, -0.083114144853592309]
		FunctionFit.__init__(self, p=p0)#[0., 0., 0.093, -0.358])
	
	def f(self, k, *p):
		r"""
		Regge fitting function:
		
		.. math:: C_1 \cdot (2\nu)^a + C_2 \cdot (2\nu)^b
		
		Parameters
		----------
		k : float
			:math:`\sqrt{2 q_1 \cdot q_2}, \; [GeV]`
		p : args, optional
			explicitely specified list of fitting parameters --- :math:`C_1`, :math:`C_2`, :math:`a`, :math:`b`
		
		Returns
		-------
		float
			fit function value at specified `k` with (optionally) specified parameters
		"""
		#k = -gg.shift_mass(-k)
		#s = gg.k2s(k)
		nu = gg.k2nu(k)
		return p[0] * (2*nu)**p[2] + p[1] * (2*nu)**p[3]



class PiPiFit(FunctionFit):
	r"""
	a fit for the lowest order :math:`\pi^+\pi^-` decay cross-section in scalar QED.
	
	Without any fitting parameters.
	
	Also note that form-factors *are not included here*.
	"""
	def __init__(self):
		FunctionFit.__init__(self)
		#self.mon_m = 0.776 # GeV
	
	#def form_factor_fraction(self):
	#	r"""
	#	squared monopole-form form-factor fraction:
	#	
	#	.. math:: \left| \frac{F(Q_1, Q_2)}{F(0, 0)} \right|^2 = \left( \frac{1}{1 + Q_1^2/\Lambda^2} \cdot \frac{1}{1 + Q_2^2/\Lambda^2} \right)^2
	#	
	#	Returns
	#	-------
	#	float
	#		squared form-factor fraction
	#	"""
	#	return ( 1. / (1 + gg.Q1/self.mon_m**2) / (1 + gg.Q2/self.mon_m**2) )**2
	
	def f(self, k):
		r"""
		.. math:: \frac{1}{2} \alpha^2 \frac{\pi}{2} \frac{s^2\nu^2}{X^3} \left(
			\sqrt{a} \left(2 - a - \left(1 - \frac{2X}{s\nu}\right)^2\right) -
			(1-a)\left(3 - \frac{4X}{s\nu} + a\right)L \right),
		.. math:: L = \ln\left(\frac{1+\sqrt{a}}{\sqrt{1-a}}\right), \;\;\;\;\;\;\;
			a = \frac{X}{\nu^2} \left(1 - \frac{4m^2}{s}\right), \;\;\;\;\;\;\;
			X = \nu^2 - Q_1^2 \cdot Q_2^2
		"""
		Q1 = gg.Q1
		Q2 = gg.Q2
		nu = gg.k2nu(k)
		X = gg.nu2X(nu)
		s = gg.k2s(k)
		s2 = s*s
		m = gg.M_pi
		m2 = m*m
		if s <= 4*m2:
			return 0.
		else:
			a = X*(1 - 4*m2/s)/nu**2
			L = log((1 + a**0.5)/(1 - a)**0.5)
			r = (a**0.5*(2 - a - (1 - 2*X/(s*nu))**2) - (1-a)*(3 + a - 4*X/(nu*s))*L) * s2*nu**3/X**3
			return r * (10**4*hc**2*alpha**2*pi/4)# * self.form_factor_fraction()
