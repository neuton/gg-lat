r"""
a common module with constants and helper functions definitions for :math:`\gamma\gamma` scattering.

Here and everywhere else:

- :math:`k = \sqrt{2 \, q_1 \cdot q_2}`.
	When both photons are *real*, it's a total energy of 2 photons in CM (mass of the system) in [GeV].

- all energy values are in :math:`[GeV]`, while cross-sections are in :math:`[\mu b]`.

Masses and photons virtualities are set here to initial *physical* values.
This values should be changed from scripts which import this module like so:

.. code-block:: python

	import gg
	gg.Q1 = 1.  # change first photon virtuality
	print gg.k2s(1.) # do something after changes

------------------------------------------------------------------------
"""


hc = 0.197326972 # [fm*GeV]
# [hc^2] = [fm^2 * GeV^2] = [10^4 mcb * GeV^2]
r""" :math:`\hbar c, \, [fm \cdot GeV]`

Used for convertion from :math:`[\mu b]` to :math:`[1/GeV^2]`:

.. math:: [\hbar c^2] = [fm^2 \cdot GeV^2] = [10^4 \mu b \cdot GeV^2]
.. math:: 1 \; \frac{1}{GeV^2} = 10^4 \, (\hbar c \; [fm \cdot GeV])^2 \, \left[\frac{1}{GeV^2}\right]
"""

alpha = 0.00729735257
r" fine structure constant: :math:`\alpha` "


# Initial physical parameters values:

f_pi = 0.0924 # [GeV]
r" pion decay constant: :math:`f_\pi, \, [GeV]` "

M_pi  = 0.13957018 # [GeV]
r" pion mass for both :math:`\pi^0` and :math:`\pi^{+/-}`: :math:`M_{\pi}, \, [GeV]` "

M_rho = 0.7755 # [GeV]
r" :math:`\rho`-meson mass: :math:`M_{\rho}, \, [GeV]` "

Q1 = 0. # [GeV^2] : -q^2 of one photon
r" virtuality of first photon: :math:`Q_1^2, \, [GeV^2]` "

Q2 = 0. # [GeV^2] : -q^2 of another photon
r" virtuality of second photon: :math:`Q_2^2, \, [GeV^2]` "


def nu2s(nu):
	r"""
	convert :math:`\nu\, [GeV^2] \rightarrow s\, [GeV^2]`:
	
	.. math:: s = 2\nu - Q_1^2 - Q_2^2
	"""
	return 2*nu - Q1 - Q2


def s2nu(s):
	r"""
	convert :math:`s\, [GeV^2] \rightarrow \nu\, [GeV^2]`:
	
	.. math:: \nu = \frac{1}{2}(s + Q_1^2 + Q_2^2)
	"""
	return 0.5*(s + Q1 + Q2)


def k2nu(k):
	r"""
	convert :math:`k\, [GeV] \rightarrow \nu\, [GeV^2]`:
	
	.. math:: \nu = \frac{k^2}{2}
	"""
	return 0.5 * k*k


def nu2k(nu):
	r"""
	convert :math:`\nu\, [GeV^2] \rightarrow k\, [GeV]`:
	
	.. math:: k = \sqrt{2\nu}
	"""
	return (2*nu)**0.5


def k2s(k):
	r"""
	convert :math:`k\, [GeV] \rightarrow s\, [GeV^2]`:
	
	.. math:: s = k^2 - Q_1^2 - Q_2^2
	"""
	return k*k - Q1 - Q2


def s2k(s):
	r"""
	convert :math:`s\, [GeV^2] \rightarrow k\, [GeV]`:
	
	.. math:: k = \sqrt{s + Q_1^2 + Q_2^2}
	"""
	return (s + Q1 + Q2)**0.5


def nu2X(nu):
	r"""
	evaluate X(nu):
	
	.. math:: X = \nu^2 - Q_1^2 Q_2^2
	"""
	return nu*nu - Q1*Q2


def shift_mass(m):
	r"""
	move a meson mass according to :math:`\rho`-meson mass shift:
	
	.. math:: \Delta M = M_\rho^{(M_\pi=0.451\,\mathrm{GeV})} - M_\rho^{(M_\pi=0.140\,\mathrm{GeV})}
	"""
	return m + M_rho - 0.7755
	#return m + (0.952 - 0.7755)*(M_pi - 0.13957018)/(0.451 - 0.13957018)


def x2x(x):
	r"""
	dummy function: f(x) = x
	"""
	return x


def reverse_convert(f):
	r"""
	a helper function to get reversed convertion function.
	
	Parameters
	----------
	f : function
		one of convertion functions defined above
	
	Returns
	-------
	function
		reversed convertion function
	"""
	if f is s2k:
		return k2s
	if f is k2s:
		return s2k
	if f is s2nu:
		return nu2s
	if f is nu2s:
		return s2nu
	if f is nu2k:
		return k2nu
	if f is k2nu:
		return nu2k
	if f is x2x:
		return f
	raise Exception('Unknown convertion function')
