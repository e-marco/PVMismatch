# -*- coding: utf-8 -*-

"""
This module contains the :class:`~pvmismatch.pvmismatch_lib.pvcell.PVcell`
object which is used by modules, strings and systems.
"""

from pvmismatch.pvmismatch_lib.pvconstants import PVconstants
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import newton,minimize
import functools
import copy

# Defaults
VOLTAGES = [0,0.3781,0.3918,0.4056,0.4193,0.4331,0.4469,0.4606,0.4744,0.4881,0.5019,
            0.5157,0.5295,0.5433,0.5571,0.571,0.5849,0.5879,0.5889,0.589,0.599,0.6083,
            0.6093,0.6172,0.6235,0.6288,0.6334,0.6375,0.6412,0.6446,0.6477,0.6505,
            0.6532,0.6558,0.6582,0.6604,0.6626,0.6647,0.6667,0.6686,0.6704,0.6722,
            0.6739,0.6756,0.6772,0.6788,0.6803,0.6818,0.6833,0.6836,0.6847,0.6861,
            0.6871,0.6875,0.6887,0.6888,0.7299]
CURRENTS = [11.25459384,11.2264,11.2254,11.2242,11.223,11.2215,11.2197,11.2174,11.214,
            11.2089,11.2007,11.1874,11.1655,11.1289,11.0676,10.9647,10.792,10.742,
            10.7238,10.7219,10.5032,10.2184,10.1801,9.8494,9.5149,9.1783,8.8406,
            8.5021,8.1631,7.8238,7.4843,7.1446,6.8048,6.4649,6.1249,5.7849,5.4448,
            5.1047,4.7646,4.4244,4.0842,3.744,3.4038,3.0636,2.7234,2.3831,2.0429,
            1.7026,1.3623,1.2977,1.022,0.6817,0.4312,0.3414,0.0506,0.0065,-13.5725]
TCELL = 298.15  # [K] cell temperature
ARBD = 1.036748445065697E-4  # reverse breakdown coefficient 1
BRBD = 0.  # reverse breakdown coefficient 2
VRBD_ = -5.527260068445654  # [V] reverse breakdown voltage
NRBD = 3.284628553041425  # reverse breakdown exponent
EG = 1.1  # [eV] band gap of cSi
ALPHA_ISC = 0.0004# [1/K] short circuit current temperature coefficient
BETA_VOC = -0.0025 # [1/K] open-circuit voltage temperature coefficient
GAMMA_PMP = -0.0034 # [1/K] Pmpp temperature coefficient
EPS = np.finfo(np.float64).eps

def cached(f):
    """
    Memoize an object's method using the _cache dictionary on the object.
    """
    @functools.wraps(f)
    def wrapper(self):
        # note:  we use self.__wrapped__ instead of just using f directly
        # so that we can spy on the original function in the test suite.
        key = wrapper.__wrapped__.__name__
        if key in self._cache:
            return self._cache[key]
        value = wrapper.__wrapped__(self)
        self._cache[key] = value
        return value
    # store the original function to be accessible by the test suite.
    # functools.wraps already sets this in python 3.2+, but for older versions:
    wrapper.__wrapped__ = f
    return wrapper

class PVcell(object):
    """
    Class for PV cells.

    :param Voltages: list of cell voltage at STC [V]
    :param Currents: list of cell current at STC [A]
    :param aRBD: reverse breakdown coefficient 1
    :param bRBD: reverse breakdown coefficient 2
    :param VRBD: reverse breakdown voltage [V]
    :param nRBD: reverse breakdown exponent
    :param Eg: band gap [eV]
    :param alpha_Isc: short circuit current temp coeff [1/K]
    :param Tcell: cell temperature [K]
    :param Ee: incident effective irradiance [suns]
    :param pvconst: configuration constants object
    :type pvconst: :class:`~pvmismatch.pvmismatch_lib.pvconstants.PVconstants`
    """

    _calc_now = False  #: if True ``calcCells()`` is called in ``__setattr__``

    def __init__(self, Voltages=VOLTAGES, Currents=CURRENTS, aRBD=ARBD, bRBD=BRBD, 
                 VRBD=VRBD_, nRBD=NRBD, Eg=EG, alpha_Isc=ALPHA_ISC, beta_Voc=BETA_VOC, 
                 gamma_Pmp=None, irradiance_correction_factor=0.06, 
                 temperature_coefficient_internal_series_resistance=5e-5,
                 internal_series_resistance=None, gamma=3.0,
                 Tcell=TCELL, Ee=1., pvconst=PVconstants()):
        # set up property cache
        self._cache = {}
        # user inputs
        self.Voltages = Voltages
        self.Currents = Currents

        if not (len(self.Voltages) == len(self.Currents)):
            raise ValueError("Current and Voltage arrays must be same length.")
        
        self.aRBD = aRBD  #: reverse breakdown coefficient 1
        self.bRBD = bRBD  #: reverse breakdown coefficient 2
        self.VRBD = VRBD  #: [V] reverse breakdown voltage
        self.nRBD = nRBD  #: reverse breakdown exponent
        
        self.Eg = Eg  #: [eV] band gap of cSi
        self.alpha_Isc = alpha_Isc  #: [1/K] short circuit temp. coeff.
        self.beta_Voc_STC = beta_Voc
        self.gamma_Pmp = gamma_Pmp
        self.Tcell = Tcell  #: [K] cell temperature
        self.Ee = Ee  #: [suns] incident effective irradiance on cell
        
        self.pvconst = pvconst  #: configuration constants
        self.Icell = None  #: cell currents on IV curve [A]
        self.Vcell = None  #: cell voltages on IV curve [V]
        self.Pcell = None  #: cell power on IV curve [W]
        
        self.irradiance_correction_factor = irradiance_correction_factor
        self.internal_series_resistance = internal_series_resistance
        self.gamma = gamma

        if internal_series_resistance is None:
            # determine internal series resistance based on slope at Voc
            self.Ee = 1.0
            self.internal_series_resistance = \
                self._estimate_internal_series_resistance()
            self.Ee = Ee
            #print("Computed internal_series_resistance "+
            #      "{:0.05f}".format(self.internal_series_resistance))
        
        self.VocSTC = self._estimate_Voc()  #: estimated Voc at STC [V]
        self.Isc0_T0 = self._estimate_Isc0_T0()
        
        if gamma_Pmp is not None: # use target gamma_Pmp to fit k
            # find k (temperature curve correction factor)
            self.Ee = 1.0
            res = minimize(self.fit_function_k, x0=100, method='Nelder-Mead')
            res = minimize(self.fit_function_k, x0=res.x, method='Powell')
            res = minimize(self.fit_function_k, x0=res.x, method='Nelder-Mead')
            self.temperature_coefficient_internal_series_resistance = res.x / 1000000.
            #print("Computed temperature_coefficient_internal_series_resistance "+
            #      "{:0.05f}".format(self.temperature_coefficient_internal_series_resistance))
            self.Tcell = Tcell
            self.Ee = Ee
        else:
            self.temperature_coefficient_internal_series_resistance = \
                temperature_coefficient_internal_series_resistance

        
        # set calculation flag
        super(PVcell, self).__setattr__('_calc_now', True)
        self._calc_now = True  # overwrites the class attribute

    def __str__(self):
        fmt = '<PVcell(Ee=%g[suns], Tcell=%g[K], Isc=%g[A], Voc=%g[V])>'
        return fmt % (self.Ee, self.Tcell, self.Isc, self.Voc)

    def __repr__(self):
        return str(self)

    def __setattr__(self, key, value):
        # check for floats
        try:
            value = np.float64(value)
        except (TypeError, ValueError):
            pass  # fail silently if not float, eg: pvconst or _calc_now
        super(PVcell, self).__setattr__(key, value)
        # recalculate IV curve
        self._cache.clear()
                           
        if self._calc_now:
            Icell, Vcell, Pcell = self.calcCell()
            self.__dict__.update(Icell=Icell, Vcell=Vcell, Pcell=Pcell)

    def clone(self):
        """
        Return a copy of this object with the same pvconst.
        """
        cloned = copy.copy(self)
        super(PVcell, cloned).__setattr__('_cache', self._cache.copy())
        return cloned
    
    def update(self, **kwargs):
        """
        Update user-defined constants.
        """
        # turn off calculation flag until all attributes are updated
        self._calc_now = False
        # don't use __dict__.update() instead use setattr() to go through
        # custom __setattr__() so that numbers are cast to floats
        for k, v in kwargs.items():
            setattr(self, k, v)
        self._calc_now = True  # recalculate

    @property
    @cached
    def Vt(self):
        """
        Thermal voltage in volts.
        """
        return self.pvconst.k * self.Tcell / self.pvconst.q

    @property
    @cached       
    def Isc(self):
        return self.Ee * self.Isc0


    @property
    @cached       
    def Isc0(self):
        """
        Short circuit current at Tcell in amps.
        """
        _delta_T = self.Tcell - self.pvconst.T0  # [K] temperature difference
        return self.Isc0_T0 * (1. + self.alpha_Isc * _delta_T)  # [A] Isc0

    @property
    @cached       
    def Voc(self):
        """
        Estimate open circuit voltage of cells.
        Returns Voc : numpy.ndarray of float, estimated open circuit voltage
        """
        _delta_T = self.Tcell - self.pvconst.T0  # [K] temperature difference

        a1 = self.VocSTC * (1. + self.beta_Voc() * _delta_T)
        #a2 = self.VocSTC + self.beta_Voc(return_absolute=True) * _delta_T


        return a1 * (1. + self.irradiance_correction_factor * np.log(self.Ee))


    def Voc2(self):
        """
        Estimate open circuit voltage of cells.
        Returns Voc : numpy.ndarray of float, estimated open circuit voltage
        """
        _delta_T = self.Tcell - self.pvconst.T0  # [K] temperature difference

        #a1 = self.VocSTC * (1. + self.beta_Voc() * _delta_T)
        a2 = self.VocSTC * (1. + self.irradiance_correction_factor * np.log(self.Ee)) \
            + self.beta_Voc(return_absolute=True) * _delta_T


        return a2 

    # def beta_Voc(self, Vin, return_absolute=False):
    #     Eg0 = 1.206
        
    #     TC_Voc_abs = (Vin / self.Tcell) - (Eg0 / self.Tcell) - \
    #         self.gamma * self.Vt / self.Tcell
        
    #     if return_absolute:
    #         return TC_Voc_abs
        
    #     TC_Voc_rel = TC_Voc_abs / Vin
    #     return TC_Voc_rel
    

    def beta_Voc(self, return_absolute=False):
        Eg0 = 1.206
        
        TC_Voc_abs = (self.VocSTC * \
                      (1. + self.irradiance_correction_factor * np.log(self.Ee)) / \
                        self.Tcell) - (Eg0 / self.Tcell) - \
            self.gamma * self.Vt / self.Tcell
        
        if return_absolute:
            return TC_Voc_abs
        
        TC_Voc_rel = TC_Voc_abs / \
            (self.VocSTC * (1. + self.irradiance_correction_factor * np.log(self.Ee)))
        return TC_Voc_rel


    def _estimate_internal_series_resistance(self):
        select_val = (abs(np.asarray(self.Currents)) <= 
                      0.1*max(np.asarray(self.Currents)))

        X = np.asarray(self.Currents)[select_val]
        Y = np.asarray(self.Voltages)[select_val]
        slope = -1 * np.polyfit(X,Y,1)[0]
        return slope

    
    def _estimate_Voc(self):
        #min_val = np.argmin(abs(np.asarray(self.Currents)))
        if np.all(np.diff(self.Currents) > 0):
            Voc = np.interp(0.0, self.Currents, self.Voltages)
        else:
            Voc = np.interp(0.0, np.flip(self.Currents), np.flip(self.Voltages))
        return Voc#self.Voltages[min_val]
        
    
    def _estimate_Isc0_T0(self):
        min_val = np.argmin(abs(np.asarray(self.Voltages)))
        return self.Currents[min_val]
    

    def fit_function_k(self,k):
        self.temperature_coefficient_internal_series_resistance = k / 1000000.
        self.Tcell = 35.+273.15
        Icell, Vcell, Pcell = self.calcCell()
        Pmax_Tcell1 = np.amax(Pcell)

        self.Tcell = 15.+273.15
        Icell, Vcell, Pcell = self.calcCell()
        Pmax_Tcell2 = np.amax(Pcell)

        self.Tcell = 25.+273.15
        Icell, Vcell, Pcell = self.calcCell()
        Pmax_STC = np.amax(Pcell)

        Gamma_Current1 = (Pmax_Tcell1 - Pmax_STC)/Pmax_STC/10.
        Gamma_Current2 = (Pmax_Tcell2 - Pmax_STC)/Pmax_STC/-10.
        Gamma_Current = np.mean([Gamma_Current1, Gamma_Current2])
        return 1000000.*abs(Gamma_Current - self.gamma_Pmp)



    def calcCell(self):
        """
        Calculate cell I-V curve at TCell and Irradiance.
        Returns (Icell, Vcell, Pcell) : tuple of numpy.ndarray of float
        """
        _delta_T = self.Tcell - self.pvconst.T0  # [K] temperature difference
        Vreverse = self.VRBD * self.pvconst.negpts
        Vff = max(self.Voltages)#self.Voc
        #Vff = self.Voc
        #delta_Voc = self.VocSTC - self.Voc
        delta_Voc = max(self.Voltages) - self.Voc
        # to make sure that the max voltage is always in the 4th quadrant, add
        # a third set of points log spaced with decreasing density, from Voc to
        # Voc @ STC unless Voc *is* Voc @ STC, then use an arbitrary voltage at
        # 80% of Voc as an estimate of Vmp assuming a fill factor of 80% and
        # Isc close to Imp, or if Voc > Voc @ STC, then use Voc as the max
        if delta_Voc == 0:
            Vff = 0.8 * max(self.Voltages)#self.Voc
            delta_Voc = 0.2 * max(self.Voltages)#self.Voc * 1.1
            #Vff = 0.8 * self.Voc
            #delta_Voc = 0.2 * self.Voc * 1.1
        elif delta_Voc < 0:
            Vff = max(self.Voltages)#self.VocSTC# * 1.1
            delta_Voc = -delta_Voc
        Vquad4 = Vff + delta_Voc * self.pvconst.Vmod_q4pts
        Vforward = Vff * self.pvconst.pts
        Vcell = np.concatenate((Vreverse, Vforward, Vquad4), axis=0)
        
        fRBD = 1. - Vcell / self.VRBD
        # use epsilon = 2.2204460492503131e-16 to avoid "divide by zero"
        fRBD[fRBD == 0] = EPS

        fRBD = self.Isc0_T0 * fRBD ** (-self.nRBD)
        Vdiode_norm = Vcell / self.Isc0_T0
        IRBD = (self.aRBD * Vdiode_norm + self.bRBD * Vdiode_norm ** 2) * fRBD
        
        IRBD = np.clip(IRBD, a_min=-3.0*self.Isc0_T0, a_max=None)
        
        Icell = np.interp(Vcell, self.Voltages, self.Currents * self.Ee, left=self.Isc, 
                          right=min(self.Currents * self.Ee))
        Icell -= IRBD
        
        if (_delta_T != 0) or (self.Ee!=1.0):
            Icell *= (1. + self.alpha_Isc * _delta_T)
            Vcell = Vcell * (1. + self.irradiance_correction_factor * np.log(self.Ee)) \
                 + self.beta_Voc(return_absolute=True) * _delta_T \
                 - (self.temperature_coefficient_internal_series_resistance \
                 * (Icell/self.Ee) * _delta_T) \
                 - (self.internal_series_resistance * ((Icell/self.Ee) - Icell))
            
        Pcell = Icell * Vcell
        
        return Icell, Vcell, Pcell


    def plot(self):
        """
        Plot cell I-V curve.
        Returns cellPlot : matplotlib.pyplot figure
        """
        cell_plot = plt.figure()
        plt.subplot(2, 2, 1)
        plt.plot(self.Vcell, self.Icell)
        plt.title('Cell Reverse I-V Characteristics')
        plt.ylabel('Cell Current, I [A]')
        plt.xlim(self.VRBD - 1, 0)
        plt.ylim(0, self.Isc + 10)
        plt.grid()
        plt.subplot(2, 2, 2)
        plt.plot(self.Vcell, self.Icell)
        plt.title('Cell Forward I-V Characteristics')
        plt.ylabel('Cell Current, I [A]')
        plt.xlim(0, self.Voc)
        plt.ylim(0, self.Isc + 1)
        plt.grid()
        plt.subplot(2, 2, 3)
        plt.plot(self.Vcell, self.Pcell)
        plt.title('Cell Reverse P-V Characteristics')
        plt.xlabel('Cell Voltage, V [V]')
        plt.ylabel('Cell Power, P [W]')
        plt.xlim(self.VRBD - 1, 0)
        plt.ylim((self.Isc + 10) * (self.VRBD - 1), -1)
        plt.grid()
        plt.subplot(2, 2, 4)
        plt.plot(self.Vcell, self.Pcell)
        plt.title('Cell Forward P-V Characteristics')
        plt.xlabel('Cell Voltage, V [V]')
        plt.ylabel('Cell Power, P [W]')
        plt.xlim(0, self.Voc)
        plt.ylim(0, (self.Isc + 1) * self.Voc)
        plt.grid()
        plt.tight_layout()
        return cell_plot
