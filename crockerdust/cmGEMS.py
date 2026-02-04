import numpy as np
import astropy.units as u
from astropy.io import ascii
from newdust.graindist.composition import Composition

COMPS = ['pyrite', 'pyrrhotite', 'troilite'] # names of compounds available
RHOS  = dict(zip(COMPS, [5.0, 4.6, 4.6])) # density of each compound, g cm^-3
# https://webmineral.com/data/Pyrite.shtml
# https://webmineral.com/data/Pyrrhotite.shtml
# https://www.webmineral.com/data/Troilite.shtml

class CmGEMS(Composition):
    def __init__(self, compstring, rho=3.0):
        #assert(compstring in COMPS)
        Composition.__init__(self)
        self.cmtype = compstring
        self.rho = rho

        data = ascii.read(compstring)
        eV = data['energy_eV']
        one_minus_n = data['delta'] # 1 - Real part
        k = data['beta'] # imaginary part

        self.wavel = eV * u.eV
        self.revals = 1 - one_minus_n
        self.imvals = k



class CmGEMS_arr(Composition):
    def __init__(self, energy, delta, beta, rho=3.0):
        #assert(compstring in COMPS)
        Composition.__init__(self)
        self.cmtype = energy
        self.rho = rho

        eV = energy
        one_minus_n = delta # 1 - Real part
        k = beta # imaginary part

        self.wavel = eV * u.eV
        self.delta = one_minus_n
        self.revals = 1 - one_minus_n
        self.imvals = k