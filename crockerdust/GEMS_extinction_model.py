#import packages
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from gastronomy.minerals import Mineral
import gastronomy
from astropy.table import Table, vstack
import astropy.units as u
import astropy.constants as c
import newdust

from .cmGEMS import CmGEMS, CmGEMS_arr


#simple effective medium theory function (f_Fei=1)
def optconst_perm(delta_m, beta_m, delta_d, beta_d, Si_Fe, dep_Si=0.1, dep_Fe=0.1, density_Si=3.6, density_Fe=5.3, 
                  mu_Si=100.3887, mu_Fe=159.6882, X_Si=1, X_Fe=2):
    '''
    Calculates GEMS optical constants

    Requires 5 arguments:
    (1) an array of deltas for the silicate
    (2) an array of betas for the silicate
    (3) an array of deltas for the nanoparticle
    (4) an array of betas for the nanoparticle
    (5) the ratio of Si/Fe in the ISM

    Optional arguments:
    (1) dep_Si: depletion of Fe in the ISM (default = 0.1)
    (2) dep_Fe: depletion of Si in the ISM (default = 0.1)
    (3) density_Si: density of the silicate (default = 3.6 g/cm^3)
    (4) density_Fe: density of the inclusion (default = 5.3 g/cm^3)
    (5) mu_Si: mean molecular weight of the silicate (default = 100.3887)
    (6) mu_Fe: mean molecular weight of the inclusion (default = 159.6882)
    (7) X_Si: number of Si atoms per silicate particle (default = 1)
    (8) X_Fe: number of Fe atoms per inclusion particle (default = 2)
    
    Returns:
    (1) an array of permitivities
    '''
    #d is the Fe, m is the silicate
    # calculate epsilons
    ed = np.array([1-2*d-d**2-b**2+2*b*(1-d)*1j for b,d in zip(beta_d, delta_d)])
    em = np.array([1-2*d-d**2-b**2+2*b*(1-d)*1j for b,d in zip(beta_m, delta_m)])
    
    #find VSi/VFe
    VSi_Fe = (mu_Si/(density_Si*X_Si))*((density_Fe*X_Fe)/mu_Fe)*Si_Fe*(dep_Si/dep_Fe)

    #find volume fractions
    cd = 1/(1+VSi_Fe)
    cm = 1-cd

    #calculate Hb
    Hb = (3*cd-1)*ed+(3*cm-1)*em

    #now find effective permitivity
    eeff = (Hb+np.sqrt(Hb**2+8*ed*em))/4

    #convert back to delta and beta
    n = np.sqrt(eeff)
    delta_comb = 1-np.real(n)
    beta_comb = np.imag(n)
    return eeff, delta_comb, beta_comb

#complex effective medium theory 

#create GEMS class
class GEMS(Mineral):
    '''
    A class for storing GEMS info for effective medium theory
    
    Requires 5 inputs:
    (1) a dictionary describing the number of atoms of each element type in the molecule
    (2) the density of the material
    (3) an array of deltas
    (4) an array of betas
    (5) an array of energies for the optical constants

    Optional inputs:
    (1) name: the name of the particle (default = None)

    Attributes:
    (1) name: the name of the particle
    (2) composition: a dictionary describing the number of atoms of each element type given in the input
    (3) weight_amu: the weight of the molecule in atomic mass units (amu)
    (4) weight: the mass of the molecule
    (5) elements: returns the elements that make up the keys for self.composition
    (6) unit_mass: eturns the mass of the mineral unit cell
    (7) rho: the density of the material in g/cm^3
    
    '''
    def __init__(self, composition, rho, delta, beta, energy_arr, name=None):

        #get mineral properties
        Mineral.__init__(self, composition, name=name)

        #add density
        self.rho = rho

        #add optical constants
        self.beta = beta
        self.delta = delta
        self.energy = energy_arr

def calc_vol_frac(GEMS_m, GEMS_i, calc_type, BSi=1, BMg=1, BFe=1, ASi=gastronomy.abundances.get_dust_abund('Si'), 
                  AMg=gastronomy.abundances.get_dust_abund('Mg'), AFe=gastronomy.abundances.get_dust_abund('Fe'), 
                  f_Fei=None):
    '''
    Calculates the volume fractions for effective medium theory calculations.

    Requires 3  arguments:
    (1) a GEMS object for the silicate
    (2) a GEMS object for the inclusion
    (3) toggle between Mg and Si constraints

    Optional arguments:
    (1) BSi: depletion of Si in the ISM (default = 1)
    (2) BMg: depletion of Mg in the ISM (default = 1)
    (3) BFe: depletion of Fe in the ISM (default = 1)
    (4) ASi: abundance of Si in the ISM (default from gastronomy)
    (5) AMg: abundance of Mg in the ISM (default from gastronomy)
    (6) AFe: abundance of Fe in the ISM (default from gastronomy)
    (7) f_Fei: fraction of Fe in the inclusion, set to 1 for simple case (default = None)

    Returns:
    (1) volume fraction of the inclusion
    (2) volume fraction of the silicate
    '''
    #check arguments
    assert calc_type == 'Mg' or calc_type == 'Si',"Please choose either 'Si' or 'Mg' for this argument"

    #for Mg constraint 
    if calc_type == 'Mg':
        #print('Mg abundance: {0:0.2e}'.format(AMg))
        if f_Fei == None:
            f_Fei = 1 - ((AMg*BMg*GEMS_m.composition['Fe']) / (AFe*BFe*GEMS_m.composition['Mg']))
        else:
            f_Fei = f_Fei

        #Calculate volume fractions
        C = (AMg*BMg*GEMS_i.composition['Fe']*GEMS_m.weight_amu) / (AFe*BFe*GEMS_m.composition['Mg']*GEMS_i.weight_amu)

        Vi_Vg = (1 + (C/f_Fei)*(GEMS_i.rho/GEMS_m.rho))**-1
        Vm_Vg = 1 -  Vi_Vg
        
    elif calc_type == 'Si':
        #print('Si abundance: {0:0.2e}'.format(ASi))
        if f_Fei == None:
            f_Fei = 1 - ((ASi*BSi*GEMS_m.composition['Fe']) / (AFe*BFe*GEMS_m.composition['Si']))
        else:
            f_Fei = f_Fei

        #Calculate volume fractions
        C = (ASi*BSi*GEMS_i.composition['Fe']*GEMS_m.weight_amu) / (AFe*BFe*GEMS_m.composition['Si']*GEMS_i.weight_amu)

        Vi_Vg = (1 + (C/f_Fei)*(GEMS_i.rho/GEMS_m.rho))**-1
        Vm_Vg = 1 -  Vi_Vg
        
    return Vi_Vg, Vm_Vg, f_Fei


def optconst_perm_all(GEMS_m, GEMS_i, calc_type, **kwargs):
    '''
    Calculates effective medium theory optical constants for both simple and complex cases.

    Requires 3 arguments
    (1) a GEMS object for the silicate
    (2) a GEMS object for the inclusion
    (3) the type of constraint used for the calculation (either Si or Mg)

    Optional arguments:
    (1) BSi: depletion of Si in the ISM (default = 1)
    (2) BMg: depletion of Mg in the ISM (default = 1)
    (3) BFe: depletion of Fe in the ISM (default = 1)
    (4) ASi: abundance of Si in the ISM (default from gastronomy)
    (5) AMg: abundance of Mg in the ISM (default from gastronomy)
    (6) AFe: abundance of Fe in the ISM (default from gastronomy)
    (7) f_Fei: fraction of Fe in the inclusion, set to 1 for simple case (default = None)
    
    Returns:
    (1) an array of permitivities
    '''

    #checks
    assert len(GEMS_m.energy) == len(GEMS_i.energy),"energy arrays for GEMS objects must be the same length"
    
    #d is the Fe, m is the silicate
    # calculate epsilons
    ed = np.array([1-2*d-d**2-b**2+2*b*(1-d)*1j for b,d in zip(GEMS_i.beta, GEMS_i.delta)])
    em = np.array([1-2*d-d**2-b**2+2*b*(1-d)*1j for b,d in zip(GEMS_m.beta, GEMS_m.delta)])
    
    #find volume fractions 
    cd, cm, f_Fei = calc_vol_frac(GEMS_m, GEMS_i, calc_type, **kwargs)

    #calculate Hb
    Hb = (3*cd-1)*ed+(3*cm-1)*em

    #now find effective permitivity
    eeff = (Hb+np.sqrt(Hb**2+8*ed*em))/4

    #convert back to delta and beta
    n = np.sqrt(eeff)
    delta_comb = 1-np.real(n)
    beta_comb = np.imag(n)
    rho = cd*GEMS_i.rho + cm*GEMS_m.rho
    energy = GEMS_m.energy 
    
    return CmGEMS_arr(energy, delta_comb, beta_comb, rho)
#generate extinction cross sections
#function for extinction monels

def extinction(fname, energy_arr, MD=1e-4, num_pts=300, size_dist='Powerlaw', scatter_model='Mie',rho=3.0, **kwargs):
    '''
    function for generating extinction models from optical constants

    requires 2 arguments:
    (1) file name of optical constant file (must have columns called "energy_eV", "delta", and "beta")
    (2) array of energies to calculate the model over

    Optional arguments
    (1) MD: line of sight column density (default = 1e-4)
    (2) num_pts: the number of points in energy array (default = 300)
    (3) size_dist: size distribution of grain population (default = 'Powerlaw')
    (4) scatter_model: scattering model used (default = 'Mie')

    Returns
    (1) array of energies
    (2) array of total extinctions
    (3) array of scattering extinctions
    (4) array of absorption extinctions
    '''
    #check function dependencies
    assert type(fname) == str,"file name must be a string"

    #set up model
    my_cm = CmGEMS(fname,rho)
    gpop = newdust.SingleGrainPop(size_dist, my_cm, scatter_model, md=MD)

    #calculate model
    gpop.calculate_ext(energy_arr)
    
    #return relevant quantities
    return gpop

def extinction_arr(my_cm, energy_arr, MD=1e-4, num_pts=300, size_dist='Powerlaw', scatter_model='Mie', **kwargs):
    '''
    function for generating extinction models from optical constants

    requires 2 arguments:
    (1) a CmGEMS object
    (2) array of energies to calculate the model over (must have units!)

    Optional arguments
    (1) MD: line of sight column density (default = 1e-4)
    (2) num_pts: the number of points in energy array (default = 300)
    (3) size_dist: size distribution of grain population (default = 'Powerlaw')
    (4) scatter_model: scattering model used (default = 'Mie')

    Returns
    (1) array of energies
    (2) array of total extinctions
    (3) array of scattering extinctions
    (4) array of absorption extinctions
    '''

    #set up model
    gpop = newdust.SingleGrainPop(size_dist, my_cm, scatter_model, md=MD)

    #calculate model
    gpop.calculate_ext(energy_arr)
    
    #return relevant quantities
    #return gpop.lam, gpop.tau_ext, gpop.tau_sca, gpop.tau_abs 
    return gpop

def calc_extinction(GEMS_m, GEMS_i, energy_arr, calc_type='Si', MD=1e-4, num_pts=300, size_dist='Powerlaw', scatter_model='Mie', save=False, fname=None, **kwargs):
    '''
    Calculates extinction profile from desired GEMS composition

    Requires 3 arguments:
    (1) a GEMS object for the silicate
    (2) a GEMS object for the inclusion
    (3) an array of energies for evaluating the extinction.

    Optional arguments:
    (1) calc_type: choose the constraint used for the volume fraction calculation (default = 'Si')
    (2) MD: line of sight column density (default = 1e-4)
    (3) num_pts: the number of points in energy array (default = 300)
    (4) size_dist: size distribution of grain population (default = 'Powerlaw')
    (5) scatter_model: scattering model used (default = 'Mie')
    (6) save: if true, save file with fname (default = False)
    (7) fname: file name for the GEMS cross section output (default = None)
    (7) BSi: depletion of Si in the ISM (default = 1)
    (8) BMg: depletion of Mg in the ISM (default = 1)
    (9) BFe: depletion of Fe in the ISM (default = 1)
    (10) ASi: abundance of Si in the ISM (default from gastronomy)
    (11) AMg: abundance of Mg in the ISM (default from gastronomy)
    (12) AFe: abundance of Fe in the ISM (default from gastronomy)
    (13) f_Fei: fraction of Fe in the inclusion, set to 1 for simple case (default = None)

    Returns:
    (1) an array of energies where the extinction was calculated
    (2) total extinction array
    (3) scattering array
    (4) absorption array
    '''
    #calculate optical constants
    my_cm = optconst_perm_all(GEMS_m, GEMS_i, calc_type, **kwargs)
    print(my_cm.rho)
    #save file
    if save == True:
        assert fname is not None,"please choose a filename"
        t = Table()
        t['energy_eV'] = my_cm.wavel
        t['delta'] = my_cm.delta
        t['beta'] = my_cm.imvals #beta=k
        t.write(fname, format='ascii', overwrite=True)
    
    #now do extinction
    result = extinction_arr(my_cm, energy_arr, MD=MD, num_pts=num_pts, size_dist=size_dist, 
                            scatter_model=scatter_model)

    return result
    
