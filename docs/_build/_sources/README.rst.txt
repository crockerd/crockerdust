======================
crockerdust Documentation
======================
This is documentation for crockerdust, found on
`Github  <https://github.com/crockerd/crockerdust>`_.

My first package! Uses effective medium theory to generate extinction models. 

This library has two core objects:

class GEMS(gastronomy.Mineral):

    A class for storing GEMS info for effective medium theory. 
    
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




CmGEMS_arr(gastronomy):

See gastronomy documentation for more details. This is a custom composition class used to interface with gastronomy code. 

This package also uses the following functions:

optconst_perm_all(GEMS_m, GEMS_i, calc_type, **kwargs):
   
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
    (1) CmGEMS object

extinction(fname, energy_arr, MD=1e-4, num_pts=300, size_dist='Powerlaw', scatter_model='Mie',rho=3.0, **kwargs):

    function for generating extinction models from optical constants

    :param fname: file name of optical constant file (must have columns called "energy_eV", "delta", and "beta")
    :type fname: str
    :param energy_arr: array of energies to calculate the model over
    :type energy_arr: np.ndarray
    :param MD: line of sight column density, defaults to 1e-4
    :type MD: float (, optional)
    :param num_pts: the number of points in energy array, defaults to 300
    :type num_pts: int (, optional)
    :param size_dist: size distribution of grain population, defaults to 'Powerlaw'
    :type size_dist: str (, optional)
    :param scatter_model: scattering model used, defaults to 'Mie'
    :type scatter_model: str (, optional)
    ...
    :raises TypeError: raises error if fname is not a str
    ...
    :return: a GrainPop object

 extinction_art(my_cm, energy_arr, MD=1e-4, num_pts=300, size_dist='Powerlaw', scatter_model='Mie', **kwargs):
    
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
    (1) grain pop object from gastronomy

calc_extinction(GEMS_m, GEMS_i, energy_arr, calc_type='Si', MD=1e-4, num_pts=300, size_dist='Powerlaw', scatter_model='Mie', save=False, fname=None, **kwargs):
   
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
    (8) BSi: depletion of Si in the ISM (default = 1)
    (9) BMg: depletion of Mg in the ISM (default = 1)
    (10) BFe: depletion of Fe in the ISM (default = 1)
    (11) ASi: abundance of Si in the ISM (default from gastronomy)
    (12) AMg: abundance of Mg in the ISM (default from gastronomy)
    (13) AFe: abundance of Fe in the ISM (default from gastronomy)
    (14) f_Fei: fraction of Fe in the inclusion, set to 1 for simple case (default = None)

    Returns:
    (1) a SingleGrainPop object
