import sys
import os
import numpy as np


import ExoPlex.minphys
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))
Earth_radius = 6371e3
Earth_mass = 5.97e24
ToPa = 100000.
ToBar = 1./ToPa
from ExoPlex import minphys as minphys
from scipy import interpolate

def initialize_by_mass(*args):
    from ExoPlex import burnman as bm

    """
   This module creates the dictionary of lists for each planetary parameter (e.g., density) for a planet of the mass
   input by user.

    Parameters
    ----------
    mass_planet: float
        input radius of planet in Earth radii

    structural_params: list
        Structural parameters of the planet; See example for description

    compositional_params: list
        Structural parameters of the planet; See example for description

    layers: list
        Number of layers for core, mantle and water
    Returns
    -------
    Planet: dictionary
        Dictionary of initial guess of pressure, temperature, expansivity, specific heat and phases for modeled planet
    """
    rock = bm.Composite([bm.minerals.SLB_2011.mg_perovskite(),
                         bm.minerals.SLB_2011.periclase()], \
                        [0.8, 0.2])

    ice = bm.minerals.other.water()
    mass_planet = args[0]
    structural_params = args[1]
    compositional_params = args[2]

    num_mantle_layers, num_core_layers, number_h2o_layers = args[3]
    core_mass_frac = args[4]
    core_grid = args[5][2]
    water_grid = args[5]
    wt_frac_water, FeMg, SiMg, CaMg, AlMg, mol_frac_Fe_mantle, wt_frac_Si_core, \
     wt_frac_O_core, wt_frac_S_core,combine_phases,use_grids,conserve_oxy = compositional_params

    core_rad_frac = structural_params[6]
    Mantle_potential_temp = structural_params[7]
    water_rad_frac = compositional_params[0]
    water_potential_temp = structural_params[9]

    Radius_planet_guess = pow(mass_planet*5.97e24 / 5500. / (4*np.pi/3.),1/3.)/6371e3

    if Radius_planet_guess > 2:
        Radius_planet_guess = 2.

    Radius_planet_guess = 1
    if number_h2o_layers > 0:
        water_thickness_guess = water_rad_frac*Radius_planet_guess*6371e3

        core_thickness_guess = core_rad_frac * (Radius_planet_guess*6371e3-water_thickness_guess)
        mantle_thickness_guess = Radius_planet_guess*6371e3 - water_thickness_guess - core_thickness_guess
    else:
        core_thickness_guess = core_rad_frac * (Radius_planet_guess*6371e3)
        mantle_thickness_guess = Radius_planet_guess*6371e3 - core_thickness_guess

    if (wt_frac_water == 0. and number_h2o_layers > 0) or (number_h2o_layers == 0 and wt_frac_water > 0):
       print ("You have layers of water but no water!")
       number_h2o_layers = 0

    num_layers = num_core_layers+num_mantle_layers + number_h2o_layers # add 50 shells if there is an h2O layer

    # these grids are initialized in this function and are then a passed around through the other routines and edited each iteration
    density_layers = np.zeros(num_layers)
    mass_layers = np.zeros(num_layers)
    Temperature_layers = np.zeros(num_layers)
    radius_layers = np.zeros(num_layers)

    # used for compressiofh2on funciton
    gravity_layers = np.zeros(num_layers)
    Pressure_layers = np.zeros(num_layers)
    alpha = np.zeros(num_layers)
    cp = np.zeros(num_layers)

    water_mass = (wt_frac_water*mass_planet)*Earth_mass
    core_mass = ((mass_planet*(1-wt_frac_water))* core_mass_frac *Earth_mass)
    mantle_mass = (mass_planet*Earth_mass)-water_mass-core_mass

    radius_planet = Radius_planet_guess
    CMB_T_guess = (4180.*(radius_planet-0)-2764.*pow(radius_planet-0,2.)+1219.*pow(radius_planet-0,3.)) + (Mantle_potential_temp-1600)*(0.82+pow(radius_planet-0,1.81))

    if CMB_T_guess < Mantle_potential_temp:
        CMB_T_guess = Mantle_potential_temp + 1000
    if CMB_T_guess > 7000:
        CMB_T_guess = 6800
    CMB_P_guess = 10000*(262.*(radius_planet)-550.*pow(radius_planet,2.) + 432.*pow(radius_planet,3.))

    dP_dr = (CMB_P_guess)/(num_mantle_layers)
    dT_dr = (CMB_T_guess-Mantle_potential_temp)/(num_mantle_layers)


    if wt_frac_water > 0:
        WMB_pres = 1e6 #in bar
    else:
        WMB_pres = 1 #in bar
    for i in range(num_layers):

            if i<number_h2o_layers:
                Pressure_layers[i] = 1 + WMB_pres*(i/number_h2o_layers)
                Temperature_layers[i] = water_potential_temp+500*(i/number_h2o_layers)

            elif i <= number_h2o_layers+num_mantle_layers-1:
                Pressure_layers[i] = WMB_pres + 1.5*dP_dr*(i-number_h2o_layers)
                Temperature_layers[i] = Mantle_potential_temp+1.5*dT_dr*(i-number_h2o_layers)

            else:
                Pressure_layers[i] = Pressure_layers[number_h2o_layers+num_mantle_layers-1] + (i-number_h2o_layers-num_mantle_layers)*5000
                Temperature_layers[i] = Temperature_layers[i-1]+0.5*(i-number_h2o_layers-num_mantle_layers)/num_core_layers




    Pressure_layers = Pressure_layers[::-1]
    Temperature_layers= Temperature_layers[::-1]

    P_core = Pressure_layers[:num_core_layers]
    T_core = Temperature_layers[:num_core_layers]
    rho_core = interpolate.griddata((core_grid['pressure'], core_grid['temperature']),
                                            core_grid['density'], (P_core, T_core), method='linear')
    for i in range(num_layers):
        if i < num_core_layers:
                mass_layers[i] = core_mass/(num_core_layers)
                radius_layers[i] =((float(i) / num_core_layers) * core_thickness_guess)
                density_layers[i] = rho_core[i]



        elif i < (num_core_layers + num_mantle_layers):

            density_layers[i]=(rock.evaluate(['density'], Pressure_layers[i]*ToPa, 300.))

            mass_layers[i] = mantle_mass/num_mantle_layers
            radius_layers[i] =  core_thickness_guess+ ((((i - num_core_layers) / num_mantle_layers) * mantle_thickness_guess))

        else:
            radius_layers[i] = core_thickness_guess + mantle_thickness_guess + \
                               ((float(
                                   i - num_core_layers - num_mantle_layers) / number_h2o_layers) * water_thickness_guess)
            mass_layers[i] = (water_mass / number_h2o_layers)
            density_layers[i] = ice.evaluate(['density'], Pressure_layers[i]*ToPa, 300.)

    mass_update = np.zeros(num_layers)

    for i in range(len(mass_layers)):
        mass_update[i] = (sum(mass_layers[:i + 1]))

    mass_layers = mass_update

    radius_layers[-1] = radius_planet*6371e3
    keys = ['mass','radius', 'density', 'temperature', 'gravity', 'pressure', \
            'alpha', 'cp']

    for i in range(len(density_layers)):
        if np.isnan(density_layers[i])==True:
            print("There is a nan in the initial density")
            print(Pressure_layers[i],Temperature_layers[i])
            sys.exit()

    return dict(zip(keys, [mass_layers, radius_layers,density_layers, Temperature_layers, gravity_layers, Pressure_layers,
                           alpha, cp]))

def compress_mass(*args):
    """
   This module iterates the density, mass within a sphere, adiabatic temperature and gravity integrals for a planet of Mass M
   until convergence is reached. Convergence is defined as the change from the previous run to the current is the
   difference in the density of all layers is <1e-6.

    Parameters
    ----------
    Planet: dictionary
        Dictionary of initial guess of pressure, temperature, expansivity, specific heat and phases for modeled planet

     grids: list of lists
        UM and LM grids containing pressure, temperature, density, expansivity, specific heat and phases

    Core_wt_per: float
        Composition of the Core

    structural_params: list
        Structural parameters of the planet; See example for description

    layers: list
        Number of layers for core, mantle and water

    verbose: bool
    Returns
    -------
    Planet: dictionary
        Dictionary of initial guess of pressure, temperature, expansivity, specific heat and phases for modeled planet
    """
    Planet = args[0]
    grids = args[1]
    Core_wt_per = args[2]
    structural_params= args[3]
    layers= args[4]
    verbose = args[5]
    n_iterations = 1
    max_iterations = 100

    old_r = [10  for i in range(len(Planet['mass']))]
    converge = False

    while n_iterations <= max_iterations and converge == False:
        if verbose == True:
            print ("iteration #",n_iterations)
        if n_iterations>1:
            converge,old_r = minphys.check_convergence(Planet['density'],old_r)

        Planet['density'] = minphys.get_rho(Planet,grids,Core_wt_per,layers)
        Planet['gravity'] = minphys.get_gravity(Planet,layers)
        Planet['pressure'] = minphys.get_pressure(Planet,layers)
        Planet['temperature'] = minphys.get_temperature(Planet, grids, structural_params, layers)
        Planet['radius'] = minphys.get_radius(Planet, layers)

        n_iterations+=1

    return Planet


