
# This file is part of ExoPlex - a self consistent planet builder
# Copyright (C) 2017 - by the ExoPlex team, released under the GNU
# GPL v2 or later.

import os
import sys
import numpy as np
# hack to allow scripts to be placed in subdirectories next to exoplex:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))
Pressure_range_mantle_UM = '1 1400000'
Temperature_range_mantle_UM = '1600 3500'

Pressure_range_mantle_LM = '1250000 40000000'
Temperature_range_mantle_LM = '1700 7000'

comp_keys = ['wt_frac_water','FeMg','SiMg','CaMg','AlMg','wt_frac_FeO_wanted','wt_frac_Si_core',
                          'wt_frac_O_core','wt_frac_S_core', 'combine_phases','use_grids','conserve_oxy']
struct_keys = ['Pressure_range_mantle_UM','Temperature_range_mantle_UM','resolution_UM',
                         'Pressure_range_mantle_LM', 'Temperature_range_mantle_LM', 'resolution_LM',
                         'Mantle_potential_temp','water_potential_temp']

water_potential_temp = 300.

combine_phases = True
use_grids = True
verbose = False

import ExoPlex as exo

if __name__ == "__main__":

    #Approximate Masses of TRAPPIST-1 Planets (Agol et al., 2021; ApJ):
    #T1b: 1.374 +/- 0.069
    #T1c: 1.308 +/- 0.056
    #T1d: 0.388 +/- 0.012
    #T1e: 0.692 +/- 0.022
    #T1f: 1.039 +/- 0.031
    #T1g: 1.321 +/- 0.038
    #T1h: 0.326 +/- 0.020

    number_of_runs = 5

    Mass_planet_avg = 1.374
    Mass_planet_sigma = 0.069

    #Need to give the run a name. This will be used as the name of the output files
    Star = 'TRAPPIST-1b'

    #Next user must input the ratios by mole (Earth is Ca/Mg = .07, Si/Mg = 0.90, Al/Mg = 0.09, Fe/Mg = 0.9)
    FeMg = 0.5
    SiMg = 0.9
    AlMg = 0.09
    CaMg = 0.07

    #How much water do you want in your planet? By mass fraction.
    wt_frac_water = 0.

    #Don't forget that if you have water you need to add water layers
    number_h2o_layers = 0

    #Now we can mix various elements into the core or mantle
    wt_frac_Si_core = 0. #by mass <1
    wt_frac_O_core = 0. #by mass
    wt_frac_S_core = 0. #by mass

    #What fraction of the mantle would you like to be made of FeO? This Fe will be pulled from the core.
    wt_frac_FeO_wanted = 0. #by mass
    conserve_oxy = False

    #What potential temperature (in K) do you want to start your mantle adiabat?
    Mantle_potential_temp = 1600.

    #Input the resolution of your upper mantle and lower mantle composition, density grids
    #These are input as number of T, P points. 50 50 = 2500 grid points, which takes about
    #5 minutes to calculate. Lower mantle resolution does not need to be higher since it's
    #mostly ppv.
    resolution_UM = '50 50'
    resolution_LM = '20 20'

    #lastly we need to decide how many layers to put in the planet. This is the resolution of
    #the mass-radius sampling.
    num_mantle_layers = 300
    num_core_layers = 300

    #create filename to store values
    Output_filename = Star + "_baseline"


    #Define output lists
    Output_radii = []
    Output_mass = []
    Output_CMF = [] # Core Mass Fraction
    Output_CRF = [] # Core Radius Fraction
    Output_CMBP = [] # Core Mantle Boundary Pressure
    Output_CMBT = [] #PCore Mantle Boundary Temperature

    #make sure to initalize your compositional parameter value or chi squared lists if
    # you want to include it in your output

    ######### Initalize and run ExoPlex
    for i in range(number_of_runs):
        Mass_planet = np.random.normal(Mass_planet_avg, Mass_planet_sigma, 1)[0]
        print ("Iteration # ", str(i) , "of ", number_of_runs)
        print ("Chosen Mass", str('%.2f'%Mass_planet), "Earth Masses")

        compositional_params = dict(
            zip(comp_keys, [wt_frac_water, FeMg, SiMg, CaMg, AlMg, wt_frac_FeO_wanted, wt_frac_Si_core, \
                            wt_frac_O_core, wt_frac_S_core, combine_phases, use_grids, conserve_oxy]))

        filename = exo.functions.find_filename(compositional_params,verbose)

        structure_params = dict(zip(struct_keys, [Pressure_range_mantle_UM, Temperature_range_mantle_UM, resolution_UM,
                                                  Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                                                  Mantle_potential_temp, water_potential_temp]))


        layers = [num_mantle_layers,num_core_layers,number_h2o_layers]

        #This is where we actually run the planet. First PerPlex grids of mineralogy, density,
        #Cp and alpha are calculated and stored in the Solutions folder. If the file already exists
        #(in name, not necessarily in composition), then PerPlex is not run again.

        Planet = exo.run_planet_mass(Mass_planet, compositional_params, structure_params, layers, filename, verbose)

        #check that everything worked
        exo.functions.check(Planet)

        # Planet is a dictionary containing many parameters of interest:
        # Planet.get('radius') = list of the radial points from calculation (m)
        # Planet.get('mass') = list of the cumulative mass at each radius point from calculation (kg)
        # Planet.get('density') = list of densities from calculation (kg/m^3)
        # Planet.get('temperature') = list of temperature points from calculation (K)
        # Planet.get('gravity') = list of gravity points from calculation (SI)
        # Planet.get('pressure') = list of pressure points from calculation (bar)
        # Planet.get('alpha') = list of values of thermal expansivity points from calculation (1/K)
        # Planet.get('cp') = list of values of specific heat points from calculation (SI)
        # Planet.get('phases') = list of phases and their molar fractions
        print ()
        print ("Mass = ", '%.3f' % (Planet['mass'][-1] / 5.97e24), "Earth masses")
        print ("Radius = ", '%.3f' % (Planet['radius'][-1] / 6371e3), "Earth radii")
        print ("Core Mass Fraction = ", '%.2f' % ( Planet['mass'][num_core_layers] / Planet['mass'][-1]))
        print ("Core Radius Fraction = ", '%.2f' % ( Planet['radius'][num_core_layers] / Planet['radius'][-1]))
        print ("CMB Pressure = ", '%.2f' % (Planet['pressure'][num_core_layers] / 10000), "GPa")
        print ()
        if number_h2o_layers >0:
            print ("WMB pressure", '%.2f' % (Planet['pressure'][num_core_layers+num_mantle_layers] / 10000), "GPa")
        Output_mass.append(Mass_planet)
        Output_radii.append(Planet['radius'][-1]/6371e3)
        Output_CMBP.append(Planet['pressure'][num_core_layers] / 10000)
        Output_CMBT.append(Planet['temperature'][num_core_layers])
        Output_CMF.append(Planet['mass'][num_core_layers] / Planet['mass'][-1])
        Output_CRF.append(Planet['radius'][num_core_layers] / Planet['radius'][-1])

        #If you want to include your compositional parameter And/Or Chi_squared in your output,
        #write the function below as above


    #Stitch together the data for output into a file, if you add more parameters(like chi squared
    # or compositional parameter), make sure to update this line!
    Data = [[Output_mass[i], Output_radii[i],Output_CMF[i],Output_CRF[i],\
             Output_CMBP[i],Output_CMBT[i]] for i in range(number_of_runs)]

    #The Header for the datafile,
    #make sure to update this if you add parameters (like compositional parameter or chi squared)
    #Separate each value by \t
    header = 'Mass\tRadius\tCMF\tCRF\tCMB_P[GPa]\tCMB_T[K]'

    #Save the file
    np.savetxt("Data/"+Output_filename + ".csv", Data, '%.3f', delimiter=',', newline='\n', comments='# ', header=header)
    print ("Basic Data File saved in Data/", Output_filename + '.csv')


    # If you want the full output, uncomment this line. Warning, your files will get cluttered quickly!
    #This file contails all of the radial profiles for density, pressure, temperature and composition
    Output_full_filename = Output_filename + '_full'
    exo.functions.write(Planet,Output_full_filename)

    #If you'd like to skip this, uncomment BOTH this next line with """ and also the very last line at bottom with """
    #"""

    #Now let us plot
    import matplotlib.pyplot as plt

    figure = plt.figure(figsize = (12,10))

    ax1 = plt.subplot2grid((6, 3), (0, 0), colspan=3, rowspan=3)
    ax2 = plt.subplot2grid((6, 3), (3, 0), colspan=3, rowspan=1)
    ax3 = plt.subplot2grid((6, 3), (4, 0), colspan=3, rowspan=1)
    ax4 = plt.subplot2grid((6, 3), (5, 0), colspan=3, rowspan=1)

    ax1.plot(Planet['radius'] / 1.e3, Planet['density'] / 1.e3, 'k', linewidth=2.)
    ax1.set_ylim(0., (max(Planet['density']) / 1.e3) + 1.)
    ax1.set_xlim(0., max(Planet['radius']) / 1.e3)
    ax1.set_ylabel("Density ( $\cdot 10^3$ kg/m$^3$)")

    # Make a subplot showing the calculated pressure profile
    ax2.plot(Planet['radius'] / 1.e3, Planet['pressure'] / 1.e4, 'b', linewidth=2.)
    ax2.set_ylim(0., (max(Planet['pressure']) / 1e4) + 10.)
    ax2.set_xlim(0., max(Planet['radius']) / 1.e3)
    ax2.set_ylabel("Pressure (GPa)")

    # Make a subplot showing the calculated gravity profile
    ax3.plot(Planet['radius'] / 1.e3, Planet['gravity'], 'r', linewidth=2.)
    ax3.set_ylabel("Gravity (m/s$^2)$")
    ax3.set_xlim(0., max(Planet['radius']) / 1.e3)
    ax3.set_ylim(0., max(Planet['gravity']) + 0.5)

    # Make a subplot showing the calculated temperature profile
    ax4.plot(Planet['radius'] / 1.e3, Planet['temperature'], 'g', linewidth=2.)
    ax4.set_ylabel("Temperature ($K$)")
    ax4.set_xlabel("Radius (km)")
    ax4.set_xlim(0., max(Planet['radius']) / 1.e3)
    ax4.set_ylim(0., max(Planet['temperature']) + 100)

    #plt.show()
    #"""