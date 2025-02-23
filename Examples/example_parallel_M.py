
# This file is part of ExoPlex - a self consistent planet builder
# Copyright (C) 2017 - by the ExoPlex team, released under the GNU
# GPL v2 or later.

"""
This example uses parallel processing to quickly calculate the best fit radius and pressure, temperature and density profiles
for a planet of a given mass, its uncertainty and a chosen composition.

The code begins by initializing the composition of the planet and retrieving the grids. In the main text code (at bottom)
one can set the number of samplings and the mass, radius, and uncertainties.

该示例使用并行处理来快速计算最佳的拟合半径和压力，温度和密度曲线
对于给定质量的行星，其不确定性和选择的组成。

该代码首先要初始化行星的组成并检索网格。在主文本代码（在底部）中
可以设置采样数量以及质量，半径和不确定性。

"""
import os
import sys
import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt

# hack to allow scripts to be placed in subdirectories next to exoplex:

if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))

from ExoPlex import functions
from ExoPlex import run_perplex as perp
from ExoPlex import make_grids

Pressure_range_mantle_UM = '1 1400000'
Temperature_range_mantle_UM = '1600 3500'

Pressure_range_mantle_LM = '1250000 40000000'
Temperature_range_mantle_LM = '1700 7000'

water_potential_temp = 300.

comp_keys = ['wt_frac_water','FeMg','SiMg','CaMg','AlMg','wt_frac_FeO_wanted','wt_frac_Si_core',
                          'wt_frac_O_core','wt_frac_S_core', 'combine_phases','use_grids','conserve_oxy']
struct_keys = ['Pressure_range_mantle_UM','Temperature_range_mantle_UM','resolution_UM',
                         'Pressure_range_mantle_LM', 'Temperature_range_mantle_LM', 'resolution_LM',
                         'Mantle_potential_temp','water_potential_temp']
combine_phases = True
use_grids = True

# To have ExoPlex to give you compositional info and status of calculation set Verbose to TRUE.
# Note: setting this to True will slightly slow down the program
verbose = True

# create filename to store values

Output_filename = 'D:\plant-sun\ExoPlex\Examples\example_parallel_M'
# Next user must input the ratios by mole (Earth is Ca/Mg = .07, Si.Mg = 0.90, Al/Mg = 0.09, Fe/Mg = 0.9)
CaMg = 0.07
SiMg = 0.9
AlMg = 0.09
FeMg = 0.9

# How much water do you want in your planet? By mass fraction.
wt_frac_water = 0.0

# Don't forget that if you have water you need to add water layers
number_h2o_layers = 0

# The potential Temperature of Water, if present
water_potential_temp = 300.

# What fraction of the mantle would you like to be made of FeO? This Fe will be pulled from the core.
wt_frac_FeO_wanted = 0.  # by mass
conserve_oxy = False

# Now we can mix various elements into the core or mantle
wt_frac_Si_core = 0.  # by mass <1, note if you conserve oxygen this is calculated for you
wt_frac_O_core = 0.  # by mass
wt_frac_S_core = 0.  # by mass

# What potential temperature (in K) do you want to start your mantle adiabat?
Mantle_potential_temp = 1600.

# Input the resolution of your upper mantle and lower mantle composition, density grids
# These are input as number of T, P points. 50 50 = 2500 grid points, which takes about
# 5 minutes to calculate. Lower mantle resolution does not need to be higher since it's
# mostly ppv.
resolution_UM = '25 75'
resolution_LM = '75 75'

# lastly we need to decide how many layers to put in the planet. This is the resolution of
# the mass-radius sampling.
num_mantle_layers = 300
num_core_layers = 300

Output_radii = []
Output_mass = []

######### Initalize and run ExoPlex


compositional_params = dict(zip(comp_keys, [wt_frac_water, FeMg, SiMg, CaMg, AlMg, wt_frac_FeO_wanted, wt_frac_Si_core, \
                                            wt_frac_O_core, wt_frac_S_core, combine_phases, use_grids, conserve_oxy]))

if use_grids == True:
    filename = functions.find_filename(compositional_params, verbose)
else:
    filename = ''

structure_params = dict(zip(struct_keys, [Pressure_range_mantle_UM, Temperature_range_mantle_UM, resolution_UM,
                                          Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                                          Mantle_potential_temp, water_potential_temp]))

layers = [num_mantle_layers, num_core_layers, number_h2o_layers]
# 获取成分的质量分数、摩尔分数和地核质量分数
Core_wt_per, Mantle_wt_per, Core_mol_per, core_mass_frac = functions.get_percents(compositional_params, verbose)
# 使用 PerPlex 计算地幔的组成和结构，生成 Mantle 的文件
Mantle_filename = perp.run_perplex(*[Mantle_wt_per,compositional_params, structure_params,filename,verbose,combine_phases])
# 基于生成的 Mantle 文件和地幔质量分数，创建地幔的网格
grids_low, names = make_grids.make_mantle_grid(Mantle_filename,Mantle_wt_per, True,use_grids)
# 将 'Fe'（铁）添加到网格的名称列表中
names.append('Fe')

# 如果水层数量大于 0，则计算水的网格和水的相态
if layers[-1] > 0:
    # 使用 make_grids 创建水的网格和水的相态
    water_grid, water_phases = make_grids.make_water_grid()
    # 将所有水的相态添加到名称列表中
    for i in water_phases:
        names.append(i)
else:
    # 如果没有水层，则设置空的水网格
    water_grid = []

# 创建高分辨率的地幔网格
grids_high = make_grids.make_mantle_grid(Mantle_filename, Mantle_wt_per, False, use_grids)[0]

# 创建地核网格
core_grid = make_grids.make_core_grid()

# 将不同分辨率的网格组合成一个网格列表
grids = [grids_low, grids_high, core_grid, water_grid]


def calc_planet(x):
    # 调用函数通过给定的参数计算行星的质量和其他属性
    Planet = functions.find_Planet_mass(x, core_mass_frac, structure_params, compositional_params, grids, Core_wt_per, layers, verbose)

    # 将行星的相态名称列表（如水、岩石等）添加到结果字典中
    Planet['phase_names'] = names

    # 获取行星的相态信息，存储在 `Planet['phases']` 中，`Planet['phase_names']` 更新为新的相态名称列表
    Planet['phases'], Planet['phase_names'] = functions.get_phases(Planet, grids, layers, combine_phases)

    # 使用 check 函数检查行星的计算结果，确保计算符合物理规律（例如，压力应随深度增加）
    functions.check(Planet)

    # 计算并保存行星的半径、质量、核心质量分数和核心半径分数
    rad = Planet['radius'][-1] / 6371e3  # 计算行星的半径，单位为地球半径
    mass = Planet['mass'][-1] / 5.97e24  # 计算行星的质量，单位为地球质量
    CMF = Planet['mass'][num_core_layers - 1] / Planet['mass'][-1]  # 计算核心质量分数（Core Mass Fraction）
    CRF = Planet['radius'][num_core_layers - 1] / Planet['radius'][-1]  # 计算核心半径分数（Core Radius Fraction）

    # 获取行星CMB（地幔-地核边界）的压力和温度
    CMB_P = Planet['pressure'][num_core_layers] / 1e4  # CMB压力，单位为 GPa
    CMB_T = Planet['temperature'][num_core_layers]  # CMB温度，单位为 K

    # 如果水层数量大于 0，则计算水层顶部的压力和温度
    if number_h2o_layers > 0:
        # 计算水层顶部的压力，单位为 GPa
        WMB_P = Planet['pressure'][num_core_layers + num_mantle_layers] / 1e4
        # 计算水层顶部的温度，单位为 K
        WMB_T = Planet['temperature'][num_core_layers + num_mantle_layers]

    # 计算地核的压力，单位为 TPa（压强单位为 Pa，因此除以 1e7）
    P_cen = Planet['pressure'][0] / 1e7
    # 计算地核的温度，单位为 K
    T_cen = Planet['temperature'][0]

    # 如果水层数量大于 0，将水层的相关参数添加到输出字典中
    if number_h2o_layers > 0:
        # 再次获取水层顶部的压力和温度
        WMB_P = Planet['pressure'][num_core_layers + num_mantle_layers] / 1e4
        WMB_T = Planet['temperature'][num_core_layers + num_mantle_layers]
        
        # 定义要返回的键（包括地球的半径、质量、核心质量分数等）
        keys = ['radius','mass','CMF','CRF','CMB_P','CMB_T','P_cen', 'T_cen','WMB_P','WMB_T']
        # 定义对应的值
        vals = [rad, mass, CMF, CRF, CMB_P, CMB_T, P_cen, T_cen, WMB_P, WMB_T]
        
        # 将键和值组成字典并返回
        return(dict(zip(keys, vals)))

    # 如果没有水层，只返回不包含水层相关数据的字典
    else:
        # 定义要返回的键（不包括水层的压力和温度）
        keys = ['radius','mass','CMF','CRF','CMB_P','CMB_T','P_cen', 'T_cen']
        # 定义对应的值
        vals = [rad, mass, CMF, CRF, CMB_P, CMB_T, P_cen, T_cen]
        
        # 将键和值组成字典并返回
        return(dict(zip(keys, vals)))

import json

if __name__ == "__main__":
    num_pts = 3  # 点的数量
    M = 1  # 行星质量（单位：地球质量）
    M_err = .1  # 质量的误差（单位：地球质量）

    # 根据正态分布生成 num_pts 个质量数据点
    Mass_planet = np.random.normal(M, M_err, num_pts)

    # 创建进程池，使用所有可用的 CPU 核心
    pool = mp.Pool(processes=mp.cpu_count())

    # 使用池中的进程并行计算每个质量点对应的行星半径
    Planets = pool.map_async(calc_planet, Mass_planet).get()

    pool.close()
    # 保存 Planets 数据到 JSON 文件
    with open('D:\plant-sun\ExoPlex\Examples\planets_data.json', 'w') as f:
        json.dump([planet.to_dict() for planet in Planets], f, indent=4)
    plt.scatter(Mass_planet, [Planets[i].get('radius') for i in range(num_pts)])
    plt.ylabel('Radius (Earth Radii)', size=20)
    plt.xlabel('Mass (Earth Masses)', size=20)
    plt.show()


