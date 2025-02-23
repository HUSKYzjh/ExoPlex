
# This file is part of ExoPlex - a self consistent planet builder
# Copyright (C) 2017 - by the ExoPlex team, released under the GNU
# GPL v2 or later.

# USE D:\plant-sun\ExoPlex\ExoPlex

import os
import sys

# hack to allow scripts to be placed in subdirectories next to exoplex:
import numpy as np

exo_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
print('ExoPlex path:', exo_path)
sys.path.insert(1, exo_path)
# print(sys.path)
# if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
#     sys.path.insert(1, os.path.abspath('..'))
#     print('Added parent directory to path')
#     print(sys.path)

# 上地幔（UM）压力范围（单位：帕斯卡，Pa）
Pressure_range_mantle_UM = '1 1400000'  
# 上地幔（UM）温度范围（单位：摄氏度）
Temperature_range_mantle_UM = '1600 3500'

# 下地幔（LM）压力范围（单位：帕斯卡，Pa）
Pressure_range_mantle_LM = '1250000 40000000'  
# 下地幔（LM）温度范围（单位：摄氏度）
Temperature_range_mantle_LM = '1700 7000'

# 水的潜在温度（单位：摄氏度）
water_potential_temp = 300.
# 组成相关的关键字（表示地幔、地核等成分和属性）
comp_keys = ['wt_frac_water','FeMg','SiMg','CaMg','AlMg','wt_frac_FeO_wanted','wt_frac_Si_core',
                          'wt_frac_O_core','wt_frac_S_core', 'combine_phases','use_grids','conserve_oxy']
# 结构相关的关键字（表示压力、温度范围和分辨率等参数）
struct_keys = ['Pressure_range_mantle_UM','Temperature_range_mantle_UM','resolution_UM',
                         'Pressure_range_mantle_LM', 'Temperature_range_mantle_LM', 'resolution_LM',
                         'Mantle_potential_temp','water_potential_temp']
# `combine_phases` 表示是否合并相，这在 `use_grids = False` 时自动设置，不能被改变
combine_phases = True #note this is automatically done when use_grids = False and cannot be changed
# `use_grids` 表示是否使用网格，默认为 `True`，表示启用网格
use_grids = True

import ExoPlex as exo

if __name__ == "__main__":
    #To have ExoPlex to give you compositional info and status of calculation set Verbose to TRUE.
    #Note: setting this to True will slightly slow down the program
    # 为了让 ExoPlex 给出成分信息和计算状态，请将 Verbose 设置为 TRUE。
    # 注意：将该设置为 True 会稍微降低程序的运行速度
    verbose = True

    Mass_planet = 1  # in Earth masses
    # 设定行星质量，单位为地球质量（Earth masses）
    #create filename to store values
    Output_filename = 'D:\plant-sun\ExoPlex\Examples\example_by_mass'
    # Output_filename = 'try'
    #Next user must input the ratios by mole (Earth is Ca/Mg = .07, Si.Mg = 0.90, Al/Mg = 0.09, Fe/Mg = 0.9)
    # 接下来，用户需要输入按摩尔比（mole ratio）设置的成分比例
    # 地球的成分比：Ca/Mg = 0.07，Si/Mg = 0.90，Al/Mg = 0.09，Fe/Mg = 0.9
    CaMg = 0.07
    SiMg = 0.9
    AlMg = 0.09
    FeMg = 0.9

    #How much water do you want in your planet? By mass fraction.
    # 你希望行星中含有多少水？以质量分数表示
    wt_frac_water = 0.0

    #Don't forget that if you have water you need to add water layers
    # 如果你有水，请记得需要添加水层
    number_h2o_layers = 0

    #The potential Temperature of Water, if present
    # 水的潜在温度，如果存在水的话（单位：摄氏度）
    water_potential_temp = 300.

    #What fraction of the mantle would you like to be made of FeO? This Fe will be pulled from the core.
    # 你希望多少质量分数的地幔由 FeO 组成？这些 Fe 会从地核中提取
    wt_frac_FeO_wanted = 0. #by mass
    # 是否保持氧的守恒，如果为 True，某些元素会自动调整以保持氧的守恒
    conserve_oxy = False

    #Now we can mix various elements into the core or mantle
    # 现在可以将各种元素混合到地核或地幔中
    wt_frac_Si_core = 0. #by mass <1, note if you conserve oxygen this is calculated for you
    wt_frac_O_core = 0. #by mass
    wt_frac_S_core = 0. #by mass

    #What potential temperature (in K) do you want to start your mantle adiabat?
    # 你希望地幔的起始潜在温度是多少？（单位：摄氏度）
    Mantle_potential_temp = 1600.

    #Input the resolution of your upper mantle and lower mantle composition, density grids
    #These are input as number of T, P points. 50 50 = 2500 grid points, which takes about
    #5 minutes to calculate. Lower mantle resolution does not need to be higher since it's
    #mostly ppv.
    # 输入上地幔和下地幔的组成、密度网格的分辨率
    # 这些输入为 T、P 点的数量。50 50 = 2500 个网格点，大约需要 5 分钟计算
    # 下地幔的分辨率不需要太高，因为它主要是 ppv（参考状态下的矿物）
    # 上地幔的温度和压力分辨率，表示 25 个温度点和 75 个压力点
    resolution_UM = '25 75'
    resolution_LM = '75 75'

    #lastly we need to decide how many layers to put in the planet. This is the resolution of
    #the mass-radius sampling.
    # 最后，我们需要决定在行星中放置多少个层。这是质量-半径采样的分辨率。
    num_mantle_layers = 300
    num_core_layers = 300

    Output_radii = []
    Output_mass = []


    ######### Initalize and run ExoPlex


    # 将成分参数转换为字典，键值为 comp_keys 列表和相应的值
    compositional_params = dict(zip(comp_keys,[wt_frac_water,FeMg,SiMg,CaMg,AlMg,wt_frac_FeO_wanted,wt_frac_Si_core, \
                          wt_frac_O_core,wt_frac_S_core, combine_phases,use_grids,conserve_oxy]))
    # 根据是否使用网格，决定是否计算文件名
    if use_grids == True:
        filename = exo.functions.find_filename(compositional_params,verbose)
    else:
        filename=''
    # 结构参数，包括地幔和地核的压力、温度范围、分辨率等
    structure_params = dict(zip(struct_keys, [Pressure_range_mantle_UM,Temperature_range_mantle_UM,resolution_UM,
                         Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                                              Mantle_potential_temp,water_potential_temp]))

    # 定义层的数量
    layers = [num_mantle_layers,num_core_layers,number_h2o_layers]

    #This is where we actually run the planet. First PerPlex grids of mineralogy, density,
    #Cp and alpha are calculated and stored in the Solutions folder. If the file already exists
    #(in name, not necessarily in composition), then PerPlex is not run again.
    # 这里开始运行行星模型。首先会计算并存储矿物学、密度、比热容（Cp）和热膨胀系数（alpha）的PerPlex网格
    # 如果文件已经存在（在名字上，而不一定是成分），那么PerPlex将不会再次运行
    Planet = exo.run_planet_mass(Mass_planet, compositional_params, structure_params, layers, filename, verbose)

    #check to see if solution works (P goes up with depth etc.)
    # 检查结果是否有效（如压力是否随着深度增加等）
    exo.functions.check(Planet)

    #Planet is a dictionary containing many parameters of interest:
    #Planet.get('radius') = list of the radial points from calculation (m)
    #Planet.get('mass') = list of the cumulative mass at each radius point from calculation (kg)
    #Planet.get('density') = list of densities from calculation (kg/m^3)
    #Planet.get('temperature') = list of temperature points from calculation (K)
    #Planet.get('gravity') = list of gravity points from calculation (SI)
    #Planet.get('pressure') = list of pressure points from calculation (bar)
    #Planet.get('phases') = list of phases and their molar fractions
    #Planet.get('alpha') = list of values of thermal expansivity points from calculation (1/K)
    #Planet.get('cp') = list of values of specific heat points from calculation (SI)
    # Planet 是一个字典，包含了许多重要的参数信息：
# Planet.get('radius') = 从计算中得到的半径点列表（单位：米）
# Planet.get('mass') = 从计算中得到的每个半径点的累计质量列表（单位：千克）
# Planet.get('density') = 从计算中得到的密度列表（单位：千克/立方米）
# Planet.get('temperature') = 从计算中得到的温度点列表（单位：开尔文）
# Planet.get('gravity') = 从计算中得到的重力点列表（单位：SI）
# Planet.get('pressure') = 从计算中得到的压力点列表（单位：巴）
# Planet.get('phases') = 阶段及其摩尔分数列表
# Planet.get('alpha') = 从计算中得到的热膨胀系数点列表（单位：1/K）
# Planet.get('cp') = 从计算中得到的比热容点列表（单位：SI）
    print()
    print("Mass = ", '%.3f'%(Planet['mass'][-1]/5.97e24), "Earth masses")
    print("Radius = ", '%.5f'%(Planet['radius'][-1]/6371e3), "Earth radii")
    print("Density = ",'%.3f'%(Planet['mass'][-1]/(4/3 * np.pi *pow(Planet['radius'][-1],3))/1000), "g/cc")
    print("Core Mass Fraction = ", '%.2f'%(100.*Planet['mass'][num_core_layers-1]/Planet['mass'][-1]))
    print("Core Radius Fraction = ", '%.2f'%(100.*Planet['radius'][num_core_layers-1]/Planet['radius'][-1]))
    print("CMB Pressure = " ,'%.2f' % (Planet['pressure'][num_core_layers]/1e4), "GPa")
    print("CMB Temperature = " ,'%.2f' % (Planet['temperature'][num_core_layers]), "K")
    if number_h2o_layers >0:
        print("WMB Pressure = ", '%.2f' % (Planet['pressure'][num_core_layers+num_mantle_layers] // 1e4), "GPa")
        print("WMB Temperature = ", '%.2f' % (Planet['temperature'][num_core_layers+num_mantle_layers]), "K")

        print("number of oceans:",'%.2f' % (wt_frac_water*Planet['mass'][-1]/1.4e21))
    print("Central core pressure",'%.2f' % (Planet['pressure'][0]/1e7),"TPa")
    print("Central core Temperature",'%.2f' % (Planet['temperature'][0]),"K")
    print()
    #If you'd like the full output, uncomment out these lines!
    Output_filename = Output_filename + '_Radius_'+ str('%.2f'%(Planet['radius'][-1]/6371e3))
    exo.functions.write(Planet,Output_filename)

    #Now let us plot
    import matplotlib.pyplot as plt

    figure = plt.figure(figsize=(10, 9))

    ax1 = plt.subplot2grid((6, 3), (0, 0), colspan=3, rowspan=3)
    ax2 = plt.subplot2grid((6, 3), (3, 0), colspan=3, rowspan=1)
    ax3 = plt.subplot2grid((6, 3), (4, 0), colspan=3, rowspan=1)
    ax4 = plt.subplot2grid((6, 3), (5, 0), colspan=3, rowspan=1)
    # 绘制地球半径和密度的关系图
    ax1.plot(Planet['radius'] / 1.e3, Planet['density'] / 1.e3, 'k', linewidth=2.)
    ax1.set_ylim(0., (max(Planet['density']) / 1.e3) + 1.)
    ax1.set_xlim(0., max(Planet['radius']) / 1.e3)
    ax1.set_ylabel("Density ( $\cdot 10^3$ kg/m$^3$)")
    ax1.minorticks_on()
    text = '%.3f' % (Planet['radius'][-1] / 6371e3) + ' Earth radii on last iteration'
    ax1.text(0.05, 0.95, text)

    # Make a subplot showing the calculated pressure profile
    # 绘制压力的变化图
    ax2.plot(Planet['radius'] / 1.e3, Planet['pressure'] / 1.e4, 'b', linewidth=2.)
    ax2.set_ylim(0., (max(Planet['pressure']) / 1e4) + 10.)
    ax2.set_xlim(0., max(Planet['radius']) / 1.e3)
    ax2.set_ylabel("Pressure (GPa)")
    ax2.minorticks_on()

    # Make a subplot showing the calculated gravity profile
    # 绘制重力的变化图
    ax3.plot(Planet['radius'] / 1.e3, Planet['gravity'], 'r', linewidth=2.)
    ax3.set_ylabel("Gravity (m/s$^2)$")
    ax3.set_xlim(0., max(Planet['radius']) / 1.e3)
    ax3.set_ylim(0., max(Planet['gravity']) + 0.5)
    ax3.minorticks_on()

    # Make a subplot showing the calculated temperature profile
    # 绘制温度的变化图
    ax4.plot(Planet['radius'] / 1.e3, Planet['temperature'], 'g', linewidth=2.)
    ax4.set_ylabel("Temperature ($K$)")
    ax4.set_xlabel("Radius (km)")
    ax4.set_xlim(0., max(Planet['radius']) / 1.e3)
    ax4.set_ylim(0., max(Planet['temperature']) + 300)
    ax4.minorticks_on()

    plt.show()
