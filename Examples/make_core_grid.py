import os
import sys
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))
import numpy as np
import ExoPlex.burnman as bm
import math
if __name__ == "__main__":

    filename_out = '2liquid_iron_grid'
    header = 'P_bar,T_K,rho_kgm3,alpha_1_K,Cp_J_K_mol'
    rock = bm.minerals.other.liquid_iron()
    output = []
    counter = 0

    Pressure = [(10+i*75)*1e9 for i in range(100)]
    Temperature = [1700+100*i for i in range(88)]

    counter = 0


    counter_full = len(Pressure)
    for i in Pressure:
        print("on P",counter," of ",counter_full)
        for j in Temperature:
            #print(i/1e9,j)
            # 调用 ExoPlex 中的 liquid_iron 类评估密度、热膨胀系数和比热容
            rho, alpha, cp = rock.evaluate(['density','thermal_expansivity','molar_heat_capacity_p'],i,j)
            # rock.evaluate 方法是 ExoPlex.burnman 模块中的一个重要方法，用于计算给定压力和温度下液态铁的物理性质，如密度（density）、热膨胀系数（thermal_expansivity）和比热容（molar_heat_capacity_p）
            output.append([i*1.e-5,j,rho,math.log10(alpha),cp*(1000./55.845)])

        # 每处理 5 个压力点，就将数据保存到文件中
        if (counter+1)%5 ==0:
            np.savetxt(filename_out + '.dat', output, delimiter=",", newline='\n', fmt = '%.5e',
                   header=header, footer='', comments='# ')

        counter+=1
