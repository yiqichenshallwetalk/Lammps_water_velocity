import pandas as pd
import numpy as np
from velocity import *

def generate_velocity(T):
    
    kb = 1.38064852e-23
    beta = 1.0/T/kb
    amu = 1.660539040e-27    
    
    with open('water.data','r') as f:
        for num, line in enumerate(f,1):
            if 'Atoms' in line:
                startl = num + 2
            if 'Velocities' in line:
                endl = num - 2 
    
    natoms = endl - startl + 1
    
    df_water = pd.read_csv('water.data', sep = ' ' ,engine= 'python', names = ["count","molid", "atomty", "charges", "x", "y", "z", "nx", "ny", "nz"],  skiprows = startl-1, nrows = natoms)
    df_water.sort_values(by=['count'], inplace =True)
    natoms = len(df_water)
    
    mass = np.array([15.9994, 1.008, 1.008])*amu
    sigma = 1.0/np.sqrt(mass*beta)
    df_water['sigma'] = list(sigma)*(len(df_water)//3)
    J = inertia_water()

    with open('water.data','r') as f: 
        with open('water_new.data','w') as fp:
            for line in f:
                if 'Velocities' in line:
                    break
                else:
                    fp.write(line)
                    
    with open('water_new.data','a') as fp:
        print("Velocities", file=fp)
        print(file=fp)
        line = "{0:6d} {1:8.10f} {2:8.10f} {3:8.10f}"
        for i in range(len(df_water)//3):
            sigma_w = list(df_water['sigma'][3*i:3*(i+1)])
            x_w, y_w, z_w = list(df_water['x'][3*i:3*(i+1)]), list(df_water['y'][3*i:3*(i+1)]), list(df_water['z'][3*i:3*(i+1)])
            vO, vH1, vH2 = water_v(sigma_w,[15.9994, 1.008, 1.008], x_w, y_w, z_w, J)
            print(line.format(3*i+1, vO[0], vO[1], vO[2]), file=fp)  
            print(line.format(3*i+2, vH1[0], vH1[1], vH1[2]), file=fp)
            print(line.format(3*i+3, vH2[0], vH2[1], vH2[2]), file=fp)
        print(file=fp)
        
    with open('water.data','r') as f:
        with open('water_new.data','a') as fp:
            for line in f:
                if 'Bonds' in line:
                    fp.write( "Bonds \n")
                    for line in f:
                        fp.write(line)
                        
T = 500.0
generate_velocity(T)
                        