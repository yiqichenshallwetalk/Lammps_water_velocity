# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 14:50:16 2019

@author: yiqichen
"""

import numpy as np

def gauss(sigma):
    r=2.
    while r > 1.0 or r == 0.0:
        n=np.random.rand(2)
        v1=2.*n[0]-1.
        v2=2.*n[1]-1.
        r=v1*v1+v2*v2
        
    gauss=v1*np.sqrt(-2.*np.log(r)/r)
    gauss=sigma*gauss
    return gauss

def normal(sigma):
    vx = gauss(sigma)*1.0e-5
    vy = gauss(sigma)*1.0e-5
    vz = gauss(sigma)*1.0e-5
    return vx, vy, vz

def inertia_water():
    mass = [15.9994, 1.008, 1.008]
    angle = 109.47
    rb = 1.0
    x_com = 2*mass[1]*rb*np.cos(angle/2*np.pi/180)/sum(mass)
    x_w = [-x_com, np.cos(angle/2*np.pi/180)-x_com, np.cos(angle/2*np.pi/180)-x_com]
    y_w = [0, np.sin(angle/2*np.pi/180), -np.sin(angle/2*np.pi/180)]
#    coord_O = [-x_com, 0, 0]
#    coord_H1 = [np.cos(angle/2*np.pi/180)-x_com,np.sin(angle/2*np.pi/180),0]
#    coord_H2 = [np.cos(angle/2*np.pi/180)-x_com,-np.sin(angle/2*np.pi/180),0]
    
    jxx = sum([mass[i]*y_w[i]*y_w[i] for i in range(3)])
    jyy = sum([mass[i]*x_w[i]*x_w[i] for i in range(3)])
    jzz = sum([mass[i]*(x_w[i]*x_w[i]+y_w[i]*y_w[i]) for i in range(3)])
    jxy = sum([mass[i]*x_w[i]*y_w[i] for i in range(3)])
    jxz = jyz = 0
    
    J = [[jxx,-jxy,-jxz], [-jxy,jyy,-jyz], [-jxz,-jyz,jzz]]
    return J

def xyz_nonp(x_w, y_w, box_l):
    x_np = []
    y_np = []
    for i in range(len(x_w)):
        dx= x_w[i]-x_w[0]
        x_np.append(x_w[i]-round(dx/box_l)*box_l) 
    for i in range(len(x_w)):
        dy= y_w[i]-y_w[0]
        y_np.append(y_w[i]-round(dy/box_l)*box_l)
    return x_np, y_np
    
def water_v(sigma_w, mass_w, x_w, y_w, z_w, J):
    vO = np.array([gauss(sigma_w[0])*1.0e-5, gauss(sigma_w[0])*1.0e-5, gauss(sigma_w[0])*1.0e-5])
    vH1 = np.array([gauss(sigma_w[1])*1.0e-5, gauss(sigma_w[1])*1.0e-5, gauss(sigma_w[1])*1.0e-5])
    vH2 = np.array([gauss(sigma_w[2])*1.0e-5, gauss(sigma_w[2])*1.0e-5, gauss(sigma_w[2])*1.0e-5])
    P_tot = vO*mass_w[0] + vH1*mass_w[1] + vH2*mass_w[2]
    v_trans = P_tot/sum(mass_w)
    
    vO_nt = vO-v_trans
    vH1_nt = vH1-v_trans
    vH2_nt = vH2-v_trans
    
    # remove period boundary effect
    x_wnp,y_wnp = xyz_nonp(x_w, y_w, 40.0)
    z_wnp = list(z_w)
    
    # compute center of mass
    com = [0.0,0.0,0.0]
    for i in range(len(x_wnp)):
        com[0] += x_wnp[i]*mass_w[i]
        com[1] += y_wnp[i]*mass_w[i]
        com[2] += z_wnp[i]*mass_w[i]
    com = np.array(com)/sum(mass_w)
    
    rO = np.array([x_wnp[0], y_wnp[0], z_wnp[0]]) - com
    rH1 = np.array([x_wnp[1], y_wnp[1], z_wnp[1]]) - com
    rH2 = np.array([x_wnp[2], y_wnp[2], z_wnp[2]]) - com
    
    # Compute angular momentum
    L = mass_w[0]*np.cross(rO,vO_nt)  + mass_w[1]*np.cross(rH1,vH1_nt) + mass_w[2]*np.cross(rH2,vH2_nt)
    
    # Construct matrix A 
    HH = [x_wnp[1]-x_wnp[2], y_wnp[1]-y_wnp[2], z_wnp[1]-z_wnp[2]]
    OH = [(x_wnp[1]+x_wnp[2])/2-x_wnp[0], (y_wnp[1]+y_wnp[2])/2-y_wnp[0], (z_wnp[1]+z_wnp[2])/2-z_wnp[0]]
    CR = np.cross(HH,OH)
    def norm(v):
        sum_v = np.sqrt(sum([i*i for i in v]))
        return v/sum_v
    HH,OH,CR = norm(HH),norm(OH),norm(CR)
    #CR = CR/sum(CR)
    A = [list(OH), list(HH), list(CR)]
    
    # Calculate angular velocity from the formula: w = A^-1 * J^-1 * A * L
    omega = np.linalg.inv(A).dot(np.linalg.inv(J)).dot(A).dot(L)
    
    # Compute final velocities
    vO_rigid = np.cross(omega, rO) + v_trans
    vH1_rigid = np.cross(omega, rH1) + v_trans
    vH2_rigid = np.cross(omega, rH2) + v_trans
    
    return vO_rigid, vH1_rigid, vH2_rigid
    
    
    
    