''' must be called from a directory 1 below date'''
from __future__ import print_function
from __future__ import division

import sys
import numpy as np
# import scipy.optimize as opt
import matplotlib.pyplot as plt
# import matplotlib.animation as animation
# import fields_020117 as fi
import matplotlib.patches as patches
# import os
import yaml
import os

work_dir = os.getcwd()
# print(work_dir)
date_dir = os.path.split(work_dir)[0]
## should eventually put procedue here to find the date folder
## regardeless of position in path.

try:
    full_path_to_yaml = os.path.join(work_dir,
        'four_param.yaml')

    stream = open(full_path_to_yaml,'r')
    parameters = yaml.load(stream)
except IOError:
    try:
        full_path_to_yaml = os.path.join(work_dir,
            'small_param.yaml')

        stream = open(full_path_to_yaml,'r')
        parameters = yaml.load(stream)
    except IOError:
        full_path_to_yaml = os.path.join(work_dir,
            '../curly_param.yaml')

        stream = open(full_path_to_yaml,'r')
        parameters = yaml.load(stream)




full_path_to_constant_yaml = os.path.join(date_dir,'physical_constants.yaml')
opened_constant_file = open(full_path_to_constant_yaml,'r')
constants = yaml.load(opened_constant_file)
e = constants['physical_constants']['e']
c = constants['physical_constants']['c']  # charge of electron in statcoloumbs
hbar =constants['physical_constants']['hbar']
nm = constants['physical_constants']['nm']


## System background
n_b = parameters['general']['background_ref_index']
eps_b = n_b**2.


def Generate_r(sensor_size, resolution, height, origin):
    '''Takes 3 inputs and generates a vector field with values r
    representing displacement from the origin on obsevation plane. 
    #1 is the size of the square sensor in meters 
    #2 grid resolution
    #3 height above dipoles
    Shape of r is (#1, #1, 3), which is handled properly by np.cross().
    Last axis is cartesien vector.
    '''
    x = np.linspace(
        -sensor_size/2, 
        sensor_size/2, 
        resolution
        )
    y = np.linspace(
        -sensor_size/2, 
        sensor_size/2, 
        resolution
        )
    x = x-origin[0]
    y = y-origin[1]

    y_33 = y[:,np.newaxis]*np.ones(np.shape(y))
    y_331 = y_33[:,:,np.newaxis]
    y_333 = y_331 * np.array([0,1,0])

    x_33 = x[:,np.newaxis]*np.ones(np.shape(y))
    x_331 = x_33.T[:,:,np.newaxis]
    x_333 = x_331 * np.array([1,0,0])

    z_333 = height * np.ones(np.shape(x_331)) * np.array([0,0,1])
    r_333 = x_333 + y_333 + z_333
    return r_333

def Field_magnitude(vec_field):
    '''Takes field of shape=(x, y, 3). 
    Returns magnitude field of shape=(x, y, 1)
    '''
    if vec_field.ndim == 3:
        mag_field = np.sqrt(np.sum(vec_field**2.,2))
        mag_field = mag_field[:,:,np.newaxis]
    elif vec_field.ndim == 2:
        mag_field = np.sqrt(np.sum(vec_field**2.,1))
        mag_field = mag_field[:,np.newaxis]
    else: 
        mag_field = np.sqrt(np.sum(vec_field**2.,-1))
        mag_field = mag_field[...,np.newaxis]
    return mag_field

def Unit_vector(vec_field):
    unit = vec_field * (1/Field_magnitude(vec_field))
    return unit


## Electromagnetic Vector Fields
def B_dipole_complex(w, dipole_magnitude, dipole_phase, dipole_ori_unit, r):
    ''' Computes the complex harmonic field E of the dipole with amplitude 
    (#2) and phase (#3) oscillating at frequency #1.
    The physical field is the real part of the returned object. 
    '''
    n = Unit_vector(r)
    r_mag = Field_magnitude(r)
    k = w * n_b / c
    phase = dipole_phase + k*r_mag
    prefactor = (k**2./r_mag)*np.exp(1j*phase)
    vector = np.cross(
        n,
        dipole_ori_unit
        )
    b_field = prefactor*vector*dipole_magnitude
    return b_field

def E_dipole_complex(w, dipole_magnitude, dipole_phase, dipole_ori_unit, r):
    ''' Computes the complex harmonic field E of the dipole with amplitude 
    (#2) and phase (#3) oscillating at frequency #1.
    The physical field is the real part of the returned object. 
    '''
    n = Unit_vector(r)
    r_mag = Field_magnitude(r)
    
    B = B_dipole_complex(
        w, 
        dipole_magnitude, 
        dipole_phase, 
        dipole_ori_unit, 
        r
        )

    E = np.cross(B, n)

    # k = w * n_b / c
    # phase = dipole_phase + k*r_mag
    # prefactor = (k**2/r_mag)*np.exp(1j*phase)
    # vector = np.cross(
    #     np.cross(
    #         n, 
    #         dipole_ori_unit
    #         ),
    #     n
    #     )
    # e_field = prefactor*vector*dipole_magnitude
    return E


## quadrapole far-fields
def B_quad_complex(
    w, quad_mom_magnitude, quad_mom_phase, dipole_ori_tensor_index_tuple, r):
    ''''''
    n = Unit_vector(r)
    r_mag = Field_magnitude(r)
    k = w * n_b / c
    phase = quad_mom_phase + k*r_mag
    prefactor = (-1j*k**3./6.)*(np.exp(1j*phase)/r_mag)
    if dipole_ori_tensor_index_tuple == (1,1):
        ## (n cross xhat )(n dort xhat)qxx
        vector = (
            np.cross(
                n, 
                np.array([1,0,0])
                )
            * 
            np.dot(
                n, 
                np.array([1,0,0])
                )[:,None]
            )

    quadrapole_moment_component = quad_mom_magnitude

    B = prefactor * vector * quadrapole_moment_component
    return B

def E_quad_complex(
    w, quad_mom_magnitude, quad_mom_phase, dipole_ori_tensor_index_tuple, r):
    ''''''
    n = Unit_vector(r)
    B = B_quad_complex(
            w, 
            quad_mom_magnitude, 
            quad_mom_phase, 
            dipole_ori_tensor_index_tuple, 
            r
            )
    E = np.cross(B, n)

    return E

def S_complex(e_field, b_field):
    prefactor = c/(8*np.pi)
    vector = np.cross(
        e_field,
        np.conj(b_field)
        )
    s_field = prefactor*vector
    return s_field










