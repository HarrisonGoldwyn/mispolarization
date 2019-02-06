from __future__ import print_function
from __future__ import division

import sys
import os
import numpy as np
import scipy.optimize as opt

## attach module folders to import search path
work_dir = os.getcwd()
date_dir = os.path.split(work_dir)[0]
fluo_drive_folder = os.path.join(date_dir, 'Fluo_drive_bens_params')
field_module_folder = os.path.join(fluo_drive_folder, 'field_functions')
sys.path.append(fluo_drive_folder)
sys.path.append(field_module_folder)
# print('\n'.join(sys.path))

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
import yaml

import dipole_mag_by_fourier as mag

## grab parameter file and appropriate field function file
stream = open('../Fluo_drive_bens_params/four_param.yaml','r')
parameters = yaml.load(stream)
if parameters['general']['fields'] == 'far':
    import far_fields as fi
elif parameters['general']['fields'] == 'full':
    import full_fields as fi

## Generate points on spherical surface
import fibonacci as fib

#### Script vaiables and constants
nm = 1e-7
# hbar = 6.582 * 10**(-16.)


# ## simulated image sensor parameters
# _sensor_size = 1000*nm
# _height = 300*nm
# _resolution = 250  # grid _resolution
lens_points = 10
radius = 300*nm  # focal length of objective lens
max_theta = np.pi/2 # defines physical aperture size

## want printed updates?
updates = True

## animation or individual figures
ani_or_figs = 'figs'  # 'ani' or 'figs'

norm_style = 'per fig'  # 'per fig' or 'across figs'

show_or_save = 'show'
# show_or_save = raw_input('( show / save ) figures? ')

file_extension = 'pdf'  # 'pdf' for mest quality

# coherent = 'y'

## plot parameters
plot_scale = 1/nm

################# functions
def row_mag(vec):
    '''Magnitude of each row of array, returned as column vector.  
    Probably should make sure the array is 2D, but I wont for now.  
    '''
    if vec.ndim == 2:
        mags = np.linalg.norm(vec, axis=(1))
        lenth = mags.size
        arranged = mags.reshape([lenth, 1])
        return arranged
    else: raise ValueError('Can\'t \'row_mag()\' array if ndim != 2')

## function for generating and fitting a 2D gaussian
def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()

# ###################

if __name__ == "__main__":  # call with either 'load' or 'ex' as sys.argv[1]
## loop through seps and amps from mode_amps.py

    ## old block of parameter loading
    # file_title = 'four_magnitudes'
    # date = sys.argv[0][-9:-3]

    ##  Load separations and calculated aplitudes 
    # data_file_name = file_title + '.out'
    # amplitude_data = np.loadtxt(data_file_name)
    seps = mag.seperation_columb(
        parameters['general']['min_sep'],
        parameters['general']['max_sep'],
        parameters['general']['sep_steps']
        )  # column in 2d
    q1x = mag.plas_osc([1,0,0])[0]
    q2x = mag.fluo_osc([1,0,0])[0]
    phase1x = mag.plas_osc([1,0,0])[1]
    phase2x = mag.fluo_osc([1,0,0])[1]

    q1y = mag.plas_osc([0,1,0])[0]
    q2y = mag.fluo_osc([0,1,0])[0]
    phase1y = mag.plas_osc([0,1,0])[1]
    phase2y = mag.fluo_osc([0,1,0])[1]
    

    omega_drive = mag.eV_to_Hz(parameters['general']['drive_energy'])  # driving frequency

    ## dipole positons
    d1_pos_unit = np.array([parameters['plasmon']['pos_unit_vec']])
    d2_pos_unit = np.array([parameters['fluorophore']['pos_unit_vec']])

    ## fix plasmon offset from center
    # d1_from_center = _sensor_size/8
    d1_from_center = 0
    d1_pos = d1_from_center * d1_pos_unit * np.ones(np.shape(seps))
    d2_pos = d2_pos_unit * (seps + d1_pos)

    # ## center fluorophore
    # d2_from_center=0
    # d2_pos = d2_from_center * d2_pos_unit * np.ones(np.shape(seps))
    # d1_pos = d1_pos_unit * (seps + d2_pos)

    # d2_pos = 0 * seps
    # # center plasmon
    # d1_from_center = 0
    # d1_pos = d1_from_center * d1_pos_unit * np.ones(np.shape(seps))
    # d2_pos = d2_pos_unit * (seps + d1_pos)


    ## Define variables used for simulated measurement grid
    # xs = np.linspace(-_sensor_size/2, _sensor_size/2, _resolution)
    # ys = np.linspace(-_sensor_size/2, _sensor_size/2, _resolution)
    # X, Y = np.meshgrid(xs, ys)
    # X_normed = X / (nm)  # assumed definition of X and Y with nm
    # Y_normed = Y / (nm)

    ## scale marking dot size to plot domain for field plots
    # dot_size = (0.025 * _sensor_size/(2*nm)) 
    
    ## points on sphere
    thetas_and_phis = fib.fib_alg_k_filter(lens_points, max_theta)
    xyzs = fib.sphere_to_cart(thetas_and_phis[:,0],thetas_and_phis[:,1],
        radius*np.ones(np.shape(thetas_and_phis[:,0])))

###### initialize arrays to hold field data while looping through seperations

    # meta_data = np.zeros((2,seps.shape[0],_resolution**2))
    field_data = np.zeros((seps.shape[0], xyzs.shape[0], 3))

    # # if what_do == 'g':  # calculate fields and fit gaussians
    #     ## initialize vector for storing maximum intensities
    # max_int = []
    # no_fit = []
    # centroids = np.hstack((seps, np.zeros((seps.shape[0], 3))))
    # mislocalization = np.zeros((seps.shape[0], 1))

###### loop through seperations and calculate fields
    for i in range(seps.shape[0]):
        print('i =',i)
        x1 = d1_pos[i]
        x2 = d2_pos[i]
        # r_d1 = fi.Generate_r(_sensor_size, _resolution, _height, x1)
        # r_d2 = fi.Generate_r(_sensor_size, _resolution, _height, x2)
        r_d1 = xyzs - x1
        r_d2 = xyzs - x2

        x_di_E_field = (
            fi.E_dipole_complex(
                w=omega_drive, 
                dipole_magnitude=q1x[i], 
                dipole_phase=phase1x[i], 
                dipole_ori_unit=[1,0,0], 
                r=r_d1
                )
            +
            fi.E_dipole_complex(
                w=omega_drive, 
                dipole_magnitude=q2x[i], 
                dipole_phase=phase2x[i], 
                dipole_ori_unit=[1,0,0], 
                r=r_d2
                )
            )

        y_di_E_field = (
            fi.E_dipole_complex(
                w=omega_drive, 
                dipole_magnitude=q1y[i], 
                dipole_phase=phase1y[i], 
                dipole_ori_unit=[0,1,0], 
                r=r_d1
                )
            +
            fi.E_dipole_complex(
                w=omega_drive, 
                dipole_magnitude=q2y[i], 
                dipole_phase=phase2y[i], 
                dipole_ori_unit=[0,1,0], 
                r=r_d2
                )
            )

        x_di_B_field = (
            fi.B_dipole_complex(
                w=omega_drive, 
                dipole_magnitude=q1x[i], 
                dipole_phase=phase1x[i], 
                dipole_ori_unit=[1,0,0], 
                r=r_d1
                )
            +
            fi.B_dipole_complex(
                w=omega_drive, 
                dipole_magnitude=q2x[i], 
                dipole_phase=phase2x[i], 
                dipole_ori_unit=[1,0,0], 
                r=r_d2
                )
            )

        y_di_B_field = (
            fi.B_dipole_complex(
                w=omega_drive, 
                dipole_magnitude=q1y[i], 
                dipole_phase=phase1y[i], 
                dipole_ori_unit=[0,1,0], 
                r=r_d1
                )
            +
            fi.B_dipole_complex(
                w=omega_drive, 
                dipole_magnitude=q2y[i], 
                dipole_phase=phase2y[i], 
                dipole_ori_unit=[0,1,0], 
                r=r_d2
                )
            )

        ## Radiation from coupled systems for each orientation
        s_x_di= fi.S_complex(
            x_di_E_field,
            x_di_B_field
            )
        s_y_di= fi.S_complex(
            y_di_E_field,
            y_di_B_field
            )

        ## average poynting vectors from orthoganal dipole orientations
        s = (s_x_di + s_y_di)/2
        # ## or specify a single orientation
        # s = s_x_di
        # s = s_y_di

        ## Define background subtracted Poynting vector
        # S = np.real(s) - np.real(background_s)
        
        # Define magnitudes of the various Poynting vectors
        # mag_bs = np.sum(np.real(background_s)**2,2)**0.5
        if s.ndim == 3: 
            mag_s = np.sum(np.real(s)**2,2)**0.5
            mag_sz = np.abs(s[:,:,2])
        if s.ndim == 2: 
            mag_s = np.sum(np.real(s)**2,1)**0.5
            mag_sz = np.abs(s[:,2])
        # mag_S = np.sum(S**2,2)**0.5
        

        ## Printing max radiation values 
        max_coupled_rad = np.max(np.sum(np.real(s)**2)**0.5)
        # max_plasmon_rad = np.max(np.sum(np.real(background_s)**2)**0.5)
        # max_fluo_rad = np.max(np.sum(np.real(fluo_s)**2)**0.5)
        print(
            # 'seperation = {:06.2f}'.format((x2[0]-x1[0])/nm)+'nm','\n',
            'spacer = {:06.2f}'.format(((x2[0]-x1[0]) 
                - parameters['plasmon']['radius'])/nm)+'nm','\n',
            # 'coupled s max = {:.5e}'.format(max_coupled_rad),'\n',
            # 'background s max = {:.5e}'.format(max_plasmon_rad),'\n',
            # 'fluorophore radiation max = {:.5e}'.format(max_fluo_rad)
            )

########## I just want to plot the magnitude of the electric field
        # decide ONE dipole orientation for now
        e_field_c = x_di_E_field
        e_field = np.real(e_field_c)
        # calculate real magnitude
        if e_field.ndim == 2: 
            mag_e = np.sum((e_field)**2,1)**0.5

        fig = plt.figure(i+1)
        ax = fig.add_subplot(121, projection='3d')
        ax2 = fig.add_subplot(122, projection='3d')

        print(e_field)

        # print(xyzs[:,0].shape, mag_e.shape)
        scatter_plot = ax.scatter(
            xyzs[:,0]/nm, 
            xyzs[:,1]/nm, 
            mag_e, 
            cmap=plt.cm.jet,
            )    

        field_norm_for_plot = mag_e*1e-2
        quiver_plot = ax2.quiver(
            xyzs[:,0]/nm, 
            xyzs[:,1]/nm,
            xyzs[:,2]/nm,
            e_field[:,0]/field_norm_for_plot,
            e_field[:,1]/field_norm_for_plot,
            e_field[:,2]/field_norm_for_plot,
            )

        for axi in [ax, ax2]:
            axi.set_xlabel('x (nm)')
            axi.set_ylabel('y (nm)')
            axi.set_zlabel('z (nm)')

########## Store field data
        field_data[i,:,:] = e_field_c

########## Next step is to perform integrals. 
##########
    # print(field_data)

###### Finish by showing plots    
    # plt.show()
