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

## colorbar stuff 
from mpl_toolkits import axes_grid1

## Integral stuff
import diffraction_int as diffi
import scipy.integrate as grate

#### Script vaiables and constants
nm = parameters['physical_constants']['nm']
mm = nm*1e6
# hbar = 6.582 * 10**(-16.)
c = parameters['physical_constants']['c']


## simulated image sensor parameters
sensor_size = 2000*nm
# sensor_size = 1
# _height = 300*nm
resolution = 100  # grid _resolution

## Build image sensor
eye = diffi.observation_points(
    x_min= -sensor_size/2, 
    x_max= sensor_size/2,
    y_min= -sensor_size/2, 
    y_max= sensor_size/2, 
    points= resolution
    )

## Experimental parameters
magnification = 1
numerical_aperture = 0.95
## numerical parameters for calculation of scattered field
lens_points = 1000
# radius = 1000*nm  # focal length of objective lens
# radius = 10*nm
max_theta = np.arcsin(numerical_aperture) # defines physical aperture size
obj_f = 2.*mm  # still dont know what this is supposed to be
tube_f = magnification * obj_f

## values required for scattered field calculation of sphere
sphere_points = fib.fib_alg_k_filter(
    num_points=lens_points, 
    max_ang=max_theta
    )

cart_points_on_sph = fib.sphere_to_cart(
    sphere_points[:,0],
    sphere_points[:,1],
    1*np.ones(np.shape(sphere_points[:,0]))
    )

refracted_coords = diffi.refract_sph_coords(sphere_points, obj_f, tube_f)
print('sphere_points= ',sphere_points)
print('refracted_coords= ',refracted_coords)

# ################# functions
# def row_mag(vec):
#     '''Magnitude of each row of array, returned as column vector.  
#     Probably should make sure the array is 2D, but I wont for now.  
#     '''
#     if vec.ndim == 2:
#         mags = np.linalg.norm(vec, axis=(1))
#         lenth = mags.size
#         arranged = mags.reshape([lenth, 1])
#         return arranged
#     else: raise ValueError('Can\'t \'row_mag()\' array if ndim != 2')

# ## function for generating and fitting a 2D gaussian
# def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
#     xo = float(xo)
#     yo = float(yo)    
#     a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
#     b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
#     c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
#     g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
#                             + c*((y-yo)**2)))
#     return g.ravel()

# # ###################

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

    ## Fix plasmon offset from center
    # d1_from_center = _sensor_size/8
    d1_from_center = 0
    d1_pos = d1_from_center * d1_pos_unit * np.ones(np.shape(seps))
    d2_pos = d2_pos_unit * (seps + d1_pos)

    ## Center fluorophore
    # d2_from_center=0
    # d2_pos = d2_from_center * d2_pos_unit * np.ones(np.shape(seps))
    # d1_pos = d1_pos_unit * (seps + d2_pos)

    # d2_pos = 0 * seps
    # # center plasmon
    # d1_from_center = 0
    # d1_pos = d1_from_center * d1_pos_unit * np.ones(np.shape(seps))
    # d2_pos = d2_pos_unit * (seps + d1_pos)


###### initialize arrays to hold field data while looping through seperations

    # meta_data = np.zeros((2,seps.shape[0],_resolution**2))
    field_data = np.zeros((seps.shape[0], cart_points_on_sph.shape[0], 3))

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
        r_d1 = cart_points_on_sph - x1
        r_d2 = cart_points_on_sph - x2

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


########## I just want to plot the magnitude of the electric field
        # decide ONE dipole orientation for now
        e_field_c = x_di_E_field
        e_field = np.real(e_field_c)

        b_field_c = x_di_B_field
        b_field = np.real(b_field_c)
        # calculate real magnitude
        if e_field.ndim == 2: 
            # mag_e = np.sum((e_field)**2,1)**0.5
            # intens = np.sum((e_field)**2,1)
            intens = np.abs(e_field[...,2])
            # print('e_field= ',e_field)
            # print('intensity= ',intens)



        diffracted_E_field = diffi.perform_integral(
            scattered_E=e_field_c, 
            scattered_sph_coords=sphere_points, 
            obser_pts=eye[0], 
            z=0, 
            obj_f=obj_f, 
            tube_f=tube_f, 
            k=omega_drive/c,
            alpha_1_max=max_theta
            )

        diffracted_H_field = diffi.perform_integral(
            scattered_E=b_field_c, 
            scattered_sph_coords=sphere_points, 
            obser_pts=eye[0], 
            z=0, 
            obj_f=obj_f, 
            tube_f=tube_f, 
            k=omega_drive/c,
            alpha_1_max=max_theta
            )

        # print('diffracted_field= ', diffracted_field)
        s = fi.S_complex(
            diffracted_E_field,
            diffracted_H_field
            )
        s_bar = np.real(s)
        image = s_bar[...,2]
        print('image.shape= ',image.shape)
        fig= plt.figure(i+1)
        # print(image)
        # print(eye[1:])
        plot = plt.pcolor(
            eye[1]/(nm * magnification), 
            eye[2]/(nm * magnification), 
            image.reshape(eye[1].shape),
            cmap=plt.cm.jet,)
        ax = plt.axes()
        sep_leg = ax.text(0.5, 0.05,
                (r'spacer = {:.5} nm'.format(
                    (d2_pos[i][0] - d1_pos[i][0]
                        -parameters['plasmon']['radius']
                    )/nm)
                    # +
                    # '\n'+
                    # r'mislocalization = {:.5}'.format(
                    #     mislocalization[i,0])+' nm'
                    ),
                horizontalalignment='center',
                transform=ax.transAxes,
                bbox={'facecolor':'white'}
                )
        ax.set_xlabel('nm in scattering plane')
        # plt.set_ylabel('y (nm)')
        # plt.set_zlabel('z (nm)')
        plt.colorbar()

        # print(e_field)

        # # print(cart_points_on_sph[:,0].shape, mag_e.shape)
        # scatter_plot = ax.scatter(
        #     cart_points_on_sph[:,0]/nm, 
        #     cart_points_on_sph[:,1]/nm, 
        #     mag_e, 
        #     cmap=plt.cm.jet,
        #     )    

        # field_norm_for_plot = mag_e*1e-2
        # quiver_plot = ax2.quiver(
        #     cart_points_on_sph[:,0]/nm, 
        #     cart_points_on_sph[:,1]/nm,
        #     cart_points_on_sph[:,2]/nm,
        #     e_field[:,0]/field_norm_for_plot,
        #     e_field[:,1]/field_norm_for_plot,
        #     e_field[:,2]/field_norm_for_plot,
        #     )

        # for axi in [ax, ax2]:
        #     axi.set_xlabel('x (nm)')
        #     axi.set_ylabel('y (nm)')
        #     axi.set_zlabel('z (nm)')

########## Store field data
        # field_data[i,:,:] = e_field_c

########## Next step is to perform integrals. 
##########
    # print(field_data)

###### Finish by showing plots    
    print('done')
    plt.show()
