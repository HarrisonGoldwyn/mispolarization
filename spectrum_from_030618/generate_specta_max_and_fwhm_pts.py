from __future__ import print_function
from __future__ import division

import sys
import os
import numpy as np
import scipy.optimize as opt
import yaml

## import diffrantion integral solver from Optics folder
work_dir = os.getcwd()
optics_folder = os.path.join(work_dir, 'Optics')
sys.path.append(optics_folder)
import diffraction_int as diffi
import fibonacci as fib

## Import field functions
field_module_folder = os.path.join(work_dir, 'field_functions')             
sys.path.append(field_module_folder)
import far_fields as fi

stream = open('param.yaml','r')
parameters = yaml.load(stream)


# elif parameters['general']['fields'] == 'full':
#     import full_fields as fi

## plotting stuff
import matplotlib.pyplot as plt
# import matplotlib.animation as animation
import matplotlib.patches as patches
# from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages



## Define function for saving all figures to one pdf
def multipage(filename, figs=None, dpi=200):
    pp = PdfPages(filename)
    if figs is None:
        figs = [plt.figure(n) for n in plt.get_fignums()]
    for fig in figs:
        fig.savefig(pp, format='pdf')
    pp.close()

## colorbar stuff 
from mpl_toolkits import axes_grid1

## Fourier solution to dipole dynamics
# import dipole_mag_by_fourier as mag
import osc_prob_slns as osc

## Load physical constants from yaml file
full_path_to_constant_yaml = os.path.join(work_dir,'physical_constants.yaml')
opened_constant_file = open(full_path_to_constant_yaml,'r')
constants = yaml.load(opened_constant_file)
e = constants['physical_constants']['e']
c = constants['physical_constants']['c']  # charge of electron in statcoloumbs
hbar =constants['physical_constants']['hbar']
nm = constants['physical_constants']['nm']
n_a = constants['physical_constants']['nA']   # Avogadro's number

mm = 1e6*nm
n_b = parameters['general']['background_ref_index']  ## not used right now
eps_b = n_b**2.

fic_driving_field_strength = parameters['general']['drive_amp'] 
inc_field_intens = (c/(8*np.pi))*fic_driving_field_strength**2.

#######################################################################
## simulated image 
sensor_size = 1000*nm

height = 2*mm  # also defines objective lens focal length 

resolution = 100  # image grid resolution

#######################################################################
## Optics stuff.  

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
numerical_aperture = 0.999
## numerical parameters for calculation of scattered field
lens_points = 300
# radius = 1000*nm  # focal length of objective lens
# radius = 10*nm
max_theta = np.arcsin(numerical_aperture) # defines physical aperture size
# obj_f = 1.*mm  # still dont know what this is supposed to be
obj_f = height
tube_f = magnification * obj_f
radius_of_integration_sphere = obj_f

dipole_orientation = parameters['general']['dipole_orientation']

print('dipole orientation set to {}'.format(dipole_orientation))

## values required for scattered field calculation of sphere
sphere_points = fib.fib_alg_k_filter(
    num_points=lens_points, 
    max_ang=max_theta
    )

cart_points_on_sph = fib.sphere_to_cart(
    sphere_points[:,0],
    sphere_points[:,1],
    obj_f*np.ones(np.shape(sphere_points[:,0]))
    )

r_hat = cart_points_on_sph / radius_of_integration_sphere

refracted_coords = diffi.refract_sph_coords(sphere_points, obj_f, tube_f)
# print('sphere_points= ',sphere_points)
# print('refracted_coords= ',refracted_coords)

## for 1D slice
obv_slice = np.linspace(-sensor_size/2, sensor_size/2, resolution)
obv_slice = np.vstack((obv_slice, np.zeros(obv_slice.shape))).T

#######################################################################

#######################################################################
## functions
#######################################################################
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
def twoD_Gaussian(
    (x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset
    ):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()

def add_colorbar(im, aspect=20, pad_fraction=0.5, **kwargs):
    """Add a vertical color bar to an image plot."""
    divider = axes_grid1.make_axes_locatable(im.axes)
    width = axes_grid1.axes_size.AxesY(im.axes, aspect=1./aspect)
    pad = axes_grid1.axes_size.Fraction(pad_fraction, width)
    current_ax = plt.gca()
    cax = divider.append_axes("right", size=width, pad=pad)
    plt.sca(current_ax)
    return im.axes.figure.colorbar(im, cax=cax, **kwargs)


#######################################################################
## Begin script
#######################################################################

# plot_mode = raw_input('want to see field plots? (\'y\') \n')
# plot_cents = raw_input('Plot centroid offsets? (\'y\') \n')


##  Load separations and calculated aplitudes 
# data_file_name = file_title + '.out'
# amplitude_data = np.loadtxt(data_file_name)
seps = osc.seperation_columb(
    parameters['general']['min_sep'],
    parameters['general']['max_sep'],
    parameters['general']['sep_step_size']
    )  # column in 2d

min_energy = parameters['general']['min_energy']
max_energy = parameters['general']['max_energy']
freq_step_size = parameters['general']['freq_step_size']
energies = np.arange(min_energy, max_energy, freq_step_size)
# energies = np.arange(2.6, 3, 1)

sprectra_per_seperation = np.zeros((seps.shape[0], energies.shape[0]))

#######################################################################
## assign dipole positons in scattering plane
#######################################################################

d1_pos_unit = np.array([-1, 0, 0])
d2_pos_unit = np.array([ 1, 0, 0])

## fix plasmon offset from center
# # d1_from_center = sensor_size/8
# d1_from_center = 0
# d1_pos = d1_from_center * d1_pos_unit * np.ones(np.shape(seps))
# d2_pos = d2_pos_unit * (seps + d1_pos)

## center fluorophore
d2_from_center = 0
d2_pos = d2_from_center * d2_pos_unit * np.ones(np.shape(seps))
d1_pos = d1_pos_unit * (seps + d2_pos)

# d2_pos = 0 * seps
# center plasmon
# d1_from_center = 0
# d1_pos = d1_from_center * d1_pos_unit * np.ones(np.shape(seps))
# d2_pos = d2_pos_unit * (seps + d1_pos)

#######################################################################
## Define variables used for simulated measurement grid
#######################################################################

xs = np.linspace(-sensor_size/2, sensor_size/2, resolution)
ys = np.linspace(-sensor_size/2, sensor_size/2, resolution)
X, Y = np.meshgrid(xs, ys)
X_normed = X / (nm)  # assumed definition of X and Y with nm
Y_normed = Y / (nm)

## scale marking dot size to plot domain for field plots
dot_size = (0.025 * sensor_size/(2*nm)) 

## initialize array to store Poynting vec (S_T) at all seperations
slice_meta_data = np.zeros((2,seps.shape[0],resolution))
gaussian_slice_meta_data = np.zeros((2,seps.shape[0],resolution))
## and again to compare simplified asymptotic expression for S_T
twoD_meta_data = np.zeros((2,seps.shape[0],resolution**2))

# if what_do == 'g':  # calculate fields and fit gaussians
    ## initialize vector for storing maximum intensities
max_int = []
diff_max_int = []
no_fit = []
diff_no_fit = []

centroids = np.hstack((seps, np.zeros((seps.shape[0], 3))))

mislocalization = np.zeros((seps.shape[0], 1))
mislocalization_map = np.zeros((energies.shape[0], seps.shape[0]))

fit_params = np.hstack((seps, np.zeros((seps.shape[0], 7))))
diff_fit_params = np.hstack((seps, np.zeros((seps.shape[0], 7))))

ini_fluo_pos = d2_pos[0][0]/nm  
sca_ini_gue = [1, 
    ini_fluo_pos, 0.0, 
    4.06241329e+04, 5.29040320e+04, 
    1e-11, 1e-01]  # used at sensor_size = 1e6*nm

diff_ini_gue = [1, 
        ini_fluo_pos, 1.87573764e-05,
        1.20641229e+02, 1.03712147e+02,
        1.04596226e-07, 2.22083446e-03]

#######################################################################
## Begin loop through seperations
#######################################################################
for j in range(energies.shape[0]):
    q1, phase1 = osc.plas_l1_osc(dipole_orientation, energies[j])
    q2, phase2 = osc.fluo_osc(dipole_orientation, energies[j])

    print('beggining seperation loop at {} eV'.format(energies[j]))

    for i in range(seps.shape[0]):
        x1 = d1_pos[i]
        x2 = d2_pos[i]
       

    #######################################################################
        ## Compute diffracted fields on GRS and energy flux on CCD
    #######################################################################
        # print(
        #     'Calculating fields at sep = {}'.format(seps.ravel()[i]))
        diff_r_d1 = cart_points_on_sph - x1
        diff_r_d2 = cart_points_on_sph - x2

        di_E_field_on_sph = (
            fi.E_dipole_complex(
                w=energies[j]/hbar, 
                dipole_magnitude=q1[i], 
                dipole_phase=phase1[i], 
                dipole_ori_unit=dipole_orientation, 
                r=diff_r_d1
                )
            +
            fi.E_dipole_complex(
                w=energies[j]/hbar, 
                dipole_magnitude=q2[i], 
                dipole_phase=phase2[i], 
                dipole_ori_unit=dipole_orientation, 
                r=diff_r_d2
                )
            )

        di_B_field_on_sph = (
            fi.B_dipole_complex(
                w=energies[j]/hbar, 
                dipole_magnitude=q1[i], 
                dipole_phase=phase1[i], 
                dipole_ori_unit=dipole_orientation, 
                r=diff_r_d1
                )
            +
            fi.B_dipole_complex(
                w=energies[j]/hbar, 
                dipole_magnitude=q2[i], 
                dipole_phase=phase2[i], 
                dipole_ori_unit=dipole_orientation, 
                r=diff_r_d2
                )
            )

        s = fi.S_complex(
                di_E_field_on_sph,
                di_B_field_on_sph
                )
        S_time_avrg = np.real(s)
        
        n_dot_S = np.sum(S_time_avrg * r_hat, axis=-1)


        #### Perform integral like in diffraction
        number_of_lens_points = cart_points_on_sph.shape[0]
        lens_surface_area = (
            2*np.pi  
            * 
            radius_of_integration_sphere**2. 
            )
        area_per_point = lens_surface_area/number_of_lens_points

        ## Perform integral
        total_power_radiated = np.sum(n_dot_S*area_per_point)

        ## Divide by incident field intensity to normamalize, units of per area
        cross_section = total_power_radiated / inc_field_intens


        # print(total_power_radiated)
        sprectra_per_seperation[i, j] = cross_section


## iniciate array to hold spectral peaks and half-max energies
peaks_and_hm_energies = np.zeros((seps.shape[0], 3))

## find max values and fwhm values
for i in range(seps.shape[0]):

    peak_pos = np.argmax(sprectra_per_seperation[i])
    peak_energy = energies[peak_pos] 

    half_max = sprectra_per_seperation[i,peak_pos]/2.
    
    first_half_spec = sprectra_per_seperation[i,:peak_pos]
    second_half_spec = sprectra_per_seperation[i,peak_pos:]
    
    first_hm_pos = np.argmin(
        (first_half_spec - half_max)**2.
        )
    second_hm_pos = np.argmin(
        (second_half_spec - half_max)**2.
        )

    second_half_energies = energies[peak_pos:]

    peaks_and_hm_energies[i] = [energies[first_hm_pos], peak_energy, second_half_energies[second_hm_pos]]




# spectra_figure = plt.figure()
# for i in range(seps.shape[0]):
    
#     plot = plt.plot(energies, sprectra_per_seperation[i],
#         label='d = {:06.0f} nm'.format(seps[i][0] /nm)
#         )
# # spectra_plot = plt.plot(drive_hbar_omega_range, 
# #     (
# #         sprectra_per_seperation[-1] *np.max(sprectra_per_seperation[0])
# #         / np.max(sprectra_per_seperation[-1])
# #         ),
# #     label='d = {:06.0f} nm, normalized'.format(seps[-1][0] /nm)
# #     )

# ## Mark plasmon resonance
# plas_resonance = parameters['plasmon']['fit_hbar_w0']
# plt.plot(
#     [plas_resonance,plas_resonance],
#     [0,np.max(sprectra_per_seperation[0])],
#     color='black',ls='--')

# ## Mark fluo resonance
# fluo_resonance = parameters['fluorophore']['res_energy']
# plt.plot(
#     [fluo_resonance,fluo_resonance],
#     [0,np.max(sprectra_per_seperation[0])],
#     color='black',ls='--')

# plt.legend()



# ## peak and half max plot
# plt.figure()
# plt.plot(seps.ravel() ,peaks_and_hm_energies[:,0], ls='--')
# plt.plot(seps.ravel() ,peaks_and_hm_energies[:,1])
# plt.plot(seps.ravel() ,peaks_and_hm_energies[:,2], ls='--')


# multipage_figure_file_name = raw_input('Name for figure pdf? ')
# if multipage_figure_file_name != '':
#     multipage(multipage_figure_file_name, figs=None, dpi=200)

# plt.show()
# print(sprectra_per_seperation)

########################
#### Save data
if dipole_orientation == [1,0,0]:
    orientation_letter = 'x'
elif dipole_orientation == [0,1,0]:
    orientation_letter = 'y'
elif dipole_orientation == [0,0,1]:
    orientation_letter = 'z'

print('nparray([0,energies].shape = ',[[[0],energies]])
print('(seps / nm).shape = ', (seps / nm).shape)
print('sprectra_per_seperation', sprectra_per_seperation.shape)
np.savetxt(
    'spectra_at_seps_{}_or.txt'.format(orientation_letter), 
    np.vstack((
        np.append([0],energies), 
        np.hstack((seps / nm, sprectra_per_seperation)) 
        )),
    header= 'first row: energies (eV), first column: separations (nm), emission crossection (cm^-2)'
    )

np.savetxt(
    'half_max_and_peak_energies_per_sep_{}_or.txt'.format(orientation_letter),
    peaks_and_hm_energies,
    header='left half max, peak energy, right half max (all in eV)'
    )




