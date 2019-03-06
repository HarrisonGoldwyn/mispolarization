"""
This file provides classes and functions for modeling and fitting the
diffraction-limited images produced by a dipole emitter coupled to a
plasmonic nanorod, modeled as a polarizable point dipole with
polarizability of a prolate ellipsoid in the modified long wavelength
approximation.
-----------------------------------------------------------------------


-----------
Patch notes
-----------

02/06/19:
    This file was renamed from
    'fitting_misLocalization_adding_noise_to_modeled_images__011619v11'.
    Currently hardcoded for a certain parameter file and fit
    parameters. This will be changed in the future, but for now I just
    want to get stuff done.

02/07/19:
    Hardcoded parameter file changed to vacuum values fit from BEM spectra
    with built in.

02/21/19:
    Removed hardcoded dependence on .yaml file for plasmon and fluo
    parameters. Leaving some bits that depend on the 'general'
    parameters for now.

03/06/19:
    Added functionality to remove interference at initialization of
    'MolCoupNanoRodExp' instance.

    Sometime before this I also added functionality to isolate mode
    effects... That happened in the last month.

----
TODO
----

- Want to eliminate hardcoded dependence to .yaml files.
    - Should load physical constants from scipy, yaml file is pretty useless.

    """
from __future__ import print_function
from __future__ import division

import pdb
import sys
import os
import numpy as np
import scipy.optimize as opt
from scipy import interpolate
import scipy.io as sio
import scipy.special as spf
import yaml

## import diffrantion integral solver from Optics folder
project_path = ('/Users/chair/Documents/Academia/SuperRes/Biteen_colab/'
    +
    'Mispolarization/python/gitted'
    )

optics_path = project_path + '/Optics'

sys.path.append(optics_path)
import diffraction_int as diffi
import fibonacci as fib



## Read parameter file to obtain fields
parameter_files_path = (
    project_path + '/parameter_files')

curly_yaml_file_name = '/curly_nrod_water_JC.yaml'
parameters = yaml.load(open(parameter_files_path+curly_yaml_file_name, 'r'))
print('reading parameters from {}'.format(
    parameter_files_path+curly_yaml_file_name
    )
)


## plotting stuff
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['text.usetex'] = True
mpl.rcParams["lines.linewidth"]

## colorbar stuff
from mpl_toolkits import axes_grid1


## analytic image fields
import anal_foc_diff_fields as afi

## solution to coupled dipole problem
modules_path = project_path + '/solving_problems/modules'
sys.path.append(modules_path)
import coupled_dipoles as cp

txt_file_path = project_path + '/solving_problems/txt_files'

## Import physical constants
phys_const_file_name = '/physical_constants.yaml'
opened_constant_file = open(
    parameter_files_path+phys_const_file_name,
    'r')

constants = yaml.load(opened_constant_file)
e = constants['physical_constants']['e']
c = constants['physical_constants']['c']  # charge of electron in statcoloumbs
hbar = constants['physical_constants']['hbar']
m_per_nm = constants['physical_constants']['nm']
n_a = constants['physical_constants']['nA']   # Avogadro's number
# Z_o = 376.7303 # impedence of free space in ohms (SI)

## STOPPED HERE 02/21/19 5:00 PM. Removing parameters dependance
## System background
n_b = parameters['general']['background_ref_index']
eps_b = n_b**2.


# a = parameters['plasmon']['radius']




#######################################################################
## Optics stuff.
sensor_size = parameters['optics']['sensor_size']*m_per_nm
# height = 2*mm  # also defines objective lens focal length
# height = parameters['optics']['obj_f_len']
resolution = parameters['optics']['sensor_pts']  # image grid resolution
## Build image sensor
eye = diffi.observation_points(
    x_min= -sensor_size/2,
    x_max= sensor_size/2,
    y_min= -sensor_size/2,
    y_max= sensor_size/2,
    points= resolution
    )

## Experimental parameters
magnification = parameters['optics']['magnification']
numerical_aperture = parameters['optics']['numerical_aperture']
max_theta = np.arcsin(numerical_aperture) # defines physical aperture size

## numerical parameters for calculation of scattered field
lens_points = parameters['optics']['lens_points']

# obj_f = 1.*mm  # still dont know what this is supposed to be
obj_f = parameters['optics']['obj_f_len']

tube_f = magnification * obj_f

## calculate dipole magnitudes

# drive_hbar_omega = parameters['general']['drive_energy']
    ## rod long mode max at 1.8578957289256757 eV
# omega_drive = drive_hbar_omega/hbar  # driving frequency





class DipoleProperties(object):
    """ Will eventually call parameter file as argument, currently (02/07/19)
        just loads relevant values from hardcoded paths. ew.
        """


    def __init__(self,
        eps_inf=parameters['plasmon']['fit_eps_inf'],
        hbar_omega_plasma=parameters['plasmon']['fit_hbar_wp'],
        hbar_gamma_drude=parameters['plasmon']['fit_hbar_gamma'],
        a_long_in_nm=parameters['plasmon']['fit_a1'],
        a_short_in_nm=parameters['plasmon']['fit_a2'],
        eps_b=parameters['general']['background_ref_index']**2.0,
        fluo_ext_coef=parameters['fluorophore']['extinction_coeff'],
        fluo_mass_hbar_gamma=parameters['fluorophore']['mass_gamma'],
        fluo_nr_hbar_gamma=parameters['fluorophore']['test_gamma'],
        fluo_quench_region_nm=10,
        isolate_mode=None,
        drive_energy_eV=parameters['general']['drive_energy'],
        ):

        self.drive_energy_eV = drive_energy_eV
        self.eps_inf = eps_inf
        self.omega_plasma = hbar_omega_plasma / hbar
        self.gamma_drude = hbar_gamma_drude / hbar
        self.a_long_meters = a_long_in_nm * m_per_nm
        self.a_short_meters = a_short_in_nm * m_per_nm

        # self.fit_result_params = [
        #     ## eps_inf, hbar*omega_p, hbar*gamma_nr, eps_b
        #     ## (not used as fit param), a_x, a_yz
        #     eps_inf, # parameters['plasmon']['fit_eps_inf'],
        #     hbar_omega_plasma / hbar, # parameters['plasmon']['fit_hbar_wp']/hbar,
        #     hbar_gamma_drude / hbar,# parameters['plasmon']['fit_hbar_gamma']/hbar,
        #     a_long_in_nm * m_per_nm, # parameters['plasmon']['fit_a1']*m_per_nm,
        #     a_short_in_nm * m_per_nm, # parameters['plasmon']['fit_a2']*m_per_nm
        #     ]

        self.eps_b = eps_b

        ## hardcoded region around nanoparticle to through out results because
        ## dipole approximation at small proximities
        self.fluo_quench_range = fluo_quench_region_nm

        self.alpha0_diag_dyad = cp.sparse_polarizability_tensor(
            mass=cp.fluorophore_mass(
                ext_coef=fluo_ext_coef, # parameters['fluorophore']['extinction_coeff'],
                gamma=fluo_mass_hbar_gamma/hbar, # parameters['fluorophore']['mass_gamma']/hbar
                ),
            w_res=drive_energy_eV/hbar,
            w=drive_energy_eV/hbar,
            gamma_nr=fluo_nr_hbar_gamma/hbar, # parameters['fluorophore']['test_gamma']/hbar,
            a=0,
            eps_inf=1,
            ebs_b=1
            )

        self.alpha1_diag_dyad = (
            cp.sparse_ret_prolate_spheroid_polarizability_Drude(
                drive_energy_eV/hbar,
                self.eps_inf,
                self.omega_plasma,
                self.gamma_drude,
                eps_b,
                self.a_long_meters,
                self.a_short_meters,
                isolate_mode=isolate_mode)
            )


class BeamSplitter(object):

    def __init__(self):
        pass

    def powers_and_angels(self,E):
        drive_I = np.abs(parameters['general']['drive_amp'])**2.

        normed_Ix = np.abs(E[0])**2. / drive_I
        normed_Iy = np.abs(E[1])**2. / drive_I

        Px_per_drive_I = np.sum(normed_Ix,axis=-1) / sensor_size**2.
        Py_per_drive_I = np.sum(normed_Iy,axis=-1) / sensor_size**2.


        angles = np.arctan(Py_per_drive_I**0.5/Px_per_drive_I**0.5)
        return [angles, Px_per_drive_I, Py_per_drive_I]

    def powers_and_angels_no_interf(self,E1,E2):
        drive_I = np.abs(parameters['general']['drive_amp'])**2.

        normed_Ix = (np.abs(E1[0])**2. + np.abs(E2[0])**2.) / drive_I
        normed_Iy = (np.abs(E1[1])**2. + np.abs(E2[1])**2.) / drive_I

        Px_per_drive_I = np.sum(normed_Ix,axis=-1) / sensor_size**2.
        Py_per_drive_I = np.sum(normed_Iy,axis=-1) / sensor_size**2.


        angles = np.arctan(Py_per_drive_I**0.5/Px_per_drive_I**0.5)
        return [angles, Px_per_drive_I, Py_per_drive_I]


class FittingTools(object):

    def __init__(self, obs_points=None):
        """
        Args:
            obs_points: 3 element list (in legacy format of eye),
            in units of nm.

                obs_points[0]: list of points as rows
                obs_points[1]: meshed X array
                obs_points[2]: meshed Y array
        """
        if obs_points==None:
            self.obs_points = eye
        else:
            self.obs_points = obs_points



    def twoD_Gaussian(self,
        X, ## tuple of meshed (x,y) values
        amplitude,
        xo,
        yo,
        sigma_x,
        sigma_y,
        theta,
        offset,
        ):

        xo = float(xo)
        yo = float(yo)
        a = (
            (np.cos(theta)**2)/(2*sigma_x**2)
            +
            (np.sin(theta)**2)/(2*sigma_y**2)
            )
        b = (
            -(np.sin(2*theta))/(4*sigma_x**2)
            +
            (np.sin(2*theta))/(4*sigma_y**2)
            )
        c = (
            (np.sin(theta)**2)/(2*sigma_x**2)
            +
            (np.cos(theta)**2)/(2*sigma_y**2)
            )
        g = (
            offset
            +
            amplitude*np.exp( - (a*((X[0]-xo)**2) + 2*b*(X[0]-xo)*(X[1]-yo)
            +
            c*((X[1]-yo)**2)))
            )
        return g.ravel()

    def misloc_data_minus_model(
        self,
        fit_params,
        *normed_raveled_image_data,
        ):
        ''' fit gaussian to data '''
        gaus = self.twoD_Gaussian(
            (self.obs_points[1]/m_per_nm, self.obs_points[2]/m_per_nm),
            *fit_params ## ( A, xo, yo, sigma_x, sigma_y, theta, offset)
            )

        return gaus - normed_raveled_image_data

    def calculate_max_xy(self, images):
        ## calculate index of maximum in each image.
        apparent_centroids_idx = images.argmax(axis=-1)
        ## define locations for each maximum in physical coordinate system

        x_cen = (self.obs_points[1]/m_per_nm).ravel()[apparent_centroids_idx]
        y_cen = (self.obs_points[2]/m_per_nm).ravel()[apparent_centroids_idx]

        return [x_cen,y_cen]

    def calculate_apparent_centroids(self, images):
        """ calculate index of maximum in each image. """
        num_of_images = images.shape[0]

        apparent_centroids_xy = np.zeros((num_of_images,2))

        max_positions = self.calculate_max_xy(images)

        for i in np.arange(num_of_images):
            x0 = max_positions[0][i]
            y0 = max_positions[1][i]
            params0 = (1,x0,y0,100, 100, 0,0)
            args=tuple(images[i]/np.max(images[i]))
            fit_gaussian = opt.least_squares(self.misloc_data_minus_model, params0, args=args)
            resulting_fit_params = fit_gaussian['x']
            fit_result = self.twoD_Gaussian(
                (self.obs_points[1]/m_per_nm, self.obs_points[2]/m_per_nm), ## tuple of meshed (x,y) values
                *resulting_fit_params
                )
            centroid_xy = resulting_fit_params[1:3]
            apparent_centroids_xy[i] = centroid_xy
        ## define locations for each maximum in physical coordinate system

        return apparent_centroids_xy.T  ## returns [x_cen(s), y_cen(s)]

    def image_from_E(self, E):
        drive_I = np.abs(parameters['general']['drive_amp'])**2.

        normed_I = np.sum(np.abs(E)**2.,axis=0) / drive_I

        return normed_I

class PlottingStuff(DipoleProperties):

    def __init__(self,
        isolate_mode=None,
        drive_energy_eV=parameters['general']['drive_energy'],
        ):
        ## Establish dipole properties as atributes for reference in plotting
        ## functions.
        DipoleProperties.__init__(self,
            isolate_mode=isolate_mode,
            drive_energy_eV=drive_energy_eV,
            )

    def connectpoints(self, cen_x, cen_y, mol_x, mol_y, p, ax=None, zorder=1):
        x1, x2 = mol_x[p], cen_x[p]
        y1, y2 = mol_y[p], cen_y[p]
        if ax == None:
            plt.plot([x1,x2],[y1,y2],'k-', linewidth=.3, zorder=zorder)
        else:
            ax.plot([x1,x2],[y1,y2],'k-', linewidth=.3, zorder=zorder)

    def scatter_centroids_wLine(self, x_mol_loc, y_mol_loc, appar_cents, ax=None):

    #     trial_images = image_from_E(E)
    #     appar_cents = calculate_apparent_centroids(trial_images)

        x, y = appar_cents

    #     ## This part doesnt work right now
    #     el_a = 19
    #     el_c = 67
    #     quel_a = el_a + 10
    #     quel_c = el_c + 10
    #     pt_is_in_ellip = np.ones(x.shape, dtype=bool)
    # #     for i in np.arange(x.shape[0]):
    # #         if (x[i]**2./quel_a**2. +  y[i]**2./quel_c**2.) < 1:
    # #             pt_is_in_ellip[i] = False
    #     ####

        x_plot = x
        y_plot = y

        if ax == None:
            plt.figure(dpi=300)
            for i in np.arange(x_plot.shape[0]):
                self.connectpoints(x_plot, y_plot, x_mol_loc, y_mol_loc, i,zorder=2)

            plt.scatter(x_plot, y_plot, s=10, c='Red', zorder=3)
            plt.tight_layout()


        else:
            for i in np.arange(x_plot.shape[0]):
                self.connectpoints(
                    x_plot, y_plot, x_mol_loc, y_mol_loc, i, ax,zorder=2
                    )
            ax.scatter(x_plot, y_plot, s=10, c='Red', zorder=3)
            return ax
    #     plt.xlim([-2000,2000])
    #     plt.ylim([-2000,2000])


    # In[41]:

    def quiver_plot(
        self,
        x_plot,
        y_plot,
        angles,
        plot_limits=[-25,550],
        title=r'Apparent pol. per mol. pos.',
        true_mol_angle=None,
        nanorod_angle=0,
        given_ax=None
        ):

        if true_mol_angle is None:
            true_mol_angle = angles

        self.el_a = self.a_long_meters / m_per_nm
        self.el_c = self.a_short_meters / m_per_nm

        # self.fluo_quench_range

        # quel_a = el_a + self.fluo_quench_range
        # quel_c = el_c + self.fluo_quench_range
        pt_is_in_ellip = np.ones(x_plot.shape, dtype=bool)
    #     for i in np.arange(x_plot.shape[0]):
    #         if (x_plot[i]**2./quel_a**2. +  y_plot[i]**2./quel_c**2.) < 1:
    #             pt_is_in_ellip[i] = False

        x_plot = x_plot[pt_is_in_ellip]
        y_plot = y_plot[pt_is_in_ellip]
        angles = angles[pt_is_in_ellip]

        if given_ax==None:
            fig, (ax0, ax_cbar) = plt.subplots(
                nrows=1,ncols=2, figsize=(3.25,3), dpi=300,
                gridspec_kw = {'width_ratios':[6, 0.5]}
                )
        else:
            ax0 = given_ax

        cmap = mpl.cm.nipy_spectral

        ## Mark molecule locations
        scat_tr = ax0.scatter(x_plot, y_plot, s=3,
                        color='black',
        #                    cmap='inferno',
        #                    clim = [0, np.pi/2],
    #                        width=0.005,
    #                        scale=20,
                #            scale_units='width',
    #                        pivot='mid',
        #                   linewidth=100.,
    #                       headaxislength=0.0,
    #                       headlength=0.0
                          )
        ## mark true orientation
        quiv_tr = ax0.quiver(
            x_plot, y_plot, np.cos(true_mol_angle),np.sin(true_mol_angle),
            color='black',
        #     cmap='inferno',
        #     clim = [0, np.pi/2],
            width=0.005,
            scale=15,
              scale_units='width',
            pivot='mid',
        #   linewidth=100.,
            headaxislength=0.0,
            headlength=0.0
            )

        ## Mark apparent orientation
        quiv_ap = ax0.quiver(
            x_plot,
            y_plot,
            np.cos(angles),
            np.sin(angles),
            angles,
            cmap=cmap,
            clim = [0, np.pi/2],
            width=0.01,
            scale=12,
            scale_units='width',
            pivot='mid',
            zorder=4,
            headaxislength=2.5,
            headlength=2.5,
            headwidth=2.5,
            )

        ax0.axis('equal')
        ax0.set_xlim(plot_limits)
        ax0.set_ylim(plot_limits)
        ax0.set_title(title)
        ax0.set_xlabel('x [nm]')
        ax0.set_ylabel('y [nm]')

        norm = mpl.colors.Normalize(vmin=0, vmax=np.pi/2)

        # Build colorbar if building single Axes figure
        if given_ax==None:
            cb1 = mpl.colorbar.ColorbarBase(ax_cbar, cmap=cmap,
                                            norm=norm,
                                            orientation='vertical')
            cb1.set_label(r'observed angle $\phi$')

            cb1.set_ticks([0, np.pi/8, np.pi/4, np.pi/8 * 3, np.pi/2])
            cb1.set_ticklabels(
                [r'$0$', r'$22.5$',r'$45$',r'$67.5$',r'$90$']
                )
        else: # Don't build colorbar
            pass

        if nanorod_angle == np.pi/2:
            #### Draw rod
            circle = mpl.patches.Circle((0, 24), 20, facecolor='Gold',
                        edgecolor='Black', linewidth=0)
            bot_circle = mpl.patches.Circle((0, -24), 20, facecolor='Gold',
                        edgecolor='Black', linewidth=0)
            rect = mpl.patches.Rectangle(
                (-20,-24), 40, 48, angle=0.0, facecolor='Gold',
                edgecolor='Black', linewidth=0)

            ax0.add_patch(circle)
            ax0.add_patch(rect)
            ax0.add_patch(bot_circle)

        ## model ellipsoid
        ellip = mpl.patches.Ellipse(
            (0,0), 2*self.el_a, 2*self.el_c, angle=nanorod_angle*180/np.pi,
            fill=False, edgecolor='Black',linestyle='--')
        ax0.add_patch(ellip)

        quiver_axis_handle = ax0
        return [quiver_axis_handle]

    def calculate_mislocalization_magnitude(self, x_cen, y_cen, x_mol, y_mol):
        misloc = ( (x_cen-x_mol)**2. + (y_cen-y_mol)**2. )**(0.5)
        return misloc


class CoupledDipoles(PlottingStuff, FittingTools):

    ## Q: do I need to manually call 'PlottingStuff.__init__'?
    def __init__(self,
        obs_points=None,
        isolate_mode=None,
        drive_energy_eV=parameters['general']['drive_energy'],
        ):
        """
        PlottingStuff.__init__(): Really just calls DipoleProperties.__init__()
            which initializes polarizabilities and maybe some other stuff

        FittingTools.__init__(obs_points): defines obs points or assigns
            default
         """
        PlottingStuff.__init__(self,
            isolate_mode=isolate_mode,
            drive_energy_eV=drive_energy_eV,
            )
        FittingTools.__init__(self, obs_points)

    def mb_p_fields(self, dipole_mag_array, dipole_coordinate_array):
        ''' As of 081418,fixing: currently only treats dipole at origin.'''
        p = dipole_mag_array
    #     print('Inside mb_p_fields, p= ',p)
        bfx = dipole_coordinate_array

        v_rel_obs_x_pts = (self.obs_points[1].ravel()[:,None] - bfx.T[0]).T
        v_rel_obs_y_pts = (self.obs_points[2].ravel()[:,None] - bfx.T[1]).T

        px_fields = np.asarray(
            afi.E_field(
                0,
                v_rel_obs_x_pts,
                v_rel_obs_y_pts,
                (self.drive_energy_eV/hbar)*np.sqrt(self.eps_b)/c
                )
            )
        py_fields = np.asarray(
            afi.E_field(
                np.pi/2,
                v_rel_obs_x_pts,
                v_rel_obs_y_pts,
                (self.drive_energy_eV/hbar)*np.sqrt(self.eps_b)/c
                )
            )
        pz_fields = np.zeros(py_fields.shape)
    #     print('px_fields.shape=',px_fields.shape)
    #     print('p.shape=',p.shape)
        ## returns [Ex, Ey, Ez] for dipoles oriented along cart units

        Ex = (
            p[:,0,None]*px_fields[0]
            +
            p[:,1,None]*py_fields[0]
            +
            p[:,2,None]*pz_fields[0]
            )
        Ey = (
            p[:,0,None]*px_fields[1]
            +
            p[:,1,None]*py_fields[1]
            +
            p[:,2,None]*pz_fields[1]
            )
        Ez = (
            p[:,0,None]*px_fields[2]
            +
            p[:,1,None]*py_fields[2]
            +
            p[:,2,None]*pz_fields[2]
            )

        return np.array([Ex,Ey,Ez])

    def dipole_fields(self, locations, mol_angle=0, plas_angle=np.pi/2):
        d = locations*m_per_nm
        p0, p1 = cp.dipole_mags_gened(
            mol_angle,
            plas_angle,
            d_col=d,
            E_d_angle=None,
            alpha0_diag=self.alpha0_diag_dyad,
            alpha1_diag=self.alpha1_diag_dyad,
            )
        mol_E = self.mb_p_fields(
            dipole_mag_array=p0,
            dipole_coordinate_array=d,
            )
        plas_E = self.mb_p_fields(
            dipole_mag_array=p1,
            dipole_coordinate_array=np.zeros(d.shape),
            )

        # p0_unc, = cp.uncoupled_p0(mol_angle=0, d_col=d[0,None], E_d_angle=None)
        p0_unc, = cp.uncoupled_p0(
            mol_angle,
            E_d_angle=None,
            drive_hbar_w=self.drive_energy_eV,
            )
    #     print('p0.shape = ',p0.shape)
    #     print('p1.shape = ',p1.shape)
    #     print('p0_unc.shape = ',p0_unc.shape)
    #     p0_unc_E = self.mb_p_fields(dipole_mag_array=p0_unc[None,:], dipole_coordinate_array=np.zeros(d[0][None,:].shape))
        if type(mol_angle)==np.ndarray and mol_angle.shape[0]>1:
            p0_unc_E = self.mb_p_fields(
                dipole_mag_array=p0_unc,
                dipole_coordinate_array=d
                )
        elif (type(mol_angle) == int or
              type(mol_angle) == float or
              type(mol_angle) == np.float64 or
              (type(mol_angle) == np.ndarray and mol_angle.shape[0]==1)
              ):
            p0_unc_E = self.mb_p_fields(
                dipole_mag_array=p0_unc[None,:],
                dipole_coordinate_array=d,
                )
#         print(type(mol_angle))
        return [mol_E, plas_E, p0_unc_E, p0, p1]


# In[44]:

class MolCoupNanoRodExp(CoupledDipoles, BeamSplitter):
    ''' Collect focused+diffracted far-field information from molecules
        nearby a nanorod.
        '''

    ## set up inverse mapping from observed -> true angle for signle molecule in the plane.
    saved_mapping = np.loadtxt(txt_file_path+'/obs_pol_vs_true_angle.txt')
    true_ord_angles, obs_ord_angles =  saved_mapping.T
    #from scipy import interpolate
    f = interpolate.interp1d(true_ord_angles,obs_ord_angles)
    f_inv = interpolate.interp1d(
        obs_ord_angles[:251],
        true_ord_angles[:251],
        bounds_error=False,
        fill_value=(0,np.pi/2)
        )

    def __init__(
        self,
        locations,
        mol_angle=0,
        plas_angle=np.pi/2,
        obs_points=None,
        for_fit=False,
        isolate_mode=None,
        drive_energy_eV=parameters['general']['drive_energy'],
        exclude_interference=False,
        ):

        CoupledDipoles.__init__(self,
            obs_points,
            isolate_mode,
            drive_energy_eV,
            )

        # Set up instance attributes
        self.exclude_interference = exclude_interference
        self.mol_locations = locations
        self.mol_angles = mol_angle
        self.rod_angle = plas_angle

        # Filter out molecules in region of fluorescence quenching
        self.el_a = self.a_long_meters / m_per_nm
        self.el_c = self.a_short_meters / m_per_nm
        self.quel_a = self.el_a + self.fluo_quench_range ## define quenching region
        self.quel_c = self.el_c + self.fluo_quench_range
        self.input_x_mol = locations[:,0]
        self.input_y_mol = locations[:,1]
        self.pt_is_in_ellip = self.mol_too_close()
        ## select molecules outside region,
        if for_fit==False:
            self.mol_locations = locations[self.pt_is_in_ellip]
            ## select molecule angles if listed per molecule,
            if type(mol_angle)==np.ndarray and mol_angle.shape[0]>1:
                self.mol_angles = mol_angle[self.pt_is_in_ellip]
            else: self.mol_angles = mol_angle
        elif for_fit==True:
            self.mol_locations = locations
            self.mol_angles = mol_angle

        # Automatically calculate fields with coupled dipoles upon
        # instance initialization.
        (
            self.mol_E,
            self.plas_E,
            self.p0_unc_E,
            self.p0,
            self.p1
            ) = self.dipole_fields(
                self.mol_locations,
                self.mol_angles,
                self.rod_angle
                )

        # Calcualte images
        self.trial_images = self.image_from_E(self.mol_E + self.plas_E )

        # Calculate plot domain from molecule locations
        self.default_plot_limits = [
            (
                np.min(self.mol_locations)
                -
                (
                    (
                        np.max(self.mol_locations)
                        -
                        np.min(self.mol_locations))*.1
                    )
                ),
            (
                np.max(self.mol_locations)
                +
                (
                    (
                        np.max(self.mol_locations)
                        -
                        np.min(self.mol_locations))*.1
                    )
                ),
            ]

    def calculate_localization(self, save_fields=True):
        """ """
        self.appar_cents = self.calculate_apparent_centroids(self.trial_images)
        self.x_cen, self.y_cen = self.appar_cents
        if save_fields == False:
            del self.mol_E
            del self.plas_E


    def calculate_polarization(self):
        """ Calculate polarization with beam splitter """

        # Calculate fields and angles and assign as instance attribute
        if hasattr(self, 'mol_E') and hasattr(self, 'plas_E'):
            if self.exclude_interference == False:
                self.angles, self.Px_per_drive_I, self.Py_per_drive_I = (
                    self.powers_and_angels(
                        self.mol_E + self.plas_E
                        )
                    )
            # For exclusion of interference, fields must be input
            # seperately into funtion 'powers_and_angels_no_interf'.
            elif self.exclude_interference == True:
                self.angles, self.Px_per_drive_I, self.Py_per_drive_I = (
                    self.powers_and_angels_no_interf(
                        self.mol_E,
                        self.plas_E
                        )
                    )

        # Simulation instance will have field namd differently and
        # must be transposed.
        elif hasattr(self, 'bem_E'):
            # Can not extract interference from simulation
            self.angles, self.Px_per_drive_I, self.Py_per_drive_I = (
                    self.powers_and_angels(
                        np.transpose(self.bem_E, (2,0,1))
                        )
                )
        self.mispol_angle= MolCoupNanoRodExp.f_inv(self.angles)


    def mol_too_close(self):
        '''Returns molecule locations that are outside the fluorescence
           quenching zone, defined as 10 nm from surface of fit spheroid
            '''
        rotated_x = (
            np.cos(self.rod_angle)*self.input_x_mol
            + np.sin(self.rod_angle)*self.input_y_mol
            )
        rotated_y = (
            -np.sin(self.rod_angle)*self.input_x_mol
            + np.cos(self.rod_angle)*self.input_y_mol
            )

        long_quench_radius = self.quel_a
        short_quench_radius = self.quel_c

        rotated_ellip_eq = (
            rotated_x**2./long_quench_radius**2
            + rotated_y**2./short_quench_radius**2
            )
        return (rotated_ellip_eq > 1)


    def plot_mispol_map(self,
        plot_limits=None,
        given_ax=None
        ):

        if plot_limits == None: plot_limits = self.default_plot_limits
        if not hasattr(self, 'mispol_angle'):
            self.calculate_polarization()
        quiv_ax, = self.quiver_plot(self.mol_locations[:,0],
            self.mol_locations[:,1],
            self.mispol_angle,
            plot_limits,
            true_mol_angle=self.mol_angles,
            nanorod_angle=self.rod_angle,
            title=r'Split Pol. and Gau. Fit Loc.',
            given_ax=given_ax
            )
        return quiv_ax


    def plot_mispol_map_wMisloc(self, plot_limits=None):
        if not hasattr(self, 'appar_cents'):
            self.calculate_localization()
        if plot_limits == None: plot_limits = self.default_plot_limits
        quiv_ax = self.plot_mispol_map(plot_limits)

        self.scatter_centroids_wLine(self.mol_locations[:,0],
            self.mol_locations[:,1],
            self.appar_cents,
            quiv_ax,
            )
        return quiv_ax

    def plot_mislocalization_magnitude_correlation(self):
        if not hasattr(self, 'appar_cents'):
            self.calculate_localization()
        self.misloc_mag = calculate_mislocalization_magnitude(
            self.x_cen,
            self.y_cen,
            self.mol_locations[:,0],
            self.mol_locations[:,1],
            )

        plt.figure(dpi=300)
        plt.scatter(self.misloc_mag, mispol_angle, s=10, c='Black', zorder=3)
        plt.tight_layout()
        plt.xlabel('Magnitude of mislocalization [nm]')
        plt.ylabel('Apparent angle [deg]')
        plt.yticks([0,  np.pi/8,  np.pi/4, np.pi/8 *3, np.pi/2],
                   ['0','22.5','45','57.5','90'])
        return plt.gca()

    def plot_fields(self, ith_molecule):
        plt.figure(figsize=(3,3),dpi=600)
        plt.pcolor(eye[1]/m_per_nm,eye[2]/m_per_nm,(self.trial_images[ith_molecule,:]).reshape(eye[1].shape))
        plt.colorbar()
        plt.title(r'$|E|^2/|E_\mathrm{inc}|^2$')
        plt.xlabel(r'x [nm]')
        plt.ylabel(r'y [nm]')
#         plt.quiver(self.mol_locations[ith_molecule, 0], self.mol_locations[ith_molecule, 1],
#                    np.cos(self.mol_angles[ith_molecule]),np.sin(self.mol_angles[ith_molecule]),
#                    color='white',pivot='middle')
        return plt.gca()

# ## Fit function for model fields.

# In[87]:

# MolCoupNanoRodExp(locations, mol_angle=0, plas_angle=np.pi/2, x_obv_grid=eye[1], y_obv_grid=eye[2])
class FitModelToData(FittingTools,PlottingStuff):
    ''' Class to contain fitting functions that act on class 'MolCoupNanoRodExp'
    as well as variables that are needed by 'MolCoupNanoRodExp'
    '''
    def __init__(self,
        image_data,
        obs_points=None,
        ini_guess=None,
        isolate_mode=None,
        drive_energy_eV=parameters['general']['drive_energy'],
        ):

        self.mol_angle=0
        self.rod_angle=np.pi/2
        self.image_data=image_data

        self.ini_guess = ini_guess

        FittingTools.__init__(self, obs_points)
        ## pointer to DipoleProperties.__init__()
        PlottingStuff.__init__(self,
            isolate_mode,
            drive_energy_eV)

    def fit_model_to_image_data(self, images=None):
        ## calculate index of maximum in each image,
        ## going to use this for the initial position guess

        if images is None:
            images = self.image_data
        ## initialize array to hold fit results for arbitrary number of images
        num_of_images = images.shape[0]
        self.model_fit_results = np.zeros((num_of_images,3))
        ## Use positions of max intensity as initial guess for molecule
        ## position.
        max_positions = self.calculate_max_xy(images)

        ## Loop through images and fit.
        for i in np.arange(num_of_images):
            ## Establish initial guesses for molecules
            if self.ini_guess is None:
                ini_x = np.round(max_positions[0][i])
                ini_y = np.round(max_positions[1][i])
            else:
                ini_x = np.round(self.ini_guess[i, 0])
                ini_y = np.round(self.ini_guess[i, 1])
            print(
                'initial guess for molecule {} location: ({},{})'.format(
                    i, ini_x, ini_y
                    )
                )
            ## Randomize initial molecule oriantation, maybe do something
            ## smarter later.
            ini_mol_orientation = np.random.random(1)*np.pi/2
            params0 = (ini_x, ini_y, ini_mol_orientation)

            ## Normalize images for fitting.
            a_raveled_normed_image =  images[i]/np.max(images[i])

            ## Place image in tuple as required by `opt.least_squares`.
            tuple_normed_image_data=tuple(a_raveled_normed_image)

            ## Run fit.
            optimized_fit = opt.least_squares(
                self._misloc_data_minus_model, ## residual
                params0, ## initial guesses
                args=tuple_normed_image_data, ## data to fit
                )
            ## Store fit result parameters as class instance attribute.
            self.model_fit_results[i] = optimized_fit['x']

        return self.model_fit_results

    def _misloc_data_minus_model(self, fit_params, *normed_raveled_image_data):
        ''' fit image model to data.
        arguments;
            fit_params = [ini_x, ini_y, ini_mol_orintation]
                'ini_x' : units of nm from plas position
                'ini_y' : units of nm from plas position
                'ini_mol_orintation' : units of radians counter-clock from +x
        '''

        raveled_model = self.raveled_model_of_params(fit_params)
        normed_raveled_model = raveled_model/np.max(raveled_model)

        return normed_raveled_model - normed_raveled_image_data

    def raveled_model_of_params(self, fit_params):
        locations = np.array([[fit_params[0], fit_params[1], 0]])
        exp_instance = MolCoupNanoRodExp(
            locations,
            mol_angle=fit_params[2],
            plas_angle=self.rod_angle,
            obs_points=self.obs_points,
            for_fit=True
            )
        raveled_model = exp_instance.trial_images[0].ravel()
        return raveled_model

    def plot_image_from_params(self, fit_params, ax=None):
        raveled_image = self.raveled_model_of_params(fit_params)
        self.plot_raveled_image(raveled_image,ax)
#         plt.quiver(self.mol_locations[ith_molecule, 0], self.mol_locations[ith_molecule, 1],
#                    np.cos(self.mol_angles[ith_molecule]),np.sin(self.mol_angles[ith_molecule]),
#                    color='white',pivot='middle')

    def plot_raveled_image(self, image, ax=None):
        if ax is None:
            plt.figure(figsize=(3,3),dpi=600)
            plt.pcolor(
                self.obs_points[-2]/m_per_nm,
                self.obs_points[-1]/m_per_nm,
                image.reshape(self.obs_points[-2].shape),
                )
            plt.colorbar()
        else:
            ax.contour(self.obs_points[-2]/m_per_nm,
                self.obs_points[-1]/m_per_nm,
                image.reshape(self.obs_points[-2].shape),
                cmap='Greys',
                linewidths=0.5,
                )
        plt.title(r'$|E|^2/|E_\mathrm{inc}|^2$')
        plt.xlabel(r'x [nm]')
        plt.ylabel(r'y [nm]')
        return plt.gca()

    def plot_fit_results_as_quiver_map(
        self,
        fitted_exp_instance,
        plot_limits=None,
        ):
        '''...'''

        if not hasattr(self, 'model_fit_results'):
            self.fit_model_to_image_data()
        if plot_limits == None:
            plot_limits = fitted_exp_instance.default_plot_limits

        quiv_ax, = self.quiver_plot(
            fitted_exp_instance.mol_locations[:,0],
            fitted_exp_instance.mol_locations[:,1],
            self.model_fit_results[:,2],
            plot_limits,
            true_mol_angle = fitted_exp_instance.mol_angles,
            nanorod_angle = fitted_exp_instance.rod_angle,
            title=r'Model Fit Pol. and Loc.',
            )
        self.scatter_centroids_wLine(fitted_exp_instance.mol_locations[:,0],
                                fitted_exp_instance.mol_locations[:,1],
                                self.model_fit_results[:,:2].T, quiv_ax)
        return quiv_ax

    def plot_contour_fit_over_data(self, image_idx):
            ax = self.plot_raveled_image(self.image_data[image_idx])
            self.plot_image_from_params(self.model_fit_results[image_idx], ax)


# class FitModelToNoisedModel(FitModelToData,PlottingStuff):

#     def __init__(self, image_or_expInstance):
#         if type(image_or_expInstance) == np.ndarray:
#             super().__init__(image_or_expInstance)
#         elif type(image_or_expInstance) == MolCoupNanoRodExp:
#             super().__init__(image_or_expInstance.trial_images)
#             self.instance_to_fit = image_or_expInstance

#     def plot_image_noised(self, image, PEAK=1):
#         noised_image_data = self.make_image_noisy(image, PEAK)
#         self.plot_raveled_image(noised_image_data)

#     def plot_image_noised_from_params(self, fit_params, PEAK=1):
#         image = self.raveled_model_of_params(fit_params)
#         self.plot_image_noised(image, PEAK)

#     def make_image_noisy(self, image, PEAK=1):
#         if image.ndim == 2:   ## set of raveled images
#             max_for_norm = image.max(axis=1)[:,None]
#         elif image.ndim == 1:
#             max_for_norm = image.max()
#         normed_image = image/max_for_norm
#         noised_normded_image_data = (
#             np.random.poisson(normed_image/255.0* PEAK) / (PEAK) *255
#             )
#         noised_image_data = noised_normded_image_data*max_for_norm
#         return noised_image_data

#     def make_noisy_image_attr(self, PEAK=1):
#         if hasattr(self, 'noised_image_data'):
#             return None
#         noised_image_attr = self.make_image_noisy(self.image_data, PEAK)
#         self.noised_image_data = noised_image_attr

#     def fit_model_to_noised_model_image_data(self, PEAK=5000):
#         if not hasattr(self, 'noised_image_data'):
#             self.make_noisy_image_attr(PEAK)
#         self.model_fit_results = self.fit_model_to_image_data(
#             images=self.noised_image_data)
#         return self.model_fit_results

# #     def plot_fit_results(self):
#     def plot_fit_results_as_quiver_map(
#         self, fitted_exp_instance, plot_limits=None):
#         '''...'''
#         if not hasattr(self, 'model_fit_results'):
#             self.fit_model_to_noised_model_image_data()
#         if plot_limits == None:
#             plot_limits = fitted_exp_instance.default_plot_limits
#         quiv_ax, = self.quiver_plot(
#             fitted_exp_instance.mol_locations[:,0],
#             fitted_exp_instance.mol_locations[:,1],
#             self.model_fit_results[:,2],
#             plot_limits,
#             true_mol_angle=fitted_exp_instance.mol_angles,
#             nanorod_angle = fitted_exp_instance.rod_angle,
#             title=r'Model Fit Pol. and Loc. w/noise',
#             )
#         self.scatter_centroids_wLine(
#             fitted_exp_instance.mol_locations[:,0],
#             fitted_exp_instance.mol_locations[:,1],
#             self.model_fit_results[:,:2].T,
#             quiv_ax,
#             )
#         return quiv_ax

#     def plot_contour_fit_over_data(self, image_idx):
#         ax = self.plot_raveled_image(self.noised_image_data[image_idx])
#         self.plot_image_from_params(self.model_fit_results[image_idx], ax)


# Testing fits

# In[89]:

def fixed_ori_mol_placement(
    x_min=0, x_max=350, y_min=0, y_max=350, mol_grid_pts_1D = 10, mol_angle=0):
    locations = diffi.observation_points(
        x_min, x_max, y_min, y_max, points=mol_grid_pts_1D
        )[0]
    locations = np.hstack((locations,np.zeros((locations.shape[0],1))))

    mol_linspace_pts = mol_grid_pts_1D
#     random_mol_angles= (np.random.random(mol_linspace_pts**2)*np.pi*2)
    return [locations, mol_angle]

def random_ori_mol_placement(
    x_min=0, x_max=350, y_min=0, y_max=350, mol_grid_pts_1D = 10):
    locations = diffi.observation_points(
        x_min, x_max, y_min, y_max, points=mol_grid_pts_1D
        )[0]
    locations = np.hstack((locations,np.zeros((locations.shape[0],1))))

    mol_linspace_pts = mol_grid_pts_1D
    random_mol_angles_0To360= (np.random.random(mol_linspace_pts**2)*np.pi*2)
    return [locations, random_mol_angles_0To360]

if __name__ == '__main__':
    '''This shit is all broken, or at least um_per_nmaintained'''

    print('Why hello there.')