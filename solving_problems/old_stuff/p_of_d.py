''' 
Would like this file just to define functions for conputation so that no
modificantions need to be made here, but only in script files which call 
functions from here. 
'''
from __future__ import print_function
from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import yaml
import sys
import os
## import diffrantion integral solver from Optics folder
work_dir = os.getcwd()
date_dir = os.path.split(work_dir)[0]
optics_folder = os.path.join(date_dir, 'Optics')
sys.path.append(optics_folder)
import diffraction_int as diffi
import fibonacci as fib

## Import field functions
field_module_folder = os.path.join(date_dir, 'field_functions')             
sys.path.append(field_module_folder)

## grab yaml file path
full_path_to_yaml = os.path.join(work_dir,'param.yaml')
print('reading parameters from {}'.format(full_path_to_yaml))
opened_param_file = open(full_path_to_yaml,'r')
parameters = yaml.load(opened_param_file)

## Load physical constants from yaml file
full_path_to_constant_yaml = os.path.join(date_dir,'physical_constants.yaml')
opened_constant_file = open(full_path_to_constant_yaml,'r')
constants = yaml.load(opened_constant_file)
e = constants['physical_constants']['e']
c = constants['physical_constants']['c']  # charge of electron in statcoloumbs
hbar =constants['physical_constants']['hbar']
nm = constants['physical_constants']['nm']
n_a = constants['physical_constants']['nA']   # Avogadro's number
# Z_o = 376.7303 # impedence of free space in ohms (SI)


## System background
n_b = parameters['general']['background_ref_index']
eps_b = n_b**2.


## Plasmon Drude properties
hbar_w_p = parameters['plasmon']['plasma_energy']
drude_damping_energy = parameters['plasmon']['drude_damping_energy']
eps_inf = parameters['plasmon']['eps_inf']

## Other plasmon properties
a = parameters['plasmon']['radius']


## Driving force
ficticious_field_amp = parameters['general']['drive_amp'] 
## this is a problem for the ss
# drive_hbar_w = parameters['general']['drive_energy']


## Print whether parameter file say include quadrapole?
if parameters['plasmon']['include_quadrupole'] == False:
    print( 'Quadrapole plasmon NOT included.' )
elif parameters['plasmon']['include_quadrupole'] == True:
    print( 'Quadrapole plasmon included.' )
else: print( 'So what are we doing here?' )


############################################################################
############################################################################
##### Functions ... ######
############################################################################
############################################################################


############################################################################
### Functions for converting parameters to computationally useful values ###
############################################################################

def unit_vector(row_vecs): 
    '''Take argument as matrix of row vectors and divide by its linalg.norm.  
    Probably should make sure the array is 2D, but I wont for now. 
    '''
    if row_vecs.ndim == 2:
        vals = 1/np.linalg.norm(row_vecs, axis=(1))  # breaks if mag == 0, ok?
        lenth = vals.size
        arranged = vals.reshape([lenth, 1])
        return np.multiply(arranged, row_vecs)
    else: raise ValueError('Can\'t \'hat()\' array if ndim != 2')

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


############################################################################
### Functions for converting parameters to computationally useful values ###
############################################################################

def eV_to_Hz(energy):
    return energy/hbar

def seperations(min_sep, max_sep, steps):
    return np.linspace(min_sep, max_sep, steps)

def seperation_columb(min_sep, max_sep, steps):
    return np.linspace(min_sep, max_sep, steps).reshape(steps,1)


############################################################################
### Functions for generating coupling constants for equations of motion  ###
############################################################################
def pp_coupling_constant(drive_hbar_w,
    d1_pos_unit, d1_or, d2_pos_unit, d2_or, seperation_col,
    ):
    ''' returns columb corresponding to seperations '''
    d1_pos_unit = np.array(d1_pos_unit)
    d2_pos_unit = np.array(d2_pos_unit)
    
    d1_pos = d1_pos_unit * (seperation_col/2.)
    d2_pos = d2_pos_unit * (seperation_col/2.)
    the_n = unit_vector(d1_pos - d2_pos)  # overkill

    d1_or = np.array([d1_or])
    d2_or = np.array([d2_or])
    # r = row_mag(d1_pos - d2_pos)
    r = seperation_col
    # print('r = ', r)
    w = eV_to_Hz(drive_hbar_w)
    k = w * n_b / c
    complex_phase_factor = np.exp(1j*k*r)
    nf_if_vecs = (
        3.*np.multiply(
            np.dot(the_n, d1_or.T), 
            np.dot(the_n, d2_or.T)
            )
        - np.dot(d1_or, d2_or.T)
        )
    # print('nf_if_vecs = ',nf_if_vecs)
    ff_vecs = (
        np.multiply(
            np.dot(the_n, d1_or.T), 
            np.dot(the_n, d2_or.T)
            )
        - np.dot(d1_or, d2_or.T)
        )
    # print('ff_vecs = ',ff_vecs)
    g_dip_dip = (
        e**2. 
        * 
        # np.real(
        (    
            complex_phase_factor*(
                (nf_if_vecs)/r**3.
                - 
                (1j * k * nf_if_vecs)/r**2.
                - 
                (k**2. * ff_vecs)/r
                )
            ) 
        )
    # radiation_backforce_coupling = 1j * 2. * w**3. * e**2. / (3 * c**3.)
    radiation_backforce_coupling = 0
    
    g = g_dip_dip + radiation_backforce_coupling
    
    # print('g = ',g)
    # g = np.zeros(seperation_col.shape)
    return g

def qxx_px_ori_specific_coupling_constant(seperation_col):
    ''' only xx to x '''
    r = seperation_col
    g = 3. * e**2. * a / r**4.

    return g 


#############################################################################
######## Determine masses based on quasistatic derivation #####
#############################################################################

def plasmon_mass_QS(
    res_freq, 
    radius=a, 
    eps_inf= eps_inf):
    ''' Old way of calculating the mass, not sure where the 4 pi comes in
    '''
    m = e**2. * (eps_inf + 2.)/(3. * res_freq**2. * radius**3.)
    return m

def plasmon_lth_resonance_QS(l):
    w_l = eV_to_Hz(hbar_w_p)*np.sqrt(
        l
        /
        (l*(eps_inf + eps_b) + eps_b)
        )
    
    in_eV = hbar * w_l
    shift_to_ret = in_eV - 0.005
    w_l = shift_to_ret / hbar

    # print('w at l={}, is {} eV'.format(l, w_l*hbar))

    return w_l 

def plasmon_lth_mass(l):
    w_l = plasmon_lth_resonance_QS(l)

    jakes_prefactor = 4 * np.pi * e**2. * eps_b / a**3. 
    arpc_prefactor = 3 * e**2. / a**3.
    elliots_prefactor =  e**2. / a**3.
    harrisons_prefactor = e**2. / (a**3. * eps_b)
   
    jakes_ratio = (
        ( l * ( eps_inf + eps_b ) + eps_b )
        /
        ( w_l * ( 2*l + 1 ) )**2.
        )
    arpc_ratio = (
        ( l * ( eps_inf + eps_b ) + eps_b )
        /
        ( w_l**2. * ( 2*l + 1 ) )
        )

    m = (
        # jakes_prefactor / ( 4 * np.pi )
        # *
        # arpc_prefactor
        # *
        elliots_prefactor
        *
        arpc_ratio
        )
    return m

#############################################################################
#### Determine fluorophore mass by relating extinction coefficient to 
#### SHO absorption, maybe should incorperate scattering...  but it doesnt
#### seem to matter if the fluo mass changes slightly.                  #####
#############################################################################

def fluorophore_mass(ext_coef, gamma):
    '''Derived at ressonance'''
    m = 4 * np.pi * e**2 * n_a  / (
            ext_coef * np.log(10) * c * n_b * gamma
            )
    return m

# def fluo_mass_w_freq_dep(ext_coef, gamma, drive_omeg, res_omega):
#     '''I think this is wrong'''
#     pref = 2 * np.pi * e**2 * n_a / (
#         3 * ext_coef * np.log(10) * c * n_b
#         )
#     frac = (gamma*drive_omeg**2)/(
#         (res_omega**2-drive_omeg**2)**2+(gamma**2*drive_omeg**2)
#         )
#     m = pref*frac
#     return m


#############################################################################
#############################################################################
############# Functions for computation of fourier amplitude coefficients ###
#############################################################################
#############################################################################

def ressonant_denom(drive_hbar_w, res_freq, non_rad_damping, alpha = 0):
    w = eV_to_Hz(drive_hbar_w)
    # print(res_freq)
    denom = (
        res_freq**2. 
        - w**2. 
        - 1j*( non_rad_damping + alpha * w**2. ) * w
        )

    return denom


#############################################################################
#### Generate complex Fourier amplitudes for the oscillators 
#### and then observables
#############################################################################


def fluo_amp_Four(fluo_mass, fluo_res,
    fluo_damp, drive_hbar_w, drive_force, sum_gSqrd_on_mD):
    ''' 
    Fourier space complex amplitude for calculation of observable
    '''
    complex_amp = (
        drive_force
        /
        (
            fluo_mass
            *
            ressonant_denom(drive_hbar_w, fluo_res, fluo_damp)
            - 
            sum_gSqrd_on_mD
            )
        )
    # print('complex_amp, ', complex_amp, '\n',
    #     'drive_force, ', drive_force, '\n'
    #     'fluo_mass, ', fluo_mass, '\n'
    #     'ressonant_denom(drive_hbar_w, fluo_res, fluo_damp), ', ressonant_denom(drive_hbar_w, fluo_res, fluo_damp), '\n',
    #     'sum_gSqrd_on_mD', sum_gSqrd_on_mD
    #      )
    return complex_amp

def fluo_osci_mag_phase(fluo_mass, fluo_res,
    fluo_damp, drive_hbar_w, drive_force, sum_gSqrd_on_mD):
    ''' Computes the harmonic oscillation amplitude of the oscillator
    described by 'this' variables coupled to the 'other' oscillator.
    '''
    complex_amp = fluo_amp_Four(fluo_mass, fluo_res,
        fluo_damp, drive_hbar_w, drive_force, sum_gSqrd_on_mD)
    
    real_dipole_amplitude = -e*(
        complex_amp.real**2. 
        + 
        complex_amp.imag**2.
        )**0.5

    phase = np.arctan2(complex_amp.imag,complex_amp.real)
    # phase = np.arccos(
    #     (complex_amp.real / (real_dipole_amplitude/e))
    #         )

    return [real_dipole_amplitude, phase, complex_amp]


## plasmon stuff

def plas_amp_Four(g_on_mD, fluo_amp_Four):
    complex_amp = g_on_mD * fluo_amp_Four
    return complex_amp


def static_drude_sphere_polarizability(
        plasmon_mass,
        plas_l1_ress_hbar_omega 
        ):
    # alpha_stat = (
    #     (e**2. * (eps_inf - eps_b))
    #     /
    #     (plasmon_mass * 3. * eps_b * eV_to_Hz(plas_l1_ress_hbar_omega)**2.)
    #     )

    alpha_stat = a**3. * (
        (eps_inf - eps_b)
        /
        (eps_inf + 2*eps_b)
        )

    return alpha_stat

def static_plas_amp(alpha_stat, g_l1, fluo_four_amp):
    q_plas_stat = alpha_stat * g_l1 * fluo_four_amp / e**2.

    return q_plas_stat


# def plas_osci_mag_phase(g_on_mD, fluo_mass, fluo_res,
def plas_osci_mag_phase(g_l1, plas_l1_mass, res_Denom, fluo_mass, fluo_res,
    fluo_damp, drive_hbar_w, drive_force, sum_gSqrd_on_mD,
    plas_l1_ress_hbar_omega,
    want_static_piece):
    ''' 
    Fourier space complex amplitude for calculation of observable
    '''

    complex_amp = plas_amp_Four(
        g_on_mD=(g_l1 / ( plas_l1_mass * res_Denom )), 
        fluo_amp_Four=fluo_amp_Four(
            fluo_mass= fluo_mass, 
            fluo_res= fluo_res,
            fluo_damp= fluo_damp, 
            drive_hbar_w= drive_hbar_w, 
            drive_force= drive_force, 
            sum_gSqrd_on_mD= sum_gSqrd_on_mD
            )
        )
        
    if want_static_piece == 'True':
        complex_amp += static_plas_amp(
            alpha_stat = static_drude_sphere_polarizability(
                plasmon_mass=plas_l1_mass,
                plas_l1_ress_hbar_omega=plas_l1_ress_hbar_omega
                ),
            g_l1=g_l1, 
            fluo_four_amp=fluo_amp_Four(
                fluo_mass= fluo_mass, 
                fluo_res= fluo_res,
                fluo_damp= fluo_damp, 
                drive_hbar_w= drive_hbar_w, 
                drive_force= drive_force, 
                sum_gSqrd_on_mD= sum_gSqrd_on_mD
                ),
            )

    real_dipole_amplitude = -e*(
        complex_amp.real**2. 
        + 
        complex_amp.imag**2.
        )**0.5

    phase = np.arctan2(complex_amp.imag, complex_amp.real)
    # phase = np.arccos(
    #     (complex_amp.real / (real_dipole_amplitude/-e))
    #         )
    

    return [real_dipole_amplitude, phase, complex_amp]




#############################################################################
#############################################################################
#### Here, we insert parameters to specifications in the parameter file
#############################################################################
#############################################################################

#############################################################################
#### define dipole and quadrapole parameters by their spectrally fit values
#### or Drude values according to specification in parameter file
#############################################################################

if parameters['plasmon']['incorperate_parameters_from_fit'] == True:
    plas_l1_ress_freq= eV_to_Hz(parameters['plasmon']['fit_hbar_w0'])
    plas_l1_mass= parameters['plasmon']['fit_mass']  
        ## with scaling to simplify entry
    plas_l1_non_rad_damp_freq = eV_to_Hz(
        parameters['plasmon']['fit_hbar_gamma'])
    

elif parameters['plasmon']['incorperate_parameters_from_fit'] == False:
    ## Drude values
    plas_l1_ress_freq= plasmon_lth_resonance_QS(1)
    plas_l1_mass= plasmon_lth_mass(1)  ## parameters default
    plas_l1_non_rad_damp_freq = eV_to_Hz(
        parameters['plasmon']['drude_damping_energy']
        )

print('plasmon_l1_mass =', plas_l1_mass)

## not fitting quadrapole yet, so ill just define it's parameters here
plas_l2_mass= parameters['plasmon']['fit_l2_mass'] 
# plas_l2_ress_freq = plasmon_lth_resonance_QS(2)
plas_l2_ress_freq = eV_to_Hz(
    parameters['plasmon']['fit_hbar_w2']
    )
plas_l2_non_rad_damp_freq = parameters['plasmon']['fit_l2_hbar_gamma']

print('plasmon_l2_mass =', plas_l2_mass)

#############################################################################
#### ressonant denominators, plasmon adjusted for radiation backforce
#############################################################################

def D_plas_l1(drive_hbar_w):
    ressonant_denominator = ressonant_denom(
        # res_freq=eV_to_Hz(parameters['plasmon']['p_res_energy']), 
        drive_hbar_w = drive_hbar_w,
        res_freq= plas_l1_ress_freq,
        non_rad_damping= plas_l1_non_rad_damp_freq,
        alpha = ( (2./3.) * e**2. / ( plas_l1_mass * c**3. ) ) 
        )
    return ressonant_denominator

def D_plas_l2(drive_hbar_w):
    ressonant_denominator = ressonant_denom(
        drive_hbar_w = drive_hbar_w,
        res_freq=plas_l2_ress_freq, 
        non_rad_damping=plas_l2_non_rad_damp_freq
        )
    return ressonant_denominator

#############################################################################
#### Asign fluorophore values
#############################################################################

if parameters['fluorophore']['make_sphere'] == False:
    fluo_mass = fluorophore_mass(
        ext_coef=parameters['fluorophore']['extinction_coeff'],
        gamma = eV_to_Hz(parameters['fluorophore']['mass_gamma']),
        # gamma= 1/parameters['fluorophore']['lifetime']
        )
    fluo_non_rad_damp_ene = parameters['fluorophore']['test_gamma']
    fluo_ress_ene = parameters['fluorophore']['res_energy']

elif parameters['fluorophore']['make_sphere'] == True:
    fluo_mass = plas_l1_mass
    fluo_non_rad_damp_ene = plas_l1_non_rad_damp_freq
    fluo_ress_ene = hbar * plas_l1_ress_freq 

# fluo_mass = fluo_mass_w_freq_dep(
#     ext_coef= parameters['fluorophore']['extinction_coeff'],
#     gamma = 1/parameters['fluorophore']['lifetime'],
#     drive_hbar_wa=eV_to_Hz(drive_hbar_w),
#     res_omega=eV_to_Hz(fluo_ress_ene),
    # )
# fluo_mass=plas_l1_mass*1e2
print('fluo mass =', fluo_mass)


#############################################################################
#############################################################################
#### Reduce functions that generate solution to only orientation dependence
#############################################################################
#############################################################################

#############################################################################
#### Dipole - dipole coupling
#############################################################################

def g_l1(mol_orientation, plas_orientation, drive_hbar_w, d):
    gp = pp_coupling_constant(
        drive_hbar_w,
        d1_pos_unit=parameters['plasmon']['pos_unit_vec'],
        d1_or=plas_orientation,
        d2_pos_unit=parameters['fluorophore']['pos_unit_vec'],
        d2_or=mol_orientation,
        seperation_col=d
        )
    # print('gp = ',gp)
    return gp


#############################################################################
#### Dipole - quadrapole coupliong
#############################################################################

def g_l2(dip_orientation,d):
    gq =  qxx_px_ori_specific_coupling_constant(
        # mol_orientation=dip_orientation,
        seperation_col=d
        )
    return gq


#############################################################################
# plasmon dipole term
#############################################################################
 
def sum_gSqrd_on_mD_terms(mol_orientation, plas_orientation, drive_hbar_w, d):
    gs_on_mds = (
        g_l1(mol_orientation, plas_orientation, drive_hbar_w,d)**2.
        /
        (plas_l1_mass*D_plas_l1(drive_hbar_w))
        )
        ## plasmon dipole term
    # print('sum_gSqrd_on_mD_terms without quad = ', gs_on_mds)
    if parameters['plasmon']['include_quadrupole']==True:
        gs_on_mds += (
        ## plasmon quadrapole term
            g_l2(mol_orientation,d)**2.
            /
            (plas_l2_mass*D_plas_l2(drive_hbar_w))
            )
        # print('sum_gSqrd_on_mD_terms WITH quad = ', gs_on_mds)
    # print(gs_on_mds)
    return gs_on_mds


def plas_l1_osc(mol_orientation, plas_orientation, drive_hbar_w, d): 
    '''returns amp and phase vectors
    '''
    plas_l1_amp, plas_l1_phase, c_amp = plas_osci_mag_phase(
        g_l1=g_l1(mol_orientation, plas_orientation, drive_hbar_w, d),
        plas_l1_mass=plas_l1_mass, 
        res_Denom=D_plas_l1(drive_hbar_w),
        fluo_mass=fluo_mass, 
        fluo_res=eV_to_Hz(fluo_ress_ene),
        fluo_damp=eV_to_Hz(fluo_non_rad_damp_ene), 
        drive_hbar_w=drive_hbar_w, 
        drive_force= -e * ficticious_field_amp, 
        sum_gSqrd_on_mD=sum_gSqrd_on_mD_terms(mol_orientation, plas_orientation, drive_hbar_w, d),
        plas_l1_ress_hbar_omega= plas_l1_ress_freq * hbar,
        want_static_piece='True'
        )
    # print('plasmon p amp = ',plas_l1_amp, '\n'
    #         'for or = ',orientation, '\n',
    #         'with quadrapole = ', parameters['plasmon']['include_quadrupole'])
    return [plas_l1_amp, plas_l1_phase, c_amp]

if parameters['plasmon']['include_quadrupole']==True:
    def plas_l2_osc(orientation, drive_hbar_w, d): 
        '''returns amp and phase vectors
        '''
        plas_l2_amp, plas_l2_phase, c_amp = plas_osci_mag_phase(
            g_l1=g_l2(drive_hbar_w, d),
            plas_l1_mass=plas_l2_mass, 
            res_Denom=D_plas_l2(drive_hbar_w),
            fluo_mass=fluo_mass, 
            fluo_res=eV_to_Hz(fluo_ress_ene),
            fluo_damp=eV_to_Hz(fluo_non_rad_damp_ene), 
            drive_hbar_w=drive_hbar_w, 
            drive_force= a * -e * ficticious_field_amp, 
            sum_gSqrd_on_mD=sum_gSqrd_on_mD_terms(orientation, drive_hbar_w, d),
            plas_l1_ress_hbar_omega= plas_l2_ress_freq * hbar,
            want_static_piece= 'False'
            )

        return [plas_l2_amp, plas_l2_phase, c_amp]

def fluo_osc(mol_orientation, plas_orientation, drive_hbar_w, d):
    '''
    '''
    fluo_amp, fluo_phase, c_amp = fluo_osci_mag_phase(
        fluo_mass=fluo_mass, 
        fluo_res=eV_to_Hz(fluo_ress_ene),
        fluo_damp=eV_to_Hz(fluo_non_rad_damp_ene), 
        drive_hbar_w=drive_hbar_w, 
        drive_force= -e * ficticious_field_amp, 
        sum_gSqrd_on_mD=sum_gSqrd_on_mD_terms(mol_orientation, plas_orientation, drive_hbar_w, d)
        )
    return [fluo_amp, fluo_phase, c_amp]



# if __name__ == "__main__":
#     print(
#         'plasmon mass = ', plas_l1_mass,'\n',
#         'fluo mass =', fluo_mass, '\n',
#         # 'plas omeg = ', eV_to_Hz(parameters['plasmon']['p_res_energy']),'\n',
#         'plas omeg = ', plas_l1_ress_freq,'\n',
#         'fluorophore mass = ', fluo_mass,'\n',
#         'fluo omeg = ', eV_to_Hz(fluo_ress_ene), '\n'
#         'drive omega = ', eV_to_Hz(drive_hbar_w))
