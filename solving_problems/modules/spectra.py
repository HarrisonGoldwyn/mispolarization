''' for python 3 '''
import numpy as np 
import yaml

opened_constant_file = open('../physical_constants.yaml','r')
constants = yaml.load(opened_constant_file)
e = constants['physical_constants']['e']
c = constants['physical_constants']['c']  # charge of electron in statcoloumbs
hbar =constants['physical_constants']['hbar']
nm = constants['physical_constants']['nm']
mm = 1e6*nm
n_a = constants['physical_constants']['nA']

path_to_yaml = '../curly_param.yaml'
print('reading parameters from {}'.format(path_to_yaml))
opened_param_file = open(path_to_yaml,'r')
parameters = yaml.load(opened_param_file)
print(path_to_yaml)
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

import coupled_dipoles as cp




import sys
import os
sys.path.append('../Optics')
import diffraction_int as diffi
import fibonacci as fib
sys.path.append('../field_functions')
import far_fields as fi





############# Defining varablies for spherical integration 

magnification = 1
numerical_aperture = 0.999
max_theta = np.arcsin(numerical_aperture) # defines physical aperture size

## numerical parameters for calculation of scattered field
lens_points = 300
obj_f = 2*mm
tube_f = magnification * obj_f

sphere_points = fib.fib_alg_k_filter(
    num_points=lens_points, 
    max_ang=max_theta
    )

cart_points_on_sph = fib.sphere_to_cart(
    sphere_points[:,0],
    sphere_points[:,1],
    obj_f*np.ones(np.shape(sphere_points[:,0]))
    )

r_hat = cart_points_on_sph / obj_f

refracted_coords = diffi.refract_sph_coords(
    sphere_points, obj_f, tube_f)




###### Need to figure out the data structure so I can get spectra as functino of dipole orientation and location. 
###### dioples come out in shape = ( number of dipoles, 3 {cartesient coords} )


min_energy = parameters['spectra']['min_energy']
max_energy = parameters['spectra']['max_energy']
freq_step_size = parameters['spectra']['energy_step_size']
energy_axis = np.arange(min_energy, max_energy, freq_step_size)
# energies = np.arange(2.6, 3, 1)


#######################################################################
## Begin loop through seperations
#######################################################################
def single_dip_scat_spec(dip_number, angle=0, energies=energy_axis):
    if dip_number == 0:
        dip = cp.uncoupled_p0
    elif dip_number == 1:
        dip = cp.uncoupled_p1
    else: raise ValueError(
        "first argument of 'single_dip_scat_spec' requires value of",
        " '0' or '1'"
        )

    number_of_energies = energies.shape[0]

    scattering_spectrum = np.zeros((number_of_energies))

    for j in range(energies.shape[0]):
        ## calculate complec magnitufe vectors 
        p = dip(angle, E_d_angle=None,
            drive_hbar_w=energies[j])

        mag_p = cp.vec_mag(p)
        p_hat = p/cp.vec_mag(p) 
        phi = np.angle(mag_p)

        # print('beggining location loop at {} eV'.format(energies[j]))

        # for i in range(number_of_locations):
        x = np.array([[0,0,0]])
           

        #######################################################################
            ## Compute diffracted fields on GRS and energy flux on CCD
        #######################################################################
            # print(
            #     'Calculating fields at sep = {}'.format(locations.ravel()[i]))
        sph_coords_p = cart_points_on_sph - x
            
        E_p_sph = fi.E_dipole_complex(
            w=energies[j]/hbar, 
            dipole_magnitude=mag_p, 
            dipole_phase=phi, 
            dipole_ori_unit=p_hat, 
            r=sph_coords_p
            )
        # E_p1_sph = fi.E_dipole_complex(
        #     w=energies[j]/hbar, 
        #     dipole_magnitude=mag_p1, 
        #     dipole_phase=phi1, 
        #     dipole_ori_unit=p1_hat, 
        #     r=sph_coords_p1
        #     )

        B_p_sph = fi.B_dipole_complex(
            w=energies[j]/hbar, 
            dipole_magnitude=mag_p, 
            dipole_phase=phi, 
            dipole_ori_unit=p_hat, 
            r=sph_coords_p
            )
            # B_p1_sph = fi.B_dipole_complex(
            #     w=energies[j]/hbar, 
            #     dipole_magnitude=mag_p1, 
            #     dipole_phase=phi1, 
            #     dipole_ori_unit=p1_hat, 
            #     r=sph_coords_p1
            #     )

        ## Calculate total Poynting vector on large sphere
        complex_s_tot = fi.S_complex(
                E_p_sph,B_p_sph
                )
        avrg_S = np.real(complex_s_tot)

        #### Perform integral of (n_hat * S dA)
        n_dot_S = np.sum(avrg_S * r_hat, axis=-1)
        number_of_lens_points = cart_points_on_sph.shape[0]
        lens_surface_area = (
            2*np.pi  
            * 
            obj_f**2. 
            )
        area_per_point = lens_surface_area/number_of_lens_points
        total_power_radiated = np.sum(n_dot_S)/area_per_point

        ## Divide by incident field intensity to normamalize, units of per area
        inc_field_intens = ( c /(8*np.pi) ) * np.abs(
            parameters['general']['drive_amp'])**2.
        cross_section = total_power_radiated / inc_field_intens

        # print(total_power_radiated)
        scattering_spectrum[j] = cross_section
    return [energies,scattering_spectrum] 

def generalized_single_dip_scat_spec(
        polarizability, energies=energy_axis, rel_drive_angle=0
        ):
    ''' first argument of 'polarizability' bust be hbar*W '''

    number_of_energies = energies.shape[0]
    scattering_spectrum = np.zeros((number_of_energies))

    for j in range(energies.shape[0]):
        ## calculate complec magnitufe vectors 
        alpha_dyad = polarizability(energies[j])
        E_drive = cp.rotation_by(rel_drive_angle) @ np.array([1,0,0])*ficticious_field_amp
        p = alpha_dyad @ E_drive
        p = p[None,:]
        # print(p)
        mag_p = cp.vec_mag(p)
        p_hat = p/cp.vec_mag(p) 
        phi = np.angle(mag_p)

        # print('beggining location loop at {} eV'.format(energies[j]))

        # for i in range(number_of_locations):
        x = np.array([[0,0,0]])
        

        #######################################################################
            ## Compute diffracted fields on GRS and energy flux on CCD
        #######################################################################
            # print(
            #     'Calculating fields at sep = {}'.format(locations.ravel()[i]))
        sph_coords_p = cart_points_on_sph - x
            
        E_p_sph = fi.E_dipole_complex(
            w=energies[j]/hbar, 
            dipole_magnitude=mag_p, 
            dipole_phase=phi, 
            dipole_ori_unit=p_hat, 
            r=sph_coords_p
            )
        # E_p1_sph = fi.E_dipole_complex(
        #     w=energies[j]/hbar, 
        #     dipole_magnitude=mag_p1, 
        #     dipole_phase=phi1, 
        #     dipole_ori_unit=p1_hat, 
        #     r=sph_coords_p1
        #     )

        B_p_sph = fi.B_dipole_complex(
            w=energies[j]/hbar, 
            dipole_magnitude=mag_p, 
            dipole_phase=phi, 
            dipole_ori_unit=p_hat, 
            r=sph_coords_p
            )
            # B_p1_sph = fi.B_dipole_complex(
            #     w=energies[j]/hbar, 
            #     dipole_magnitude=mag_p1, 
            #     dipole_phase=phi1, 
            #     dipole_ori_unit=p1_hat, 
            #     r=sph_coords_p1
            #     )

        ## Calculate total Poynting vector on large sphere
        complex_s_tot = fi.S_complex(
                E_p_sph,B_p_sph
                )
        avrg_S = np.real(complex_s_tot)

        #### Perform integral of (n_hat * S dA)
        n_dot_S = np.sum(avrg_S * r_hat, axis=-1)
        number_of_lens_points = cart_points_on_sph.shape[0]
        lens_surface_area = (
            2*np.pi  
            * 
            obj_f**2. 
            )
        area_per_point = lens_surface_area/number_of_lens_points
        total_power_radiated = np.sum(n_dot_S)/area_per_point

        ## Divide by incident field intensity to normamalize, units of per area
        inc_field_intens = ( c /(8*np.pi) ) * np.abs(
            parameters['general']['drive_amp'])**2.
        cross_section = total_power_radiated / inc_field_intens

        # print(total_power_radiated)
        scattering_spectrum[j] = cross_section
    return [energies,scattering_spectrum] 






    spectrum = dip(0, E_d_angle=None, 
        drive_hbar_w=parameters['general']['drive_energy'])

def calculate_scatt_spectra(mol_angle, plas_angle, locations, energies=energy_axis):
    
    number_of_locations = locations.shape[0]
    number_of_energies = energies.shape[0]
    ## Initialize output spectra
    totral_scattering_spectra_per_seperation = np.zeros(
        (number_of_locations, number_of_energies))

    for j in range(energies.shape[0]):
        ## calculate complec magnitufe vectors 
        p0, p1 = cp.dipole_magnitudes(
            mol_angle, plas_angle, d_col=locations, E_d_angle=None,
            drive_hbar_w=energies[j])

        mag_p0 = cp.vec_mag(p0)
        p0_hat = p0/cp.vec_mag(p0) 
        phi0 = np.angle(mag_p0)

        mag_p1 = cp.vec_mag(p1) 
        p1_hat = p1/cp.vec_mag(p1) 
        phi1 = np.angle(mag_p1)

        # print('beggining location loop at {} eV'.format(energies[j]))

        for i in range(number_of_locations):
            x0 = locations[i]
            x1 = np.array([[0,0,0]])
           

        #######################################################################
            ## Compute diffracted fields on GRS and energy flux on CCD
        #######################################################################
            # print(
            #     'Calculating fields at sep = {}'.format(locations.ravel()[i]))
            sph_coords_p0 = cart_points_on_sph - x0
            sph_coords_p1 = cart_points_on_sph - x1

            # di_E_field_on_sph = (
            #     fi.E_dipole_complex(
            #         w=energies[j]/hbar, 
            #         dipole_magnitude=mag_p0[i], 
            #         dipole_phase=phase1[i], 
            #         dipole_ori_unit=dipole_orientation, 
            #         r=diff_r_d1
            #         )
            #     +
            #     fi.E_dipole_complex(
            #         w=energies[j]/hbar, 
            #         dipole_magnitude=q2[i], 
            #         dipole_phase=phase2[i], 
            #         dipole_ori_unit=dipole_orientation, 
            #         r=diff_r_d2
            #         )
            #     )
            E_p0_sph = fi.E_dipole_complex(
                w=energies[j]/hbar, 
                dipole_magnitude=mag_p0[i], 
                dipole_phase=phi0[i], 
                dipole_ori_unit=p0_hat[i], 
                r=sph_coords_p0
                )
            B_p0_sph = fi.B_dipole_complex(
                w=energies[j]/hbar, 
                dipole_magnitude=mag_p0[i], 
                dipole_phase=phi0[i], 
                dipole_ori_unit=p0_hat[i], 
                r=sph_coords_p0
                )


            E_p1_sph = fi.E_dipole_complex(
                w=energies[j]/hbar, 
                dipole_magnitude=mag_p1[i], 
                dipole_phase=phi1[i], 
                dipole_ori_unit=p1_hat[i], 
                r=sph_coords_p1
                )
            B_p1_sph = fi.B_dipole_complex(
                w=energies[j]/hbar, 
                dipole_magnitude=mag_p1[i], 
                dipole_phase=phi1[i], 
                dipole_ori_unit=p1_hat[i], 
                r=sph_coords_p1
                )

            ## Calculate total Poynting vector on large sphere
            complex_s_tot = fi.S_complex(
                    E_p0_sph+E_p1_sph,
                    B_p0_sph+B_p1_sph
                    )
            avrg_S = np.real(complex_s_tot)

            #### Perform integral of (n_hat * S dA)
            n_dot_S = np.sum(avrg_S * r_hat, axis=-1)
            number_of_lens_points = cart_points_on_sph.shape[0]
            lens_surface_area = (
                2*np.pi  
                * 
                obj_f**2. 
                )
            area_per_point = lens_surface_area/number_of_lens_points
            total_power_radiated = np.sum(n_dot_S)/area_per_point

            ## Divide by incident field intensity to normamalize, units of per area
            inc_field_intens = ( c /(8*np.pi) ) * np.abs(
                parameters['general']['drive_amp'])**2.
            cross_section = total_power_radiated / inc_field_intens


            # print(total_power_radiated)
            totral_scattering_spectra_per_seperation[i, j] = cross_section
    return [energies,totral_scattering_spectra_per_seperation] 




