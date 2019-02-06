''' '''
import numpy as np
from scipy import interpolate
import matplotlib as mpl

saved_mapping = np.loadtxt('obs_pol_vs_true_angle.txt')
true_ord_angles, obs_ord_angles =  saved_mapping.T

f = interpolate.interp1d(true_ord_angles,obs_ord_angles)
f_inv = interpolate.interp1d(
    obs_ord_angles[:251],
    true_ord_angles[:251],
    bounds_error=False,
    fill_value=(0,np.pi/2)
    )

def mb_p_fields(dipole_mag_array, dipole_coordinate_array):
    ''' As of 081418,fixing: currently only treats dipole at origin.'''    
    p = dipole_mag_array
#     print('Inside mb_p_fields, p= ',p)
    bfx = dipole_coordinate_array

    v_rel_obs_x_pts = (eye[1].ravel()[:,None] - bfx.T[0]).T
    v_rel_obs_y_pts = (eye[2].ravel()[:,None] - bfx.T[1]).T
    
    px_fields = np.asarray(afi.E_field(0, v_rel_obs_x_pts, v_rel_obs_y_pts, omega_drive*n_b/c))
    py_fields = np.asarray(afi.E_field(np.pi/2, v_rel_obs_x_pts, v_rel_obs_y_pts, omega_drive*n_b/c))
    pz_fields = np.zeros(py_fields.shape)
#     print('px_fields.shape=',px_fields.shape)
    ## returns [Ex, Ey, Ez] for dipoles oriented along cart units
    
    Ex = p[:,0,None]*px_fields[0] + p[:,1,None]*py_fields[0] + p[:,2,None]*pz_fields[0]
    Ey = p[:,0,None]*px_fields[1] + p[:,1,None]*py_fields[1] + p[:,2,None]*pz_fields[1]
    Ez = p[:,0,None]*px_fields[2] + p[:,1,None]*py_fields[2] + p[:,2,None]*pz_fields[2]

    return np.array([Ex,Ey,Ez])

def dipole_fields(locations, mol_angle=0, plas_angle=np.pi/2):
    d = locations*nm
    p0, p1 = cp.dipole_magnitudes(mol_angle, plas_angle, d_col=d, E_d_angle=None)
    mol_E = mb_p_fields(dipole_mag_array=p0, dipole_coordinate_array=d) 
    plas_E = mb_p_fields(dipole_mag_array=p1, dipole_coordinate_array=np.zeros(d.shape))

    # p0_unc, = cp.uncoupled_p0(mol_angle=0, d_col=d[0,None], E_d_angle=None)
    p0_unc, = cp.uncoupled_p0(mol_angle, E_d_angle=None)
#     p0_unc_E = mb_p_fields(dipole_mag_array=p0_unc[None,:], dipole_coordinate_array=np.zeros(d[0][None,:].shape))
    p0_unc_E = mb_p_fields(dipole_mag_array=p0_unc[None,:], dipole_coordinate_array=d) 
    
    return [mol_E, plas_E, p0_unc_E, p0, p1]

def powers_and_angels(E):
    drive_I = np.abs(parameters['general']['drive_amp'])**2.
    
    normed_Ix = np.abs(E[0])**2. / drive_I
    normed_Iy = np.abs(E[1])**2. / drive_I

    Px_per_drive_I = np.sum(normed_Ix,axis=-1) / sensor_size**2.
    Py_per_drive_I = np.sum(normed_Iy,axis=-1) / sensor_size**2.
    

    angles = np.arctan(Py_per_drive_I**0.5/Px_per_drive_I**0.5)
    return [angles, Px_per_drive_I, Py_per_drive_I]

def powers_and_angels_no_interf(E1,E2):
    drive_I = np.abs(parameters['general']['drive_amp'])**2.
    
    normed_Ix = (np.abs(E1[0])**2. + np.abs(E2[0])**2.) / drive_I
    normed_Iy = (np.abs(E1[1])**2. + np.abs(E2[1])**2.) / drive_I

    Px_per_drive_I = np.sum(normed_Ix,axis=-1) / sensor_size**2.
    Py_per_drive_I = np.sum(normed_Iy,axis=-1) / sensor_size**2.
    

    angles = np.arctan(Py_per_drive_I**0.5/Px_per_drive_I**0.5)
    return [angles, Px_per_drive_I, Py_per_drive_I]

def rotating_dipole_fields(locations):
    d = locations*nm
    p0, p1 = cp.dipole_magnitudes(mol_angle=np.linspace(0,2*np.pi,locations.shape[0]), 
                                  plas_angle=np.pi/2, 
                                  d_col=d, 
                                  E_d_angle=None)
    mol_E = mb_p_fields(dipole_mag_array=p0, dipole_coordinate_array=d) 
    plas_E = mb_p_fields(dipole_mag_array=p1, dipole_coordinate_array=np.zeros(d.shape))

    # p0_unc, = cp.uncoupled_p0(mol_angle=0, d_col=d[0,None], E_d_angle=None)
#     p0_unc, = cp.uncoupled_p0(mol_angle=np.linspace(0,2*np.pi,d.shape[0]), d_col=d, E_d_angle=None)
    p0_unc, = cp.uncoupled_p0(mol_angle=np.linspace(0,2*np.pi,d.shape[0]), E_d_angle=None)
#     p0_unc_E = mb_p_fields(dipole_mag_array=p0_unc[None,:], dipole_coordinate_array=np.zeros(d[0][None,:].shape))
    p0_unc_E = mb_p_fields(dipole_mag_array=p0_unc, dipole_coordinate_array=d) 
    
    return [mol_E, plas_E, p0_unc_E]


def quiver_plot(x_plot, y_plot, angles, plot_limits=[-25,550],
               title=r'Apparent pol. per mol. pos.'):

    fig, (ax0, ax_cbar) = plt.subplots(nrows=1,ncols=2, figsize=(3.25,3), dpi=300, 
                                       gridspec_kw = {'width_ratios':[6, 0.5]}
                                      )

    cmap = mpl.cm.nipy_spectral
    quiv = ax0.quiver(x_plot, y_plot, 
                      np.cos(angles),
                      np.sin(angles), 
                      angles,
                      cmap=cmap,
                      clim = [0, np.pi/2], 
                      width=0.01,
                      scale=15,
            #            scale_units='width',
                      pivot='mid'
                      )
#     qk = ax_cbar.quiverkey(quiv, 0.5, 0.105, 1, r'molecule orientation', labelpos='E',
#                        coordinates='figure')

    plas_mark = ax0.quiver(0, 0, 0,.3, color='grey',
    #                        marker='o',
    #                        s=2000
                           width=0.03,
                           scale=1,
                           units='width',
                           pivot='mid'
                          )
    ax0.axis('equal')
    ax0.set_xlim(plot_limits)
    ax0.set_ylim(plot_limits)
    ax0.set_title(title)
    ax0.set_xlabel('x [nm]')
    ax0.set_ylabel('y [nm]')

    norm = mpl.colors.Normalize(vmin=0, vmax=np.pi/2)

    cb1 = mpl.colorbar.ColorbarBase(ax_cbar, cmap=cmap,
                                    norm=norm,
                                    orientation='vertical')
    cb1.set_label(r'observed angle $\phi$')
    plas_dot = ax0.scatter(0,0,color='k',s=30)
    cb1.set_ticks([0, np.pi/8, np.pi/4, np.pi/8 * 3, np.pi/2])
    cb1.set_ticklabels([r'$0$', r'$\pi/8$',r'$\pi/4$',r'$3\pi/8$',r'$\pi/2$'])

#     fig.tight_layout()
    
    quiver_axis_handle = ax0
    return [quiver_axis_handle]

def plot_phase_quiver(x_plot, y_plot, angles, p0, p1, plot_limits=[-25,550],
                      title=r'Apparent pol. w/ dipole phase'):

    phase_diff = np.arccos(np.real(p0/p1)/np.abs(p0/p1))
    
    fig, (ax0, ax_cbar) = plt.subplots(nrows=1,ncols=2, figsize=(3.25,3), dpi=300, 
                                       gridspec_kw = {'width_ratios':[6, 0.5]}
                                      )
    scat_color = phase_diff[:,0]
    cmap = 'brg'

    ax0.scatter(x_plot, y_plot, c=scat_color, s=100,
               cmap=cmap,
               clim = [0, np.pi], 
    #            width=0.01,
    #            scale=10,
    #            scale_units='width',
    #            pivot='mid'
              )

    ax0.set_xlim(plot_limits)
    ax0.set_ylim(plot_limits)
    ax0.set_xlabel('x [nm]')
    ax0.set_ylabel('y [nm]')

    norm = mpl.colors.Normalize(vmin=0, vmax=np.pi)
    cb1 = mpl.colorbar.ColorbarBase(ax_cbar, cmap=cmap,
                                    norm=norm,
                                    orientation='vertical',
                                    ticks=[0,np.pi/4,np.pi/2,3*np.pi/4,np.pi],
                                    )
    # ax_cbar.set_yticks([0,np.pi/4,np.pi/2,3*np.pi/4,np.pi])
    # cb1.set_ticks([0,np.pi/4,np.pi/2,3*np.pi/4,np.pi])
    cb1.set_ticklabels([r'$0$',r'$\pi/4$',r'$\pi/2$',r'$3\pi/4$',r'$\pi$'])
    cb1.set_label(r'difference in dipole phase')
    ax0.set_title(title)


    plas_arrow = ax0.quiver(0, 0, 0,.3, color='grey',
    #                        marker='o',
    #                        s=2000
                           width=0.03,
                           scale=1,
                           units='width',
                           pivot='mid',
                           headlength=3,
                           headaxislength=3
                          )
    plas_dot = ax0.scatter(0,0,color='k',s=200)

    quiv = ax0.quiver(x_plot, y_plot, np.cos(angles),np.sin(angles), 
                    color='white',
    #                    cmap='inferno',
    #                    clim = [0, np.pi/2], 
                       width=0.01,
                       scale=12,
            #            scale_units='width',
                       pivot='mid',
    #                   linewidth=100.,
                      headaxislength=0.0,
                      headlength=0.0
                      )
    # qk = ax_cbar.quiverkey(quiv, 0.5, 0.105, 1, r'molecule orientation', labelpos='E',
    #                    coordinates='figure', color='k' )

    # fig.tight_layout()

    # plt.colorbar(ticks=[0,np.pi/4,np.pi/2,3*np.pi/4,np.pi])
