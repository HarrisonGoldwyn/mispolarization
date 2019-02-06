def make_fields(
    molecule_locations_in_plane=(np.array([[1,0]])*np.ones((3,2))*150*nm),
    molecule_angles=np.linspace(0,2*np.pi,3), ## default
    plas_orientation=[0,1,0] #pi
    ):
    

    if molecule_locations_in_plane.shape[0] != molecule_angles.shape[0]:
        print('first two inputs must have same length: shape[0] = number of molecules')  
        return None
    
    plas_orientation = np.array(plas_orientations).shape
    if plas_orientations != (3,):
        print('not ready for multiple plasmon orientatnions')  
        return None
    # define angles
    # number_of_molecules = 5
    # ordered_angles = np.linspace(0,2*np.pi,number_of_molecules)

    ## plas-mol separation
    d = dnm*nm
    ## make seperation vectors
    p0_hat = make_in_plane_cartesien_angles_array(molecule_angles)

    bf_x0 = molecule_locations_in_plane
    # np.array([[1,0,0]])*np.ones((number_of_molecules,2))*d
    bf_x1 = np.zeros((number_of_molecules,2)) ## probably not needed

    ## assuming plasmon at orgin: d is equal to the molecule mosition bf_x0, 
    ##... but we need to add a z component to the array
    bf_d = np.hstack((bf_x0,np.zeros((molecule_locations_in_plane.shape[0],1))))

    # Calcualte p0 from coupled equations of motion
    bf_p0 = define_p0_array(p0_hat, plas_orientation, bf_d)
    bf_p1 = define_p1_array(p0_hat, plas_orientation, bf_d)
    
    E0 = mb_p_fields(bf_p0, bf_x0)
    E1 = mb_p_fields(bf_p1, bf_x1)



    mol_images = (np.abs(E0[0])**2.+np.abs(E0[1])**2.+np.abs(E0[2])**2.)
    plas_images = (np.abs(E1[0])**2.+np.abs(E1[1])**2.+np.abs(E1[2])**2.)
    images = (
        np.abs(E1[0]+E0[0])**2.
        +
        np.abs(E1[1]+E0[1])**2.
        +
        np.abs(E1[2]+E0[2])**2.
        )



    Ix = np.abs(E1[0]+E0[0])**2.    
    Iy = np.abs(E1[1]+E0[1])**2.    

    Px = np.sum(Ix,axis=-1)
    Py = np.sum(Iy,axis=-1)

    normed_Px = Px - np.min(Px)*.99
    normed_Py = Py - np.min(Py)*.99

    observed_angles = np.arctan(normed_Py**0.5/normed_Px**0.5)
    
    plt.figure()
    plt.hist(observed_angles,25)
    ticks = [0,np.pi/4,np.pi/2]
    tick_labels = [0, r'$\pi/4$',r'$\pi/2$']
    plt.xticks(ticks, tick_labels)
#     plt.yticks([0,np.pi/16,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2], 
#               [r'$0$', r'$\pi/16$',r'$\pi/8$',r'$\pi/4$',r'$3\pi/8$',r'$\pi/2$'])
    plt.ylabel('molecule images')
    plt.xlabel('actual angle')
    if np.all(plas_orientation==[1,0,0]):
        plas_unit_sym = r'\hat{x}'
    elif np.all(plas_orientation==[0,1,0]):
        plas_unit_sym = r'$\hat{y}$'
    else: plas_unit_sym = '{}'.format(plas_orientation)
    plt.title('d = {0} nm, '.format(dnm)+r'$\hat{p}_\mathrm{plas} = $'+plas_unit_sym)
    plt.show()