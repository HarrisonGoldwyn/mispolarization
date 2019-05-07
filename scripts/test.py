import os
import sys
# module_path = os.path.abspath(os.path.join('..'))
# sys.path.append(module_path)

import numpy as np

project_path = os.path.abspath(os.path.join('..'))
sys.path.append(project_path)

from misloc_mispol_package.calc import BEM_simulation_wrapper as bem
from misloc_mispol_package.calc import fitting_misLocalization as fit

# locations_0, angles_0 = fit.fixed_ori_mol_placement(
#     mol_grid_pts_1D=2,
#     x_max=300,
#     y_max=300,
#     mol_angle=0
#     )

locations_0 = np.array([[300, 300, 0]])
angles_0 = np.array([0])
sim_inst_even_angl = bem.SimulatedExperiment(
    locations_0,
    mol_angle=angles_0,
    )
sim_inst_even_angl.calculate_BEM_fields()

sim_inst_even_angl.plot_mispol_map_wMisloc()

import matplotlib.pyplot as plt
# plt.show()

plt.savefig('test.pdf', dpi=500, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None, metadata=None)