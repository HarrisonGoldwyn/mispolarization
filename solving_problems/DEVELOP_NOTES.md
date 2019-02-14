# Notes for development of modules

## 02/07/19 changed param file to 'curly_nrod_vacuum.yaml'

```Bash
chair@D-10-19-28-143:~/Documents/Academia/SuperRes/Biteen_colab/Mispolarization/python/gitted/solving_problems$ grep -rn 'modules/.' -e 'curly'
modules/./fitting_misLocalization_adding_noise_to_modeled_images__020619.py:40:stream = open('../parameter_files/curly_param.yaml','r')
modules/./coupled_dipoles.py:29:curly_yaml_file_name = '/curly_param.yaml'
modules/./coupled_dipoles.py:31:    parameter_files_path+curly_yaml_file_name
modules/./coupled_dipoles.py:36:    parameter_files_path+curly_yaml_file_name,'r'
modules/./coupled_dipoles.py:39:# print(curly_yaml_file_name)
modules/./fitting_misLocalization_adding_noise_to_modeled_images__012419.py:37:stream = open('../curly_param.yaml','r')
Binary file modules/./__pycache__/coupled_dipoles.cpython-36.pyc matches
modules/./_fitting_misLocalization_adding_noise_to_modeled_images__012419.py:40:stream = open('../curly_param.yaml','r')
modules/./fitting_misLocalization_adding_noise_to_modeled_images__011619v11.py:34:stream = open('../curly_param.yaml','r')
modules/./curly_quiver_plots.py:34:stream = open('../curly_param.yaml','r')
modules/./spectra.py:14:path_to_yaml = '../curly_param.yaml'
modules/./fitting_misLocalization.py:49:stream = open('../curly_param.yaml','r')
```

next TODO, 
	- change hardcoded plasmon parameters in 'fitting_misLocalization.py' to parameters reference

### ^ finished at 8:17 PM

## 02/08/19
Talked to Curly this morning. Continued to tidy up, merged git repositories that were split for no reason. Reduced prevelance of hard-coded dependencies in the ellipsoid radii. 

## 02/13/19
Back in the office after some snow days. 
- in `fitting_misLocalization.py`;
	- in `fitting_misLocalization.py`, changed the confusingly named variable `	nm` to `m_per_nm`
	- fixed units on the ellipse radii, which __fixed drawing of the ellipse in quiver plots__.
	- added documentation to `class FitModelToData`

### 10:49 AM
Model not doing any better than gaussian for localization in BEM images. Something is wrong, results look very similar to Gaussian mislocalization. Even symmetry axis studied in mislocalization paper look bad. Fits should work here. 

TODO: 
	- in `fit_localization_from_BEM_fields.ipynb`, Should figure out why I am getting complex casting to real errors in 