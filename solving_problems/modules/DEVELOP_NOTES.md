# Notes for development of modules

## 02/07/19 changed param file to 'curly_nrod_vacuum.yaml'

'''Bash
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
'''

next TODO, 
	- change hardcoded plasmon parameters in 'fitting_misLocalization.py' to parameters reference

### ^ finished at 8:17 PM