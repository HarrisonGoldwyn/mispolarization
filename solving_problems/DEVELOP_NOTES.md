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

## 02/14/19
Might have just fixed the bug! had duplicates of eps_inf in `curly_nrod_vacuum.yaml`.

## 02/20/19: Moving back to water
First thing I want to do is change param file back to a default name to avoid confusion. 

This is going to involve grepping 
```Bash
chair@Harrisons-MacBook:~/Documents/Academia/SuperRes/Biteen_colab/Mispolarization/python/gitted$ grep -rn 'solving_problems/.' -e 'curly_nrod_vacuum'
solving_problems/./DEVELOP_NOTES.md:3:## 02/07/19 changed param file to 'curly_nrod_vacuum.yaml'
solving_problems/./DEVELOP_NOTES.md:43:Might have just fixed the bug! had duplicates of eps_inf in `curly_nrod_vacuum.yaml`.
solving_problems/./modules/coupled_dipoles.py:31:curly_yaml_file_name = '/curly_nrod_vacuum.yaml'
Binary file solving_problems/./modules/__pycache__/coupled_dipoles.cpython-36.pyc matches
Binary file solving_problems/./modules/__pycache__/fitting_misLocalization.cpython-36.pyc matches
solving_problems/./modules/fitting_misLocalization.py:65:curly_yaml_file_name = '/curly_nrod_vacuum.yaml'
solving_problems/./notebooks/fit_localization_from_BEM_fields.ipynb:24:      "reading parameters from /Users/chair/Documents/Academia/SuperRes/Biteen_colab/Mispolarization/python/gitted/parameter_files/curly_nrod_vacuum.yaml\n",
solving_problems/./notebooks/fit_localization_from_BEM_fields.ipynb:25:      "reading parameters from /Users/chair/Documents/Academia/SuperRes/Biteen_colab/Mispolarization/python/gitted/parameter_files/curly_nrod_vacuum.yaml\n",
solving_problems/./notebooks/fit_localization_from_BEM_fields.ipynb:63:    "curly_yaml_file_name = '/curly_nrod_vacuum.yaml'\n",
solving_problems/./notebooks/fit_ellipsoid_parameters_from_spectra.ipynb:199:      "reading parameters from /Users/chair/Documents/Academia/SuperRes/Biteen_colab/Mispolarization/python/gitted/parameter_files/curly_nrod_vacuum.yaml\n",
solving_problems/./notebooks/.ipynb_checkpoints/fit_localization_from_BEM_fields-checkpoint.ipynb:24:      "reading parameters from /Users/chair/Documents/Academia/SuperRes/Biteen_colab/Mispolarization/python/gitted/parameter_files/curly_nrod_vacuum.yaml\n",
solving_problems/./notebooks/.ipynb_checkpoints/fit_localization_from_BEM_fields-checkpoint.ipynb:25:      "reading parameters from /Users/chair/Documents/Academia/SuperRes/Biteen_colab/Mispolarization/python/gitted/parameter_files/curly_nrod_vacuum.yaml\n",
solving_problems/./notebooks/.ipynb_checkpoints/fit_localization_from_BEM_fields-checkpoint.ipynb:63:    "curly_yaml_file_name = '/curly_nrod_vacuum.yaml'\n",
solving_problems/./notebooks/.ipynb_checkpoints/fit_ellipsoid_parameters_from_spectra-checkpoint.ipynb:199:      "reading parameters from /Users/chair/Documents/Academia/SuperRes/Biteen_colab/Mispolarization/python/gitted/parameter_files/curly_nrod_vacuum.yaml\n",
solving_problems/./notebooks/.ipynb_checkpoints/fit_localization_from_BEM_fields__true_pos_ini_model_guess-checkpoint.ipynb:24:      "reading parameters from /Users/chair/Documents/Academia/SuperRes/Biteen_colab/Mispolarization/python/gitted/parameter_files/curly_nrod_vacuum.yaml\n",
solving_problems/./notebooks/.ipynb_checkpoints/fit_localization_from_BEM_fields__true_pos_ini_model_guess-checkpoint.ipynb:25:      "reading parameters from /Users/chair/Documents/Academia/SuperRes/Biteen_colab/Mispolarization/python/gitted/parameter_files/curly_nrod_vacuum.yaml\n",
solving_problems/./notebooks/.ipynb_checkpoints/fit_localization_from_BEM_fields__true_pos_ini_model_guess-checkpoint.ipynb:63:    "curly_yaml_file_name = '/curly_nrod_vacuum.yaml'\n",
solving_problems/./notebooks/BEM_fit_nanrod_in_vacuum__012419.ipynb:43:    "stream = open('../curly_nrod_vacuum.yaml','r')\n",
```
So the files to find and replace are 
- modules/coupled_dipoles.py
- modules/fitting_misLocalization.py
- notebooks/fit_localization_from_BEM_fields.ipynb
and maybe 
- fit_ellipsoid_parameters_from_spectra.ipynb

But before I go about changing any of that, I need to find a spectra of the rod in water and try and fit it. 

Found the file `~/Documents/MATLAB/102418/curly_rod_spectra_JC_epsb1p778.mat`, will try loading and fitting that.

Resultsing fit parameters without changing anything but the sprectra file are:

```
	old: array([15.74962612, 10.02510511,  0.10268783, 67.16871096, 20.79635188])
	new: array([15.74962612, 10.02510511,  0.10268783, 67.16871096, 20.79635188])
```
Nothing changed. Let me look back at `coupled_dipoles.py` and see what's up.

### 12:26 PM
`eps_b` was being called directly in `fit_ellipsoid_parameters_from_spectra_in_water.ipynb` but thats fixed now. More reasonable parameters except for eps_inf:
	array([31.72908818, 14.21215911,  0.08327776, 55.42533033, 16.72842459])
Not sure what that is about. `coupled_dipoles.py` seems clean now, all `k`'s are now `= w*n/c`. 

Might have to look in MNPBEM and see how they calculate the scattered power with a dielectric background. I'll do that after lunch. 
- if that doesnt work I'll rerun the spectra (with Drude) just to be sure.

### 01:28 PM 
Looking through the BEM source code, it looks like the reported __power scattered by a plane wave is just the integral of vacuum Poynting vector (E x B) divided by n_b__...

from file `/Users/chair/Documents/MATLAB/MNPBEM17/Simulation/retarded/@planewaveret/scattering.m`:
```Matlab. 
function [ sca, dsca ] = scattering( obj, sig )
%  SCATTERING - Scattering cross section for plane wave excitation.
%
%  Usage for obj = planewaveret :
%    [ sca, dsca ] = scattering( obj, sig )
%  Input
%    sig        :  compstruct object containing surface currents
%  Output
%     sca       :  scattering cross section
%    dsca       :  differential cross section

%  total and differential radiated power
[ sca, dsca ] = scattering( obj.spec, sig );

%  refractive index of embedding medium
nb = sqrt( sig.p.eps{ 1 }( sig.enei ) );
%  The scattering cross section is the radiated power normalized to the
%  incoming power, the latter being proportional to 0.5 * epsb * (clight / nb)
[ sca, dsca.dsca ] = deal( sca / ( 0.5 * nb ), dsca.dsca / ( 0.5 * nb ) );
```
where the scattered power is divided by ( 0.5 * nb ).
where it seems to reference the higher up file `/Users/chair/Documents/MATLAB/MNPBEM17/Simulation/retarded/scattering.m`
```matlab
function [ sca, dsca ] = scattering( field, medium )
%  SCATTERING - Radiated power for electromagnetic fields.
%
%  Usage  :
%    [ sca, dsca ] = scattering( field )
%  Input
%    field      :  compstruct object with scattered electromagnetic fields
%    medium     :  compute total radiated pwoer only in given medium
%  Output
%    sca        :  total radiated power
%    dsca       :  differential radiated power for particle surface
%                  given in FIELD

%  particle surface for fields at 'infinity'
pinfty = field.p;
%  scattered electric and magnetic fields
e = field.e; 
h = field.h; 
%  Poynting vector in direction of outer surface normal
dsca = 0.5 * real( inner( pinfty.nvec, cross( e, conj( h ), 2 ) ) );
  
%  area of boundary elements
area = pinfty.area;
%  total cross section only in given medium
if exist( 'medium', 'var' ) && ~isempty( medium )
  area( pinfty.expand( num2cell( pinfty.inout( :, end ) ) ) == medium ) = 0;
end
%  total radiated power
sca = squeeze( matmul( reshape( area, 1, [] ), dsca ) );
%  convert differential radiated power to compstruct object
dsca = compstruct( pinfty, field.enei, 'dsca', dsca );

```
which seems to define the scattered power as the surface integral of 
$$ \mathbf{S} = 1/2 \mathbf{E} \times \mathbf{B}$$
in line 20.

### 02:19 PM
Fitting a Drude spectra in water, I get much more reasonable fit results:
__Drude fits__
```
	water: array([11.6817827 ,  9.56022562,  0.08554153, 45.24119972, 18.7831301 ])
	vacuu: array([11.5936854 ,  9.65413421,  0.06807552, 48.04291550, 19.9337624 ])
```
A little better actually (particularly the radii), probably because the trans peaks split.  

JC still looks bad, but maybe its passable
__JC fits__
```
	water: array([28.98570393, 13.58850188,  0.09766546, 51.77701186, 15.96014311])
```
and this is all with $k = \omega n_b / c$ and 
$$ \sigma_{sca} = (8* \pi / 3) * (\omega * n_b / c)^4 /(0.5 * n_b) |\alpha|^2 $$

## 02/21/19

Water Drude long mode peak ~ 1.85 eV 
	short peak ~ 2.49 eV

Water JC long mode peak ~ 1.86 eV
	short peak ~ 2.4 eV

### 07:38 PM 
Note getting good localization fits for drude simulations in water. Looking for missing factors of n_b but they seem to be in the right places for 
- fitting_misLocalization.py 
- coupled_dipoles.py
- anal_foc_diff_fields.py
Note sure what else could be wrong. It looks like the same error I had last week in vacuum, which I fixed by removing a duplicate input of eps_inf I think. * yup, see note above under 02/14/19 heading. 
