# ElettraSPEM
Elettra Spectromicroscopy Igor Pro routines

## Installation

1. Download this repository as zip package:
    - Code -> Downilad ZIP

1. Unpack it into the User Procedures directory.

On MS Windows User Procedures are stored in: 
for example: 'C:\Users\user.name\Documents\WaveMetrics\Igor Pro 9 User Files\User Procedures'

1. Create an alias of Menus_SPEM.ipf file and paste it in Igor Procedures directory:
for example: 'C:\Users\user.name\Documents\WaveMetrics\Igor Pro 9 User Files\Igor Procedures'

1. After launching Igor Pro SPEM procedures could be loaded using "Menu -> Data -> Load SPEM"



### HDF5 library on Igor Pro
ElettraSPEM procedures require hdf5 libraries. On Igor Pro 8 those libraries should be already present if no, follow those instructions:

Install hdf5 library with hdf5 browser on Igor Pro

1. make an alias (mac os x ) or shortcut (windows ) for the following file: igor pro folder\more extensions\file loaders\hdf5.xop. move the alias or shortcut into the "igor pro folder\igor extensions" folder.
2. make an alias (mac os x ) or shortcut (windows ) for the following file: igor pro folder\wavemetrics procedures\file input output\hdf5 browser.ipf. move the alias or shortcut into the "igor pro folder\igor procedures" folder. 

this activates igor procedures that implement an hdf5 browser and add a "new hdf5 browser" menu item in the data->load waves menu.
restart igor pro.



===============================
# basic usage of spem software
===============================

===============================
# load procedures: data -> load new spem
===============================

===============================
# load data: spectro image, 4d image or 2d scan/1d spectrum (smi*, smm*, smp*)

- spem->load data (or press f2 key) 
- choose a file in file browser 
- the data will be automatically loaded and shown (in case of 4d images it might take a while) 

===============================
# load a sequence forming 3d scan

- spem->load sequence of data
- select any file in d3_x directory (the name of the file forming a sequence of 3d scan starts with smpm*)
- all files in that directory are loaded and the data are automatically converted to the 3d wave
- the software automatically detect the type of scan region/line 
- data are presented as a constant energy maps vs p/nu angles

===============================
# browse loaded data

if some data were already loaded to igor but the graph was closed it could be recreated using:
spem-> browse loaded spectro images
spem-> browse loaded 2d maps
spem-> browse loaded 3d maps


===============================
# additional features
===============================

#### spectro image

"del" - remove completely data form igor (the source hdf5 file is not removed)

"add xz" - draws a point at x and z coordinates of selected 2d/3d arpes measurement

"info" - print information about selected measurement in command window

"spectra" - browse 1d (for spectro image) or 2d (for 4d image) spectra


#### arpes data normalization / detector correction

select 2d or 3d arpes data and press norm. button.
new panel "data normalization" is created.

### integration
## first tab allows to adjust integral normalization. 

it is calculated for each nu/p emission angle separately as follows:

new_data = data - mean_of_first_n_points * (background_correction/100)
new_data = data / mean_of_points_over_energy_range ^ (integral_norm/100)

- background num. of points - number of points used for background calculation (should be over the fermi energy).

- background correction (%) - the "power" of this correction. if it is less than 100% the calculated background signal over the fermi energy is reduced by this number, if it is set to 0% no background correction is done.

- norm energy range (%) - 100% means that we use the whole energy range to calculate the integrated signal. 50% means that we use only the half of the points (starting from the lowest energy). 

- integral norm (%) - the "power" of integral normalization, 0% means no correction.


## second tab allows to correct the detector inhomogeneity

first select the wave that will be used for correction.

two types of normalization are available: by image or by 1d profile.
in first case the whole data wave is divided by the detector background wave image. this applies to snap measurements.
second type (by 1d profile) could be applied even for swep scans in wider range of kinetic energy. the detector background image is first converted to 1d profile (summation over energies) and then each line of data is divided by the profile.

warning! 
if during the beamtime the detector active area was changed (different number of energy/angular channels) it might happen that data measured for detector inhomogeneity correction will have slightly different dimensions size. in such a case the user is asked to confirm the reshape procedure of detector correction wave.

before the normalization panel is created the data are copied to separate data directory and in case of a "revert" button is pressed the data are restored.













