# ElettraSPEM
Elettra Spectromicroscopy Igor Pro routines

## Requiements
Igor Pro from WaveMetrics. Version 8.0 or newer (tested on 64-bit version of Igor Pro 8.04 and 9.01).

For older version of Igor Pro see our old routines: 

https://www.elettra.eu/lightsources/elettra/elettra-beamlines/spectro-microscopy/info-for-users/page-3.html?showall=

## Installation

1. Download this repository as a zip package:
    - Code -> Download ZIP
2. Unpack it into the User Procedures directory.
On MS Windows User Procedures are stored in (for example): 
```
C:\Users\user.name\Documents\WaveMetrics\Igor Pro 9 User Files\User Procedures
```
3. Create an alias/shortcut of Menu_SPEM.ipf file and paste it in Igor Procedures directory (for example):
```
C:\Users\user.name\Documents\WaveMetrics\Igor Pro 9 User Files\Igor Procedures
```
4. After launching Igor Pro the SPEM procedures could be loaded using "Menu -> Data -> Load SPEM"

### HDF5 library on Igor Pro
ElettraSPEM procedures require hdf5 libraries. On Igor Pro 8 those libraries should be already present or if no, follow the instruction included at at the end of this README file.

## Basic usage of SPEM routines

### Load procedures
In Igor Pro menus open: data -> load SPEM

### Load data
The data in a form of:
- Spectroimage (SMI),
- 4D image (SMM),
- 2D scan (SMP),
- 1D spectrum (SMP),

Could be loaded using F2 key or select: SPEM -> Load data
In case of 4d images it might take a while.

### Load 3D Map
3D map are stored in separate directories as a serries of SMPM files. 
To load this kind of data press F3 or select: SPEM -> Load sequence of data.

- You can select any file in D3_X directory.
- All files in that directory are loaded and the data are automatically converted to one or many 3d waves (depending on measurements type).
- The software automatically detect the type of scan region/line.
- If the region contains many line scans only the last is shown (see below how to open all of them).
- Data are presented as a constant energy maps vs p/nu angles.

### Browse loaded data
If some data were already loaded to Igor but the graph was closed it could be recreated using:
- SPEM -> browse loaded spectroimages
- SPEM -> browse loaded 2d maps
- SPEM -> browse loaded 3d maps
The last browser (3D) might be used also to browse all ares of a region scan (after the data are loaded only the last line is snown). 

### Additional features depending on the data type

#### Buttons
- "Del" - remove completely data form Igor (the source hdf5 file is not removed).
- "Add xz" - draws a point at x and z coordinates of selected 2d/3d arpes measurement.
- "Info" - print information about selected measurement in command window.
- "Spectra" - browse 1D (for spectroimage) or 2D (for 4D image) spectra.

##### 1D and 2D spectra browser
For each single pixel in a spectroimage in the data the whole energy (and also momentum for 4D images) spectrum is stored. When the "Spectra" button is pressed a new window called BSI_ (Browse spectroimage) is created. The spectra are shown there as blue and green curves for two points selected by cursors A and B on a corresponding spectroimage. The spectra are automatically updated when at least one of the cursor is moved. The area of integration is by default set to 1 pixel but could be changed using "Probe area". The spectra are stored in two separated waves called: SMI****_SA and _SB.



#### ARPES data normalization / detector correction
Select 2d or 3d arpes data and press norm. button. New panel "data normalization" is created. Several options are available.

#### Integral normalization
First tab allows to adjust the integral normalization. 
It is calculated for each nu/p emission angle separately as follows:

```
new_data = data - mean_of_first_n_points * (background_correction/100)
new_data = data / mean_of_points_over_energy_range ^ (integral_norm/100)
```
Where:
- `background num. of points` - number of points used for background calculation (should be over the fermi energy).
- `background correction (%)` - the "power" of this correction. if it is less than 100% the calculated background signal over the fermi energy is reduced by this number, if it is set to 0% no background correction is done.
- `norm energy range (%)` - 100% means that we use the whole energy range to calculate the integrated signal. 50% means that we use only the half of the points (starting from the lowest energy). 
- `integral norm (%)` - the "power" of integral normalization, 0% means no correction.


#### Correct the detector inhomogeneity
(Second tab in Data normalization panel.
First select the wave that will be used for correction. Two types of normalization are available: by an image or by 1d profile.
In first case the whole data wave is divided by the detector background wave image. This applies best to snap measurements.
Second type (division by 1d profile) could be applied even for swep scans in wider range of kinetic energy. The detector background image is first converted to 1d profile (summation over energies) and then each line of data is divided by the profile.

*****Warning!**
If during the beamtime the detector active area was changed (different number of energy/angular channels) it might happen that data measured for detector inhomogeneity correction will have slightly different dimensions size. In such a case the user is asked to confirm the reshape procedure of detector correction wave.

Before the normalization panel is created the data are copied to separate data directory and in case of a "revert" button is pressed the data are restored.


### Appendix HDF5 libraries on Igor Pro
1. Make an alias/shortcut for the following file: igor pro folder\more extensions\file loaders\hdf5.xop. 
2. Move the alias or shortcut into the "igor pro folder\igor extensions" folder.
3. Make an alias/shortcut for the following file: igor pro folder\wavemetrics procedures\file input output\hdf5 browser.ipf. 
4. Move the alias or shortcut into the "igor pro folder\igor procedures" folder. 

This activates Igor procedures that implement an hdf5 browser and add a "new hdf5 browser" menu item in the data->load waves menu.

