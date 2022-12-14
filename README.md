# ElettraSPEM
Elettra Spectromicroscopy Igor Pro routines

## Description
This routines could be used to load and process data measured at Spectromicsorcopy beamline at Elettra Trieste.

![screenshot](/README_Screenshot.png)

## Requiements
Igor Pro from WaveMetrics. Version 7.0 or newer, preferably 64-bit version. (Routines were tested on the 64-bit version of Igor Pro 8.04 and 9.01).

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
3. Check if our procedures are directly in the "...\User Procedures\ElettraSPEM\" directory. If it is necessary remove "-main" form the name or remove additional folder structure that might be created during unpacking.
4. Create an alias/shortcut of Menu_SPEM.ipf file and paste it in Igor Procedures directory (for example):
```
C:\Users\user.name\Documents\WaveMetrics\Igor Pro 9 User Files\Igor Procedures
```
5. After launching Igor Pro the SPEM procedures could be loaded using "Menu -> Data -> Load SPEM"

### HDF5 library on Igor Pro
ElettraSPEM procedures require hdf5 libraries. On Igor Pro 8 those libraries should be already present or if no, follow the instruction included at at the end of this README file.

## Basic usage of SPEM routines

### Load procedures
In Igor Pro menus open: Data -> Load SPEM

### Load data
The data in a form of:
- Spectroimage (SMI),
- 4D image (SMM),
- 2D scan (SMP),
- 1D spectrum (SMP),

Could be loaded using F2 key shortcut or using menu: SPEM -> Load data
In case of 4D images it might take a while.

### Load 3D Map
3D maps are stored in a separate directories as a serries of SMPM files. 
To load this kind of data press F3 key shortcut or select in menu: SPEM -> Load sequence of data.

- You can select any file in D3_X directory.
- All files in that directory are loaded and the data are automatically converted to one or many 3D waves (depending on a type of measurement).
- The software automatically detect the type of scan region or line.
- If the region contains many line scans, only the last one is shown (see below how to open all of them).
- Data are presented as a constant energy maps vs p/nu angles.

### Browse loaded data
If some data were already loaded to Igor but the graph was closed it could be recreated using:
- SPEM -> Browse loaded spectroimages
- SPEM -> Browse loaded 2d maps
- SPEM -> Browse loaded 3d maps
The last browser (3D) might be used also to browse all the areas of a region scan (after the data are loaded only the last line is shown). 

### Additional features
They appear on certain data type. For example "To K" is visible only on a 2D scan.

#### Buttons
- "Del" - remove completely data form Igor (the source hdf5 file is not removed).
- "Add XY" - draw a point at x and z coordinates of selected 2d/3d ARPES measurement.
- "Info" - print information about selected measurement in command window.
- "Spectra" - browse 1D (for spectroimage) or 2D (for 4D image) spectra of.
- "Norm" - look below "ARPES data normalization"
- "To K" - it is available on 2D ARPES spectrum and it allows rough transformation to K space (it is accurate only in the vicinity of the Gamma point).
- "Trans. to K" - this is available on 3D scans and allows the transformation to K space. Some parameters like Fermi Energy and P and T values for normal emission could be set in a dialog window.
- "Browse E" - browse the E(ang) or E(k) cuts of the 3D scans.
- "Save as new plot / Make plot / Make E(k) plot" - all of those allow to make a copy of data shown in the current window and create a separate plot.

#### 1D and 2D spectra browser
For each single pixel in a spectroimage the spectrum of whole energy range of the detector is stored (and also momentum for 4D images). When the "Spectra" button is pressed a new window called BSI_ (Browse spectroimage) is created. The spectra are shown there as red and green curves for two points selected by cursors A and B on a corresponding spectroimage. The spectra are automatically updated when at least one of the cursor is moved. The area of integration is by default set to 1 pixel but could be changed using "Probe area". The spectra are stored in two separated waves called: SMI****_SA and _SB.

If "Select range" is checked the spectroimage is integrated only in the range selected using Ek start/stop (and Ang. start/stop in case of 4D images).


#### ARPES data normalization / detector correction
Select 2D or 3D ARPES data and press "Norm" button. New panel "data normalization" is created. Several options are available.

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

## Transformation to K space
The procedures responsible for the transformation to the K space are stored in SPEM_Kspace.ipf file. 

Three different transformation procedures are available. Two of them (the name starting with 'R') are reverse transformations. This means that the angles P and Nu are calculated for a given point of predefined K space. In case of "Direct" transformation the Kx and Ky for each measured point are calculated directly using P, T, Nu angles (also taking into account possible tilt of the sample). 

1. `R-simple` This transformation is the most efficient and it is used by default. However it is not exact, especially for tilted (cleaved) samples. We recommend to use it only to preview the results.
2. `R-complete` This transformation uses numerical solution and is more exact than the "R-simple".
3. `Direct` This transformation is the most time consuming but should be the most accurate. Unfortunately in case of low count/high noise data it may produce odd looking maps.



### Appendix HDF5 libraries on Igor Pro
1. Make an alias/shortcut for the following file: igor pro folder\more extensions\file loaders\hdf5.xop. 
2. Move the alias or shortcut into the "igor pro folder\igor extensions" folder.
3. Make an alias/shortcut for the following file: igor pro folder\wavemetrics procedures\file input output\hdf5 browser.ipf. 
4. Move the alias or shortcut into the "igor pro folder\igor procedures" folder. 

This activates Igor procedures that implement an hdf5 browser and add a "new hdf5 browser" menu item in the data->load waves menu.

