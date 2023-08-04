# fmri-preproc
Contains functions to preprocess fMRI data. Functions include motion correction, distortion correction, despiking, and more.

Definitions of each function are shown below:



Additionally, there are batch functions that integrate several of the core functions above to complete common tasks. For instance, distortionCorrection_mean.m uses fslmerge, topup, applytopup, fslmaths and nifti_circshift to complete distortion correction in one function.

