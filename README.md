# fmri-preproc
Contains functions to preprocess fMRI data. Functions include motion correction, distortion correction, despiking, and more.

Definitions of each function are shown below:
- myfslmaths.m      Runs FSL's fslmaths function from within MATLAB.
- myfslmerge.m      Runs FSL's fslmerge function from within MATLAB. 
- myfsltopup.m      Runs FSL's topup function from within MATLAB.
- myfslapplytopup.m Runs FSL's applytopup function from within MATLAB.
- nifti_circshift.m Aligns a given volume's position in space to a reference volume and circularly shifts the volume, so the brain maintains it's

Additionally, there are batch functions that integrate several of the core functions above to complete common tasks. For instance, distortionCorrection_mean.m uses fslmerge, topup, applytopup, fslmaths and nifti_circshift to complete distortion correction, all packed into one function.

### Requirements
Have FSL installed and add to MATLAB path.
Do so by adding the following code to the startup.m file:
```
setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin'])
```


