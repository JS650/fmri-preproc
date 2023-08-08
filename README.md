# fmri-preproc
Contains functions to preprocess fMRI data. Functions include motion correction, distortion correction, despiking, and more.

Definitions of each function are shown below:
- myfslmaths.m      - Runs FSL's fslmaths function from within MATLAB.
- myfslmerge.m      - Runs FSL's fslmerge function from within MATLAB. 
- myfsltopup.m      - Runs FSL's topup function from within MATLAB.
- myfslapplytopup.m - Runs FSL's applytopup function from within MATLAB.
- nifti_circshift.m - Aligns a given volume's position in space to a reference volume and circularly shifts the volume so the brain maintains it's original position.

Additionally, there are batch functions that integrate several of the core functions above to complete common tasks. For instance, distortionCorrection_mean.m uses fslmerge, topup, applytopup, fslmaths and nifti_circshift to complete distortion correction, all packed into one function.

### Requirements
Have FSL and AFNI installed and add to MATLAB path.
Do so by adding the following code to the startup.m file:

FSL
```
setenv('PATH', [getenv('PATH') ':<path_to_fsl_binaries>']);
```

AFNI
```
setenv('PATH', [getenv('PATH') ':<path_to_afni_binaries>']);
```

Note, you might have to replace the ':/usr/local/fsl/bin' with the path where fsl was installed on your computer.

For more information on FSL, visit: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSL
For more information on AFNI, visit: https://afni.nimh.nih.gov/
