%% Distortion Correction Script
%
% This method will go through each pair of AP and PA volumes in the
% respective functional runs and correct them using the field generated
% from the MEAN images of both runs (i.e. the same distortion field is used
% to correct all images).
%
%
% PREREQUISITE:
%   - Do NOT include spaces in directory and file names
%
% Created December 6, 2022 by Sam Laxer

% close all
% clear variables

function distortionCorrectionMean(AP_RUN_PATH, PA_RUN_PATH, TOPUP_CONFIG_FILE)

% Add current location to matlab path
filePath = matlab.desktop.editor.getActiveFilename;
path_to_add = fileparts(filePath);
addpath(path_to_add);


%% Read in specified NIfTI files

AP_INFO = niftiinfo(AP_RUN_PATH);
AP_DATA = niftiread(AP_INFO);

PA_INFO = niftiinfo(PA_RUN_PATH);
PA_DATA = niftiread(PA_INFO);


%% Set up stage

NUM_AP_VOLS = AP_INFO.raw.dim(5);
NUM_PA_VOLS = PA_INFO.raw.dim(5);

% % SANITY CHECK: make sure the number of volumes is the same between AP
% % and PA runs:
% if ~(NUM_AP_VOLS == NUM_PA_VOLS)
%     error("The number of volumes does not match between the two runs.")
% end

% Separate the file paths into cell arrays consisting of each directory in
% order
AP_PATH_DIRS = regexp(AP_RUN_PATH, filesep, 'split');
PA_PATH_DIRS = regexp(PA_RUN_PATH, filesep, 'split');

% SANITY CHECK: make sure the runs have the same path up until their
% respective direct parent directories (i.e. the grandparent directory is
% the same directory for both files)
[AP_PATH, AP_FILENAME, AP_EXT] = fileparts(AP_RUN_PATH);
[AP_PARENT_PATH, AP_PARENT_DIR] = fileparts(AP_PATH);
[PA_PATH, PA_FILENAME, ~] = fileparts(PA_RUN_PATH);
[PA_PARENT_PATH, PA_PARENT_DIR] = fileparts(PA_PATH);
if ~strcmp(AP_PARENT_PATH, PA_PARENT_PATH)
    error("The grandparent directory is not the same for the AP and PA files.");
end

% Now that we know the files have the same grandparent directory, change 
% working directory to that directory:
cd(AP_PARENT_PATH)



%% Circular shift both volumes

% Change the data matrix and header of the mean PA volume so it aligns with
% the AP volume
%----------------------------------------------------------------------
% Circularly shift and switch header information
%MEAN_PA_RUN_PATH_CIRCSHIFT = nifti_circshift(MEAN_PA_PATH, MEAN_AP_PATH);
% Circularly shift, flip axes, and change header info for all 400 volumes
% too
PA_RUN_PATH_CIRCSHIFT = nifti_circshift(PA_RUN_PATH, AP_RUN_PATH);
AP_RUN_PATH_CIRCSHIFT = nifti_circshift(AP_RUN_PATH, AP_RUN_PATH); % Shouldn't change anything but creates same file with the _circshift extension

% If nifti_circshift runs and sees that the output file already exists,
% then it will return false. We can thus just continue with the program
% setting the variable to the appropriate variable name:
if islogical(PA_RUN_PATH_CIRCSHIFT)
    PA_RUN_PATH_CIRCSHIFT = append(PA_RUN_PATH(1:end-4), '_circshift.nii.gz');
end
if  islogical(AP_RUN_PATH_CIRCSHIFT)
    AP_RUN_PATH_CIRCSHIFT = append(AP_RUN_PATH(1:end-4), '_circshift.nii.gz');
end

disp(PA_RUN_PATH_CIRCSHIFT)
disp(AP_RUN_PATH_CIRCSHIFT)
[PA_CIRCSHIFT_PATH, PA_CIRCSHIFT_FILENAME, PA_CIRCSHIFT_EXT] = fileparts(PA_RUN_PATH_CIRCSHIFT);
[AP_CIRCSHIFT_PATH, AP_CIRCSHIFT_FILENAME, AP_CIRCSHIFT_EXT] = fileparts(AP_RUN_PATH_CIRCSHIFT);
% Deal with case where extension is .nii.gz
if strcmp(PA_CIRCSHIFT_EXT, '.gz')
    PA_CIRCSHIFT_FILENAME = PA_CIRCSHIFT_FILENAME(1:end-4);
    PA_CIRCSHIFT_EXT = '.nii.gz';
end
if strcmp(AP_CIRCSHIFT_EXT, '.gz')
    AP_CIRCSHIFT_FILENAME = AP_CIRCSHIFT_FILENAME(1:end-4);
    AP_CIRCSHIFT_EXT = '.nii.gz';
end

% Change data matrix orientation so it matches new header *** Not, the
% above function was changed so it completes the matrix orientation too.
% Thus, no need to run nifti_voxmatrix.
%PA_RUN_PATH_VOXMATRIX = nifti_voxmatrix(MEAN_PA_RUN_PATH_CIRCSHIFT);
%PA_RUN_PATH_VOXMATRIX = MEAN_PA_RUN_PATH_CIRCSHIFT;


%% Generate mean volumes from each run

% cd(AP_PARENT_PATH)
cd(AP_PATH)

myfslmaths(AP_CIRCSHIFT_PATH, [AP_CIRCSHIFT_FILENAME, AP_CIRCSHIFT_EXT ' -Tmean ', AP_CIRCSHIFT_FILENAME, '_mean']);
MEAN_AP_CIRCSHIFT_PATH = fullfile(AP_CIRCSHIFT_PATH, filesep, [AP_CIRCSHIFT_FILENAME, '_mean.nii.gz']);
myfslmaths(PA_CIRCSHIFT_PATH, [PA_CIRCSHIFT_FILENAME, PA_CIRCSHIFT_EXT ' -Tmean ', PA_CIRCSHIFT_FILENAME, '_mean']);
MEAN_PA_CIRCSHIFT_PATH = fullfile(PA_CIRCSHIFT_PATH, filesep, [PA_CIRCSHIFT_FILENAME, '_mean.nii.gz']);
cd ..

%% Create new directory to hold merged mean volumes and change that to current working directory
% OUTPATH = [AP_PARENT_PATH, filesep, 'merged_mean_', AP_PARENT_DIR, '_and_', PA_PARENT_DIR];
OUTPATH = [AP_PATH, filesep, 'merged_mean_', AP_PARENT_DIR, '_and_', PA_PARENT_DIR];
if ~exist(OUTPATH, "dir") 
    mkdir(OUTPATH);
    disp(['Created directory: ', OUTPATH]);
else
    disp(['Directory already exists - not creating. Directory: ', OUTPATH]);
end
cd(OUTPATH);

%% Merge the volumes
%----------------------------------------------------------------------
% Create merged files and put them in merged pairs directory
FILE1 = MEAN_AP_CIRCSHIFT_PATH;
FILE2 = MEAN_PA_CIRCSHIFT_PATH;
FILES = [FILE1, ' ', FILE2];
MERGETYPE = '-t';
OUTFILEPATH_MERGE = [OUTPATH, filesep, 'merged_mean_volumes_AP_PA'];
if ~exist([OUTFILEPATH_MERGE, '.nii'], "file") && ~exist([OUTFILEPATH_MERGE, '.nii.gz'], "file")
    myfslmerge(FILES, MERGETYPE, OUTFILEPATH_MERGE);
else
    disp([OUTFILEPATH_MERGE, ' file already exists.'])
end



%% Generate field estimates using topup
% Create acqparams.txt file in OUTFILEPATH directory
if ~exist([OUTPATH, filesep, 'acqparams.txt'], 'file')
    fid = fopen([OUTPATH, filesep, 'acqparams.txt'], 'a');
    fprintf(fid, "0 -1 0 0.1\n0 1 0 0.1\n");
    fclose(fid);

    disp('No acqparams.txt file found. Creating one with default values.');
else
    disp('acqparams.txt file already exists. Using that file for topup.');
end

% Run topup with the merged volumes to get a field map for the mean volumes
MERGEDVOLS = OUTFILEPATH_MERGE;
ACQPARAMS = [OUTPATH, filesep, 'acqparams.txt'];
OUTFILEPATH_TOPUP = [OUTPATH, filesep, 'topup_fieldmap'];
CONFIGFILE = TOPUP_CONFIG_FILE;
if ~exist([OUTFILEPATH_TOPUP, '_fieldcoef.nii'], 'file') && ~exist([OUTFILEPATH_TOPUP, '_fieldcoef.nii.gz'], 'file')
    myfsltopup(MERGEDVOLS, ACQPARAMS, OUTFILEPATH_TOPUP, CONFIGFILE);
else
    disp([OUTFILEPATH_TOPUP, ' file already exists.'])
end


%% Apply the estimated field from the mean images to all the original volumes

DISTORTEDVOLS = {AP_RUN_PATH_CIRCSHIFT, PA_RUN_PATH_CIRCSHIFT};
TOPUPVOL = OUTFILEPATH_TOPUP;
APPAINDEX = {'1', '2'};
ACQPARAMS_APPLYTOPUP = ACQPARAMS;
myfslapplytopup(DISTORTEDVOLS{1}, TOPUPVOL, APPAINDEX{1}, ACQPARAMS_APPLYTOPUP, 'jac')
myfslapplytopup(DISTORTEDVOLS{2}, TOPUPVOL, APPAINDEX{2}, ACQPARAMS_APPLYTOPUP, 'jac')

%% Remove directory added to path at beginning of function

rmpath(path_to_add);


end




