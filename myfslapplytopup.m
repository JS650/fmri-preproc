%% myfslapplytopup.m
% myfslapplytopup.m runs the FSL function applytopup to apply the estimated
% field maps from topup to the distorted images for correction.
%
% Created by Sam Laxer on December 7, 2022

% Usage:
% myfslapplytopup(distortedVolsPath, topupFields, APPAindex, acqparamsFile, outfilePath, correctionMethod)
%
% -> distortedVolsPath is a Str with the full file path of the distorted 4D
% file to be corrected (note we will correct all volumes in this function
% rather than just one at a time).
% -> topupFields is a Str with the full path of all the individual
% estimated frames from topup.
% -> APPAindex is a Str representing the index (e.g., 1 or 2) of the
% undistorted files phase-encode direction in the acqparams.txt file.
% -> acqparamsFile is a Str with the full file path of the acquired
% parameters file specifying the direction and readout time of the
% different phase-encode volumes.
% -> correctionMethod is a Str representing the method to be used to
% correct the field (i.e. Jac for correcting individual volumes).
%
% Example input to applytopup:
% applytopup    --imain=/Users/samlaxer/Documents/Projects/Project_SL_001/first_level/sub-007/20221125/MRI/T2star_FID_EPI_sat/run_40001/run_40001_T2star_FID_EPI_sat_20221125091611_40001_d_mc_circshift.nii
%               --topup=/Users/samlaxer/Documents/Projects/Project_SL_001/first_level/sub-007/20221125/MRI/T2star_FID_EPI_sat/merged_mean_run_170001_and_run_40001/topup_fieldmap
%               --datain=/Users/samlaxer/Documents/Projects/Project_SL_001/first_level/sub-007/20221125/MRI/T2star_FID_EPI_sat/merged_mean_run_170001_and_run_40001/acqparams.txt
%               --inindex=2
%               --out=run_40001_T2star_FID_EPI_sat_20221125091611_40001_d_mc_circshift_topupCorrected
%               --method=jac
%
% Created by Sam Laxer on Nov 29, 2022


function myfslapplytopup(distortedVolsPath, topupFields, APPAindex, acqparamsFile, correctionMethod)
    
    % Get path of output file - for applytopup, the output path will be the
    % same as the volume to be corrected:
    [outpath,filename,~] = fileparts(distortedVolsPath);
    % Deal with case where reor_vol has .nii.gz extension:
    if strcmp(filename(end-3:end), '.nii')
        filename = filename(1:end-4);
    end
    outfilename = [filename, '_tc'];
    
    % Go to output path where topup file will be saved. If it doesn't
    % exist, send error:
    if ~exist(outpath, "dir")
        error(['Specified output directory does not exist. Directory: ', outpath]);
    end

    % Go into outpath directory:
    cd(outpath);

    % Make sure topup corrected file doesn't already exist
    if exist([outfilename, '.nii'], 'file') || exist([outfilename, '.nii.gz'], 'file')
        disp([outfilename, ' already exists. Terminating program here to prevent overwrite and wasting time.'])
        return
    end

    % Run applytopup command
    % E.x.: applytopup --imain=run_120001_T2star_FID_EPI_sat_20221126144142_120001_MeanVolume_circshift_voxmatrix.nii
    % --topup=TOPUP_NEWPARAMS --datain=acqparams.txt --inindex=2 --out=APPLYTOPUP_NEWPARAMS12 --method=jac
    status = system([   'applytopup --imain=', distortedVolsPath, ' ',...
                        '--topup=', topupFields, ' ',...
                        '--datain=', acqparamsFile, ' ',...
                        '--inindex=', APPAindex, ' ',...
                        '--out=', outfilename, ' ',...
                        '--method=', correctionMethod]);
    
    % Make sure that applytopup ran properly:
    if ~status == 0
        error('applytopup failed. Check command and try again.');
    end

    %% Update documentation file

    % Check if a .txt file exists with the name AnalysisStepDoc.txt.
    % If not, initiate it:
    if ~exist([outpath, filesep, 'AnalysisDoc.txt'], 'file')
        fid = fopen([outpath, filesep, 'AnalysisDoc.txt'], 'a');
        fprintf(fid,    append('Date: ', string(datetime), ...
                        '\nAction: This file was created.',...
                        "\nCommand: fopen('AnalysisDoc.txt', 'a')\n\n"));
        fclose(fid);
        disp(['Created AnalysisDoc.txt in current directory: ', pwd])
    else
        disp('AnalysisDoc.txt already exists in current directory. Not creating.')
    end
    

    % Update documentation with this command that was just run:
    fid = fopen([outpath, filesep, 'AnalysisDoc.txt'], 'a');
    fprintf(fid,    append('Date: ', string(datetime), ...
                    '\nAction: Ran applytopup (through myfslapplytopup MATLAB function) command with output directory: ', outpath,...
                    '\nCommand: ', [   'applytopup --imain=', distortedVolsPath, ' ',...
                        '--topup=', topupFields, ' ',...
                        '--datain=', acqparamsFile, ' ',...
                        '--inindex=', APPAindex, ' ',...
                        '--out=', outfilename, ' ',...
                        '--method=', correctionMethod],...
                    '\nOutput File: ', outfilename, '\n\n'));
    fclose(fid);

end