%% myfsltopup.m
% myfsltopup.m runs the FSL function topup to create a field map based on
% inputted volumes of different phase-encode directions.
%
% Created by Sam Laxer on December 5, 2022
%
% Usage:
% myfsltopup(mergedVols, acqparamsFile, outfilePath, configFile)
%
% -> mergedVols is a Str with the full file path of the merged file
% containing the different phase-encode acquisitions.
% -> acqparamsFile is a Str with the full file path of the acquired
% parameters file specifying the direction and readout time of the
% different phase-encode volumes.
% -> outfilePath is a Str with the full path (including file name - Don't
% % specify an extension) where the field map and movement parameters
% should be saved.
% -> configFile is a Str with the full file path of the configuration file
% to be used with topup.
%
% Example usage:
% myfsltopup('/Users/samlaxer/Documents/Projects/Project_SL_001/first_level/sub-016/20230317/MRI/T2star_FID_EPI_sat/all_runs_merged/all_runs_merged.nii.gz',...
% '/Users/samlaxer/Documents/Projects/Project_SL_001/first_level/sub-016/20230317/MRI/T2star_FID_EPI_sat/all_runs_merged/acqparams.txt',...
% '/Users/samlaxer/Documents/Projects/Project_SL_001/first_level/sub-016/20230317/MRI/T2star_FID_EPI_sat/all_runs_merged/allruns_topup',...
% '/Users/samlaxer/Documents/Projects/Project_SL_000/blipup_blipdown_distortion_test/sub-009_20221126/trial_20221205/b02b0_human.cnf');
%
% Created by Sam Laxer on Nov 29, 2022


function myfsltopup(mergedVols, acqparamsFile, outfilePath, configFile)
    
    % Get path of output file:
    [outpath,outfilename,~] = fileparts(outfilePath);
    
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
    
    % Run topup command
    % E.x.: topup --imain=merged_volumes_AP_PA.nii.gz --datain=acqparams.txt
    % --config=b02b0_human.cnf --out=topup_humanconfig_nm --iout=unwarpedvolume_human_nm --estmov=0
    status = system([   'topup --imain=', mergedVols, ' --datain=', acqparamsFile,...
                        ' --config=', configFile, ' --out=', outfilename, ' --estmov=0']);
    
    % Make sure that fslmerge ran properly:
    if ~status == 0
        error('topup failed. Check command and try again.');
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
                    '\nAction: Ran topup (through myfsltopup MATLAB function) command with output directory: ', outpath,...
                    '\nCommand: ', ['topup --imain=', mergedVols, ' --datain=', acqparamsFile, ' --config=', configFile, ' --out=', outfilename],...
                    '\nOutput File: ', outfilename, '\n\n'));
    fclose(fid);

end