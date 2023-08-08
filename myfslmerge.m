%% myfslmerge.m
% myfslmerge.m runs the system command fslerge and also records the
% process in the corresponding AnalysisDoc.txt file (there likely isn't an
% existing AnalysisDoc.txt file since we're merging to a new directory so
% it will create one).
%
% Usage:
% myfslmerge(files, mergetype, outfilepath)
%
% -> files is a Str with the full file path of the files to be merged (separated by a space).
% -> mergetype is a Str that tells fslmerge how to merge the two files:
%       (from the documentation of fslmerge):
%       '-t' : concatenate images in time
%       '-x' : concatenate images in the x direction
%       '-y' : concatenate images in the y direction
%       '-z' : concatenate images in the z direction
%       '-a' : auto-choose: single slices -> volume, volumes -> 4D (time series)
%       '-tr' : concatenate images in time and set the output image tr to the final option value
%       '-n <N>' : only use volume <N> from each input file (first volume is 0 not 1)
% -> outfilepath is a Str with the full path where the new merged file should be
% saved.
%
% Created by Sam Laxer on Nov 29, 2022

function myfslmerge(files, mergetype, outfilepath)
    
    % Get path of output file:
    [outpath,outfilename,~] = fileparts(outfilepath);
    
    % Output name:
    %outfilename = 'merged_volumes_AP_PA';

    % Go to output path where merged file will be saved. If it doesn't
    % exist, create it:
    if ~exist(outpath, "dir")
        mkdir(outpath);
        disp(['Created directory: ', outpath]);
    else
        disp(['Directory ', outpath, ' already exists. Did not create new directory.']);
    end

    % Go into outpath directory:
    cd(outpath);

    % Run fslmerge command
    status = system(['fslmerge ', mergetype, ' ', [outpath, filesep, outfilename], ' ', files]);
    
    % Make sure that fslmerge ran properly:
    if ~status == 0
        error('fslmerge failed. Check command and try again.');
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
                    '\nAction: Ran fslmerge (through myfslerge MATLAB function) command with output directory: ', outpath,...
                    '\nCommand: ', ['fslmerge ', mergetype, ' ', outfilepath, ' ', files],...
                    '\nOutput File: ', outfilename, '\n\n'));
    fclose(fid);

end