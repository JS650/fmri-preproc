%% myfslmaths.m
% myfslmaths.m runs the system command fslmaths and also records the
% process in the corresponding AnalysisDoc.txt file.
%
% Usage:
% myfslmaths(workingDir, command)
% -> command is a Str with the command to be passed into fslmaths (i.e.,
% "fslmaths ...command...")
% -> workingDir is a Str with the working directory where the files to be
% worked with are, and where the AnalysisDoc.txt file is.
%
% Example:
% myfslmaths('/Users/me/Desktop', '/Users/me/Desktop/file1.nii -mul
% /Users/me/Desktop/mask.nii /Users/me/Desktop/file1_BrainExtracted.nii')
%
% Created by Sam Laxer on Nov 29, 2022

function myfslmaths(workingDir, command)
    
    % Go to the working directory
    cd(workingDir);

    % Run fslmaths command
    errorStatus = system(['fslmaths ', command]);

    if errorStatus
        error('fslmaths did not run correctly. Check command and try again.')
    end

    %% Update documentation file

    % Check if a .txt file exists with the name AnalysisStepDoc.txt.
    % If not, initiate it:
    if ~exist('AnalysisDoc.txt', 'file')
        fid = fopen('AnalysisDoc.txt', 'a');
        fprintf(fid,    append('Date: ', string(datetime), ...
                        '\nAction: This file was created.',...
                        "\nCommand: fopen('AnalysisDoc.txt', 'a')\n\n"));
        fclose(fid);
        disp('Created AnalysisDoc.txt in current directory.')
    else
        disp('AnalysisDoc.txt already exists in current directory. Not creating.')
    end

    % Update documentation with this command that was just run:
    fid = fopen('AnalysisDoc.txt', 'a');
    fprintf(fid,    append('Date: ', string(datetime), ...
                    '\nAction: Ran fslmaths (through myfslmaths MATLAB function) command with working directory: ', workingDir,...
                    '\nCommand: fslmaths ', command, '\n\n'));
    fclose(fid);

end