%% nifti_circshift
% nifti_circshift(reor_vol, ref_vol) takes in 4D Dataset, reor_vol, and
% circularly shifts it such that it has the same voxel matrix orientation
% as the opposite phase encode volume, ref_vol. It returns the path of the
% outputted file.
%
% reor_vol -> full path of volume to be circularly shifted.
% ref_vol -> full path of volume with opposite phase encode direction used
% as reference and for header information.
%
%
% Created by Sam Laxer on November 29, 2022

% For analyzing volumes with opposite phase encodes, it's important to
% make sure that the sets of volumes are aligned BOTH in terms of their
% matrices in space (defined by header of NIfTI images) and in terms of the
% actual brains as acquired from the scanner. Due to how we have to change
% the phase encode direction on the scanner, we are left with a one-voxel
% shift in the AP and LR dimensions. Thus, we can circularly shift the
% volumes and adjust the header so that everything aligns.

function out_vol = nifti_circshift(reor_vol, ref_vol)

    %topupfiles_path = '/Users/samlaxer/Documents/Projects/Project_SL_000/blipup_blipdown_distortion_test/sub-009_20221126/topupfiles/';
    [workingDir_reor, fileName_reor, ext] = fileparts(reor_vol);
    % Deal with .nii.gz extension
    if strcmp(ext, '.gz')
        fileName_reor = fileName_reor(1:end-4);
    end

    outvol_filepath = [workingDir_reor, filesep, fileName_reor, '_circshift'];

    % If the output file already exists, terminate the program:
    if exist([outvol_filepath, '.nii.gz'], 'file') || exist([outvol_filepath, '.nii'], 'file')
        disp('Output volume already exists. Terminating to avoid overwrite.')
        out_vol = false;
        return
    end

    % Deal with case where reor_vol has .nii.gz extension:
    if strcmp(fileName_reor(end-3:end), '.nii')
        %ext = '.nii.gz';
        fileName_reor = fileName_reor(1:end-4);
    end

    [workingDir_ref, fileName_ref, ~] = fileparts(ref_vol);
    
    % Deal with case where ref_vol has .nii.gz extension:
    if strcmp(fileName_ref(end-3:end), '.nii')
        %ext = '.nii.gz';
        fileName_ref = fileName_ref(1:end-4);
    end


    %% Read in reference NIfTI and NIfTI to be reoriented:
    
    

    % Reference file:
    refvol_info = niftiinfo(ref_vol);
    %refvol_data = niftiread(refvol_info);

    % Pick a random volume
    %rand_vol_idx = round(rand()*size(refvol_data,4));
    %figure; imshow3D(refvol_data(:,:,:,rand_vol_idx)); title('Initial Ref Volume');

    % File to be reoriented:
    reorvol_info = niftiinfo(reor_vol);
    reorvol_data = niftiread(reorvol_info);
    
    %figure; imshow3D(reorvol_data(:,:,:,rand_vol_idx)); title('Initial Reor Volume');
    
    
    %% Change matrix
    
    % NOTE, AS DISCUSSED WITH MARTYN AND RAVI, WE SHOULD BE ABLE TO FIX THE ONE
    % VOXEL ROW OFFSET USING A CIRCULAR SHIFT.
    
    % 2023-07-20:
    %   - There has been several situations where the circular shifts have
    %   shifted in the wrong way. I believe this is because we have not
    %   been taking into account the direction that the data is entered
    %   into the data matrix. We can compensate for this by taking into
    %   account the orientation (e.g., RAS, LPI, etc.)

    circshift_reorvol_data = reorvol_data;

    % Get the ranges along each of the axes:
    [r_status_ref, r_extent_ref] = system(['3dinfo -Rextent ', ref_vol]);
    [l_status_ref, l_extent_ref] = system(['3dinfo -Lextent ', ref_vol]);
    [a_status_ref, a_extent_ref] = system(['3dinfo -Aextent ', ref_vol]);
    [p_status_ref, p_extent_ref] = system(['3dinfo -Pextent ', ref_vol]);
    [i_status_ref, i_extent_ref] = system(['3dinfo -Iextent ', ref_vol]);
    [s_status_ref, s_extent_ref] = system(['3dinfo -Sextent ', ref_vol]);
    r_extent_ref = str2num(r_extent_ref);% When comparing AFNI 3dinfo function and FSLeyes, this values magnitude is the same except it is negative...
    l_extent_ref = str2num(l_extent_ref);% When comparing AFNI 3dinfo function and FSLeyes, this values magnitude is the same except it is negative...
    a_extent_ref = str2num(a_extent_ref);% When comparing AFNI 3dinfo function and FSLeyes, this values magnitude is the same except it is negative...
    p_extent_ref = str2num(p_extent_ref);% When comparing AFNI 3dinfo function and FSLeyes, this values magnitude is the same except it is negative...
    i_extent_ref = str2num(i_extent_ref);
    s_extent_ref = str2num(s_extent_ref);
    
    % Print out position (coordinates) of acquisition volume
    disp('-----------------------------------------')
    disp(['REFERENCE Acquisition Volume Coordinates (', fileName_ref, '):']);
    fprintf(['(Right)\t\t', num2str(r_extent_ref), '\t->\t', num2str(l_extent_ref), '\t(Left)\n']);
    fprintf(['(Anterior)\t', num2str(a_extent_ref), '\t->\t', num2str(p_extent_ref), '\t(Posterior)\n']);
    fprintf(['(Inferior)\t', num2str(i_extent_ref), '\t->\t', num2str(s_extent_ref), '\t(Superior)\n']);
    disp('-----------------------------------------')

    
    [r_status_reor, r_extent_reor] = system(['3dinfo -Rextent ', reor_vol]);
    [l_status_reor, l_extent_reor] = system(['3dinfo -Lextent ', reor_vol]);
    [a_status_reor, a_extent_reor] = system(['3dinfo -Aextent ', reor_vol]);
    [p_status_reor, p_extent_reor] = system(['3dinfo -Pextent ', reor_vol]);
    [i_status_reor, i_extent_reor] = system(['3dinfo -Iextent ', reor_vol]);
    [s_status_reor, s_extent_reor] = system(['3dinfo -Sextent ', reor_vol]);
    r_extent_reor = str2num(r_extent_reor);% When comparing AFNI 3dinfo function and FSLeyes, this values magnitude is the same except it is negative...
    l_extent_reor = str2num(l_extent_reor);% When comparing AFNI 3dinfo function and FSLeyes, this values magnitude is the same except it is negative...
    a_extent_reor = str2num(a_extent_reor);% When comparing AFNI 3dinfo function and FSLeyes, this values magnitude is the same except it is negative...
    p_extent_reor = str2num(p_extent_reor);% When comparing AFNI 3dinfo function and FSLeyes, this values magnitude is the same except it is negative...
    i_extent_reor = str2num(i_extent_reor);
    s_extent_reor = str2num(s_extent_reor);

    % Print out position (coordinates) of acquisition volume
    disp('-----------------------------------------')
    disp(['REORIENT Acquisition Volume Coordinates (', fileName_reor, '):']);
    fprintf(['(Right)\t\t', num2str(r_extent_reor), '\t->\t', num2str(l_extent_reor), '\t(Left)\n']);
    fprintf(['(Anterior)\t', num2str(a_extent_reor), '\t->\t', num2str(p_extent_reor), '\t(Posterior)\n']);
    fprintf(['(Inferior)\t', num2str(i_extent_reor), '\t->\t', num2str(s_extent_reor), '\t(Superior)\n']);
    disp('-----------------------------------------')

    % Get voxel size in all 3 dimensions:

    minLength = 20;

    % Reference voxel size:
    %----------------------
    [vox_i_ref_status, vox_i_ref] = system(['3dinfo -adi ', ref_vol]); % voxel size in first dimension
    [vox_j_ref_status, vox_j_ref] = system(['3dinfo -adj ', ref_vol]); % voxel size in second dimension
    [vox_k_ref_status, vox_k_ref] = system(['3dinfo -adk ', ref_vol]); % voxel size in third dimension
    
    if length(vox_i_ref) > minLength
        vox_i_ref_newline_idxs = regexp(vox_i_ref,'\n');
        vox_i_ref_start_idx = vox_i_ref_newline_idxs(end-1);
        vox_i_ref_num = vox_i_ref(vox_i_ref_start_idx:end);
        vox_i_ref = str2num(vox_i_ref_num);
    else
        vox_i_ref = str2num(vox_i_ref);
    end

    if length(vox_j_ref) > minLength
        vox_j_ref_newline_idxs = regexp(vox_j_ref,'\n');
        vox_j_ref_start_idx = vox_j_ref_newline_idxs(end-1);
        vox_j_ref_num = vox_j_ref(vox_j_ref_start_idx:end);
        vox_j_ref = str2num(vox_j_ref_num);
    else
        vox_j_ref = str2num(vox_j_ref);
    end

    if length(vox_k_ref) > minLength
        vox_k_ref_newline_idxs = regexp(vox_k_ref,'\n');
        vox_k_ref_start_idx = vox_k_ref_newline_idxs(end-1);
        vox_k_ref_num = vox_k_ref(vox_k_ref_start_idx:end);
        vox_k_ref = str2num(vox_k_ref_num);
    else
        vox_k_ref = str2num(vox_k_ref);
    end
    %----------------------

    % Reorient voxel size:
    %----------------------
    [vox_i_reor_status, vox_i_reor] = system(['3dinfo -adi ', reor_vol]); % voxel size in first dimension
    [vox_j_reor_status, vox_j_reor] = system(['3dinfo -adj ', reor_vol]); % voxel size in second dimension
    [vox_k_reor_status, vox_k_reor] = system(['3dinfo -adk ', reor_vol]); % voxel size in third dimension

    if length(vox_i_reor) > minLength
        vox_i_reor_newline_idxs = regexp(vox_i_reor,'\n');
        vox_i_reor_start_idx = vox_i_reor_newline_idxs(end-1);
        vox_i_reor_num = vox_i_reor(vox_i_reor_start_idx:end);
        vox_i_reor = str2num(vox_i_reor_num);
    else
        vox_i_reor = str2num(vox_i_reor);
    end

    if length(vox_j_reor) > minLength
        vox_j_reor_newline_idxs = regexp(vox_j_reor,'\n');
        vox_j_reor_start_idx = vox_j_reor_newline_idxs(end-1);
        vox_j_reor_num = vox_j_reor(vox_j_reor_start_idx:end);
        vox_j_reor = str2num(vox_j_reor_num);
    else
        vox_j_reor = str2num(vox_j_reor);
    end

    if length(vox_k_reor) > minLength
        vox_k_reor_newline_idxs = regexp(vox_k_reor,'\n');
        vox_k_reor_start_idx = vox_k_reor_newline_idxs(end-1);
        vox_k_reor_num = vox_k_reor(vox_k_reor_start_idx:end);
        vox_k_reor = str2num(vox_k_reor_num);
    else
        vox_k_reor = str2num(vox_k_reor);
    end
    %----------------------

    % SANITY CHECK: Check to make sure all the above NIfTI acquisition area positions
    % were successfully saved:
    if r_status_ref || l_status_ref || a_status_ref || p_status_ref || i_status_ref || s_status_ref...
             || r_status_reor || l_status_reor || a_status_reor || p_status_reor || i_status_reor || s_status_reor...
             || vox_i_ref_status || vox_j_ref_status || vox_k_ref_status ||...
             vox_i_reor_status || vox_j_reor_status || vox_k_reor_status
        error('Could not extract full NIfTI information. Check the NIfTI files.')
    end

    % SANITY CHECK: Check to make sure the voxel size is the same between
    % the reference and reoriented volumes:
    if ~(vox_i_ref == vox_i_reor) || ~(vox_j_ref == vox_j_reor)  || ~(vox_k_ref == vox_k_reor) 
        error(['Voxel size differs between the two volumes.',...
            'This program cannot currently correct orientations for runs of different voxel sizes.'])
    end
    

    %b = [zeros(size(a,1),1), a(:, 1:end-1)]
    % Dim1 = 64 voxels = Along saggital axis
    % Dim2 = 32 voxels = Along transverse axis
    % Dim3 = 30 voxels = Along coronal axis
    
    
    % Now that we know the position of the acquired volume in scanner
    % space, we can determine how much and in which direction we need to
    % circularly shift the reverse phase encode volumes so they align with
    % the normal phase encode volumes.

    % Determine how displaced the reor_vol is from the ref_vol in terms of number of voxels:
    r_shift = round((abs(r_extent_reor) - abs(r_extent_ref)) / vox_i_ref); % +ve shifts brain towards right, negative shifts towards left
    a_shift = round((abs(a_extent_reor) - abs(a_extent_ref)) / vox_j_ref); % +ve shifts brain towards anterior, negative shifts towards posterior
    i_shift = round((abs(i_extent_reor) - abs(i_extent_ref)) / vox_k_ref); % +ve shifts brain towards inferior, negative shifts towards superior
    
    % Print out displacement information:
    disp('-----------------------------------------')
    disp('Extent of displacement between the reorient and reference volumes:');
    fprintf(['Reorient volume is currently ', num2str(r_shift), ' voxels to the right of the reference volume.\n']);
    fprintf(['Reorient volume is currently ', num2str(a_shift), ' voxels to anterior of the reference volume.\n']);
    fprintf(['Reorient volume is currently ', num2str(i_shift), ' voxels to inferior of the reference volume.\n']);
    disp('-----------------------------------------')


 
    
    



    %% Get orientation and acquisition information from both datasets
    
    % Get the offsets. These will be used to compare with the other run
    % to see
    % if the other run was diplaced
%     ref_transform = round(refvol_info.Transform.T, 3);
%     ref_qoffset_x = refvol_info.raw.qoffset_x;
%     ref_qoffset_y = refvol_info.raw.qoffset_y;
%     ref_qoffset_z = refvol_info.raw.qoffset_z;

%     reor_transform = round(reorvol_info.Transform.T, 3);
%     reor_qoffset_x = reorvol_info.raw.qoffset_x;
%     reor_qoffset_y = reorvol_info.raw.qoffset_y;
%     reor_qoffset_z = reorvol_info.raw.qoffset_z;
    
    
    % NOTES:
    % When the transform contains negative elements for voxel size
    % i.e.:
    % -0.3000         0         0         0
    %       0   -0.3000         0         0
    %       0         0   -0.5000         0
    %  9.6000    6.5000   14.2500    1.0000
    %
    %   A.K.A. If (qoffset_x > 0) && (qoffset_y > 0) && (qoffset_z > 0)
    % ^^ This orientation is: RAS
    %
    %
    % and in this case:
    %  0.3000         0         0         0
    %       0    0.3000         0         0
    %       0         0   -0.5000         0
    % -9.6000   -2.5000   14.2500    1.0000
    %
    %   A.K.A. If (qoffset_x < 0) && (qoffset_y < 0) && (qoffset_z > 0)
    % ^^ This orientation is: LPS
    %
    %
    % and in this case:
    % -0.3000         0         0         0
    % -0.0000    0.3000         0         0
    %       0         0    0.5000         0
    %  9.6000   -4.5000   -6.2500    1.0000
    %
    %   A.K.A. If (qoffset_x < 0) && (qoffset_y < 0) && (qoffset_z > 0)
    % ^^ This orientation is: RPI
    %
    % R-L -> 64 voxels
    % A-P -> 32 voxels
    % I-S -> 30 voxels

    [status_ref, ref_orientation] = system(['3dinfo -orient ', ref_vol]);
    ref_orientation = ref_orientation(1:3);
    [status_reor, reor_orientation] = system(['3dinfo -orient ', reor_vol]);
    reor_orientation=reor_orientation(1:3);
    % SANITY CHECK - make sure the above system commands ran without error:
    if ~(status_reor == 0) || ~(status_ref == 0)
        error("The '3dinfo -orient ...' command did not run successfully.")
    end
    disp(['Reference volume orientation: ', ref_orientation])
    disp(['Reorient volume orientation: ', reor_orientation])
    

    %% Circularly shift if appropriate

    % circshift(M, -1, 1) to circularly shift M in the SAGGITAL direction
    % (shifting towards the left)
    if r_shift
        % Determine which direction should be circularly shifted:
        if strcmp(reor_orientation(1), 'L')
            circshift_reorvol_data = circshift(circshift_reorvol_data, r_shift, 1);
            disp(['Volume was circularly shifted by ', num2str(r_shift), ' voxels in the 1st dimension.'])
        elseif strcmp(reor_orientation(1), 'R')
            circshift_reorvol_data = circshift(circshift_reorvol_data, -r_shift, 1);
            disp(['Volume was circularly shifted by ', num2str(-r_shift), ' voxels in the 1st dimension.'])
        end
        %figure; imshow3D(circshift_reorvol_data(:,:,:,1)); title('Circ Shift 1st Dim');
    end
    
    % circshift(M, 1, 2) to circularly shift M in the TRANSVERSE direction
    % (shifting towards the anterior)
    if a_shift
        if strcmp(reor_orientation(2), 'P')
            circshift_reorvol_data = circshift(circshift_reorvol_data, a_shift, 2);
            disp(['Volume was circularly shifted by ', num2str(a_shift), ' voxels in the 2nd dimension.'])
        elseif strcmp(reor_orientation(2), 'A')
            circshift_reorvol_data = circshift(circshift_reorvol_data, -a_shift, 2);
            disp(['Volume was circularly shifted by ', num2str(-a_shift), ' voxels in the 2nd dimension.'])
        end
        %figure; imshow3D(circshift_reorvol_data(:,:,:,1)); title('Circ Shift 2nd Dim');
    end
    
    % circshift(M, 1, 3) to circularly shift M in the CORONAL direction
    % (shifting towards the anterior)
    if i_shift
        if strcmp(reor_orientaion(3), 'S')
            circshift_reorvol_data = circshift(circshift_reorvol_data, i_shift, 3);
            disp(['Volume was circularly shifted by ', num2str(i_shift), ' voxels in the 3rd dimension.'])
        elseif strcmp(reor_orientation(3), 'I')
            circshift_reorvol_data = circshift(circshift_reorvol_data, -i_shift, 3);
            disp(['Volume was circularly shifted by ', num2str(-i_shift), ' voxels in the 3rd dimension.'])
        end
        %figure; imshow3D(circshift_reorvol_data(:,:,:,1)); title('Circ Shift 3rd Dim');
    end

    %% Flip data matrix if appropriate

    flip_reorvol_data = circshift_reorvol_data;
    % Flip the data matrix so the raw data matrix is in the same
    % orientation as the reference volume.
    if ~strcmp(reor_orientation(1), ref_orientation(1))
        flip_reorvol_data = flip(flip_reorvol_data, 1); % Flip Left to Right
        disp('Data flipped along 1st dimension.')
    end
    
    if ~strcmp(reor_orientation(2), ref_orientation(2))
        flip_reorvol_data = flip(flip_reorvol_data, 2); % Flip Posterior to Anterior
        disp('Data flipped along 2nd dimension.')
    end

    if ~strcmp(reor_orientation(3), ref_orientation(3))
        flip_reorvol_data = flip(flip_reorvol_data, 3); % Flip Superior to Inferior
        disp('Data flipped along 3rd dimension.')
    end

    
    
    %% Write out new NIfTI:
    % Change header to correspond to matrix (i.e. make header the same as
    % reference volume)
    
    if ~exist([outvol_filepath, '.nii.gz'], 'file') && ~exist([outvol_filepath, '.nii'], 'file')
        % Just change the Sform and Qform matrices and use the orginal
        % image's header info:
        newHeader = reorvol_info;
        newHeader.Transform = refvol_info.Transform;
        niftiwrite(flip_reorvol_data, [outvol_filepath, '.'], newHeader, 'Compressed', true);%reorvol_info);
        
        %% Update documentation file

        % Check if a .txt file exists with the name AnalysisStepDoc.txt.
        % If not, initiate it:
        if ~exist([workingDir_reor, filesep, 'AnalysisDoc.txt'], 'file')
            fid = fopen([workingDir_reor, filesep, 'AnalysisDoc.txt'], 'a');
            fprintf(fid,    append('Date: ', string(datetime), ...
                '\nAction: This file was created.',...
                "\nCommand: fopen('AnalysisDoc.txt', 'a')\n\n"));
            fclose(fid);
            disp('Created AnalysisDoc.txt in current directory.')
        else
            disp('AnalysisDoc.txt already exists in current directory. Not creating.')
        end

        % Update documentation with this command that was just run:
        fid = fopen([workingDir_reor, filesep, 'AnalysisDoc.txt'], 'a');
        fprintf(fid,    append('Date: ', string(datetime),...
            '\nAction: Ran custom MATLAB function nifti_circshift.m with working directory: ', workingDir_reor,...
            '\nCommand: nifti_circshift(', reor_vol, ', ', ref_vol, ')',...
            '\nOutput File: ', [workingDir_reor, filesep, fileName_reor, '_circshift.nii'], '\n\n'));

        fclose(fid);
    
    else
        disp(['FILE ALREADY EXISTS IN CURRENT DIRECTORY. Terminating function to avoid override. File name: ', workingDir_reor, filesep, fileName_reor, '_circshift.nii']);
    end

    % Whether or not the file has already been circular shifted, output the
    % circular shifted volume file path:
    out_vol = [outvol_filepath, '.nii.gz'];


end
