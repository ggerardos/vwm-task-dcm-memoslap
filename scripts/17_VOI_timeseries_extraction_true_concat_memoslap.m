%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: VOI time-series extraction for DCM from concatenated first-level GLMs
%
% Author: Gerasimos Gerardos
% Project: MeMoSLAP Project (P04–Leipzig)
% Paper: Task load modulates network interactions between bilateral fronto-
%        parietal and cerebellar areas during verbal working memory
%
% Description:
% This script extracts VOI (Volume of Interest) time series from subject-level
% SPM12 first-level models that were estimated on a *true concatenated* run
% (nback1 + nback2 combined into a single session). The extracted VOIs are
% intended for downstream DCM analyses.
%
% The script supports two ROI modes:
% - ROI_type = 'group'      : use one common ROI mask for all subjects
% - ROI_type = 'individual' : use subject-specific ROI masks (e.g., top 10%
%                             activated voxels within watershed-derived fROIs)
%
% For each ROI and each subject, the script:
% 1) Loads the subject’s first-level SPM.mat from the concatenated model
%    folder (e.g., <sub>/concat_ses/SPM.mat).
% 2) Defines the ROI mask for VOI extraction (group or individual).
% 3) Runs SPM’s VOI extraction (spm.util.voi), adjusting the VOI time series
%    by the specified Effect-of-Interest (EOI) F-contrast (adjust = EOI_con_number).
% 4) Saves VOI outputs and then renames/moves them into an ROI-specific output
%    directory for clean organization.
%
% Inputs:
% - First-level concatenated SPM models:
%     <path>/sub-*/concat_ses/SPM.mat
% - ROI masks:
%     * Group ROI (if ROI_type = 'group'): a single ROI NIfTI file
%     * Individual ROIs (if ROI_type = 'individual'):
%         ROI_folder/<subject>.nii
% - EOI_con_number:
%     Index of the F-contrast used for “adjusted” VOI time series (SPM VOI adjust)
%
% Outputs (per subject × ROI):
% - VOI_<sub>_<ROI>_1.mat         (VOI structure incl. adjusted eigenvariate)
% - VOI_<sub>_<ROI>_mask.nii      (final ROI mask used by SPM)
% - VOI_<sub>_<ROI>_1_eigen.nii   (eigenvariate image)
% All saved into:
%   outdir/<ROI_name>/
%
% Notes:
% - This script assumes the first-level models were built on a single
%   concatenated session (session = 1 in the VOI batch).
% - ROI masks must be aligned with the subject’s functional data space used
%   in the first-level GLM (same voxel grid).
% - The mask threshold is set to 0.9 to include voxels labeled as 1 in the ROI.
% - Update path variables (SPM path, first-level path, ROI folders, outdir) to
%   match your environment before running.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all

%% Add SPM12 path and initialize
addpath('/MATLAB/spm12');
spm fmri

%% Define paths and variables
path = '/results/fMRI_analysis/first_level_simple';
subj_folders = dir(fullfile(path, 'sub-sham40*'));
subj_folders = subj_folders(cellfun(@(x) ~endsWith(x,'04'), {subj_folders.name}));

first_level_folder = 'concat_ses';

EOI_con_number = 21; % F-contrast number for Effect of Interest


% Choose ROI type: 'group' for one ROI for entire group, 'individual' for individual ROIs
ROI_type = 'individual';

ROI_names = {'left_IPL','right_IPL','left_IFG','right_IFG','left_Cereb_Crus1','right_Cereb_Crus1'}; % Define the name of the ROI
ROI_id = [1,2,7,8,10,11];

for iROI = 1:numel(ROI_id)
    if strcmp(ROI_type, 'group')
        ROI_path = '/data/ROI.nii'; % Path to group ROI
    elseif strcmp(ROI_type, 'individual')
        ROI_folder = ['/results/fROI_Analysis/fROI_p_00001_s4to6fmap_brain_mask/mask_CER_PAR_IFG_ROI_s4fmap_brain_first_level_concatenated/selected_individual_activations/123_vs_0back/Watershed_min2subj_123_vs_0back_0.0001_mask_CER_PAR_IFG_smooth6mm_ROI_' num2str(ROI_id(iROI)) '/top10pct_above0'];
        % Folder containing individual ROIs, assuming named 'SubjectID.nii'
    end
    
    ROI_name = ROI_names{iROI};
    outdir = ['/results/fROI_Analysis/VOI_timeseries/VOI_mask_CER_PAR_IFG_ROI_s4fmap_brain_first_level_true_concat/' ROI_name]; % Output directory for VOIs
    if isempty(dir(outdir))
        mkdir(outdir)
    end
    
    
    %% Loop through subjects and extract VOIs
    for iSubject = 1:numel(subj_folders)
        curr_folder = subj_folders(iSubject).name;
        curr_first_level_path = fullfile(path, curr_folder, first_level_folder);
    
        %%% get SPM.mat
        SPM_mat = dir(fullfile(curr_first_level_path, 'SPM.mat'));
        SPM_mat_path = fullfile(SPM_mat.folder, SPM_mat.name);
    
        %%% get ROI image 
        if strcmp(ROI_type, 'individual')
            % individual ROI for each subject
            ROI_path = fullfile(ROI_folder, [curr_folder '.nii']);
        end
    
        %%% VOI extraction batch
        clear matlabbatch
        matlabbatch{1}.spm.util.voi.spmmat = {SPM_mat_path};
        matlabbatch{1}.spm.util.voi.adjust = EOI_con_number;
        matlabbatch{1}.spm.util.voi.session = 1;
        matlabbatch{1}.spm.util.voi.name = curr_folder;
        matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {ROI_path};
        matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.9; % set to >0.9 as in mask image any included voxel is 1 (setting to 1 excludes the '1' voxels)
        matlabbatch{1}.spm.util.voi.expression = 'i1'; % take all voxels in mask
    
        spm_jobman('run', matlabbatch);
    
        %% Move the output files to outdir
    
        % Define the paths for the current VOI file and the new destination
    VOI_1_file_path = fullfile(curr_first_level_path, ['VOI_' curr_folder '_1.mat']);
    VOI_1_file_path_new = fullfile(outdir, ['VOI_' curr_folder '_' ROI_name '_1.mat']);
    
    % Display the paths for debugging purposes
    fprintf('Original file path: %s\n', VOI_1_file_path);
    fprintf('New file path: %s\n', VOI_1_file_path_new);
    
    % Check if exactly one VOI file exists at the original path
    file_info = dir(VOI_1_file_path);
    if numel(file_info) ~= 1
        error('Error: Could not find exactly 1 VOI_SubjID_1.mat file');
    else
        % Attempt to move the VOI file to the new location
        [status, msg, msgID] = movefile(VOI_1_file_path, VOI_1_file_path_new);
        if status == 0
            % If movefile failed, display the error message
            error('Failed to move file: %s\nMessage ID: %s', msg, msgID);
        else
            fprintf('File moved successfully to %s\n', VOI_1_file_path_new);
        end
    end
    
        % move VOI_SubjID_mask.nii file
        VOI_mask_file_path = [curr_first_level_path '/VOI_' curr_folder '_mask.nii'];
        VOI_mask_file_path_new = [outdir '/VOI_' curr_folder '_' ROI_name '_mask.nii'];
    
        if numel(dir(VOI_mask_file_path)) ~= 1
            error('Error: Could not find exactly 1 VOI_SubjID_mask.nii file')
        end
    
        movefile(VOI_mask_file_path, VOI_mask_file_path_new)
    
    
        % Move VOI_SubjID_1_eigen.nii file
        % Dynamically find the generated VOI file
        VOI_files = dir(fullfile(curr_first_level_path, 'VOI_*_1_eigen.nii')); % Search for VOI files
        if numel(VOI_files) ~= 1
            error(['Error: Could not find exactly 1 VOI file for ' curr_folder]);
        end
    
        VOI_eigen_file_path = fullfile(VOI_files.folder, VOI_files.name);
        VOI_eigen_file_path_new = fullfile(outdir, ['VOI_' curr_folder '_' ROI_name '_1_eigen.nii']);
        
        movefile(VOI_eigen_file_path, VOI_eigen_file_path_new);
    
        pause(2); % Optional pause between subjects
    end
end
