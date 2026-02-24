%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Extract mean activation (beta/contrast values) from subject fROIs
%
% Author: Gerasimos Gerardos
% Project: MeMoSLAP Project (P04–Leipzig)
% Paper: Task load modulates network interactions between bilateral fronto-
%        parietal and cerebellar areas during verbal working memory
%
% Description:
% This script extracts subject-level activation magnitudes from predefined
% functional ROIs (fROIs) by computing the mean value of first-level contrast
% images within each ROI mask. The ROI masks are subject-specific (e.g., top
% 10% most activated voxels within watershed-derived parcels), and the script
% summarizes mean contrast values per subject, ROI, and contrast.
%
% Workflow:
% 1) Load subject list and ROI directories.
% 2) For each ROI:
%    - For each subject:
%      a) Load the subject-specific ROI mask (binary NIfTI; e.g., top10pct).
%      b) For each selected first-level contrast (e.g., 1back>0back, 2back>0back,
%         3back>0back):
%         - Load the subject’s con_*.nii image from the first-level directory.
%         - Compute the mean contrast value within ROI voxels (Y_ROI == 1),
%           excluding NaN values.
%         - Store results in:
%             * a per-ROI matrix (subjects × contrasts)
%             * a full 3D array (subjects × ROIs × contrasts)
%             * a wide table (one row per subject; one column per ROI×contrast)
% 3) Save:
%    - One .mat file per ROI containing the subjects × contrasts matrix.
%    - A CSV file containing all extracted beta/contrast means across ROIs.
%    - A .mat file storing the contrast name list for reference.
%
% Inputs:
% - Subject-specific ROI masks:
%     ROI_dir/<ROI_folder>/top10pct/<subj>.nii
% - First-level contrast images for each subject:
%     <path>/sub-*/combine_ses/con_*.nii
%
% Outputs:
% - Per-ROI .mat files: output_dir/<ROI_name>.mat   (variable: betas)
% - Combined CSV table: output_dir/betas_all.csv
% - Contrast name list: output_dir/condition_names.mat
%
% Notes:
% - ROI masks and con images must be voxel-aligned (same space/resolution).
% - If a con image is missing or duplicated, the script stops with an error.
% - Update `path`, `ROI_dir`, and `output_dir` to match your environment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all

%% Set paths
addpath('/MATLAB/spm12');
spm fmri

%% Setup
path = '/results/fMRI_analysis/first_level_simple';
subj_folders = dir([path '/sub-*']);
subj_folders = subj_folders(cellfun(@(x) ~endsWith(x,'04'), {subj_folders.name}));


ROI_dir = '/results/fROI_Analysis/fROI_p_00001_s4to6fmap_brain_mask/mask_CER_PAR_IFG_ROI_s4fmap_brain_first_level_concatenated/selected_individual_activations/123_vs_0back';
ROIs = dir([ROI_dir '/Water*']);
ROI_names = {'l_IPL'; 'l_Cer'; 'r_Cer';'r_IPL';'l_IFG';'r_IFG'};

output_dir = '/results/fROI_Analysis/activation_magnitudes/top10pct';

% Check if the directory exists before creating it
if ~exist(output_dir, 'dir')
    mkdir(output_dir)
end

con_numbers = {'con_0005'; 'con_0006'; 'con_0007'};
con_names = {'1back_vs_0back'; '2back_vs_0back'; '3back_vs_0back'};
         
%% Pre-allocate arrays for storing results
nSubjects = numel(subj_folders);
nROIs = numel(ROIs);
nCon = numel(con_names);
betas_all = nan(nSubjects, nROIs, nCon); % Pre-allocating space for betas
betas_table = table('Size', [nSubjects, 2 + nROIs*nCon], 'VariableTypes', repmat({'double'}, 1, 2 + nROIs*nCon));

%% Main Loop
for iROI = 1:nROIs
    curr_ROI = ROIs(iROI);
    curr_ROI_name = curr_ROI.name;
    
           
    betas = nan(nSubjects, nCon); % Pre-allocating beta array per ROI
    
    for iSubject = 1:nSubjects

        
        curr_subj = subj_folders(iSubject);
        curr_subj_name = curr_subj.name;
        curr_subj_path = [curr_subj.folder '/' curr_subj.name];

        curr_ROI_path = [curr_ROI.folder '/' curr_ROI.name '/top10pct/' curr_subj.name '.nii'];
    
        % Load the ROI image
        V_ROI = spm_vol(curr_ROI_path);
        [Y_ROI, XYZ_ROI] = spm_read_vols(V_ROI);

        % Define path to first_level subdirectory for the con images
        first_level_dir = [curr_subj_path '/combine_ses'];

        for iCon = 1:nCon
            curr_con_num = con_numbers{iCon};
            curr_con_name = con_names{iCon};

            % Get the con image from the first_level directory
            con_image = dir([first_level_dir '/' curr_con_num '*nii']);
            if numel(con_image) ~= 1
                error('Error: Could not find exactly 1 con image for subject %s', curr_subj_name);
            end
            con_image_path = [con_image.folder '/' con_image.name];

            % Load the con image
            V_con = spm_vol(con_image_path);
            [Y_con, XYZ_con] = spm_read_vols(V_con);

            % Compute mean beta in ROI, excluding NaN values
            beta_in_ROI = mean(Y_con(Y_ROI == 1 & ~isnan(Y_con)));

            % Store the beta value
            betas(iSubject, iCon) = beta_in_ROI;
            betas_all(iSubject, iROI, iCon) = beta_in_ROI;
            
            % Fill the table with results
            betas_table.Subject(iSubject) = {curr_subj_name};
            betas_table.([ROI_names{iROI} '_' curr_con_name])(iSubject) = beta_in_ROI;
        end
    end
    
    % Save the betas for the current ROI
    save([output_dir '/' ROI_names{iROI} '.mat'], 'betas');
end

% Write the table to CSV and save condition names
writetable(betas_table, [output_dir '/betas_all.csv']);
save([output_dir '/condition_names.mat'], 'con_names');
