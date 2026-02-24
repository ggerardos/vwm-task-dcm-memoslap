%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Watershed parcellation of group overlap map + ROI selection (SPM/SS)
%
% Author: Gerasimos Gerardos
% Project: MeMoSLAP Project (P04–Leipzig)
% Paper: Task load modulates network interactions between bilateral fronto-
%        parietal and cerebellar areas during verbal working memory
%
% Description:
% This script performs data-driven fROI definition using a watershed
% parcellation of a group activation-overlap map, and then selects parcels
% (ROIs) based on how consistently they overlap with individual subject
% activation maps.
%
% It uses:
% - SPM12 for NIfTI I/O (spm_vol, spm_read_vols, spm_write_vol)
% - spm_ss for the watershed algorithm (spm_ss_watershed)
%
% Workflow:
% 1) Load a smoothed inter-subject overlap map for a chosen contrast
%    (e.g., 123vs0back) thresholded at a chosen uncorrected p-value.
% 2) Optionally apply an explicit binary mask to restrict parcellation to a
%    predefined region set (e.g., cerebellum + fronto-parietal regions).
% 3) Run the watershed algorithm on the (negated) overlap image, restricted
%    to voxels where at least `overlap_thr_vox` subjects show activation.
%    This yields a labeled parcel image (each ROI has an integer label).
% 4) Save the full watershed parcellation as a single labeled NIfTI file.
% 5) For each subject:
%    - Load the subject’s individual thresholded activation map (same p-thresh).
%    - Compute voxel-level intersection between activated voxels and
%      watershed parcels.
%    - Mark an ROI as “present” in that subject if the overlap contains at
%      least `nVoxels_thresh` activated voxels within the ROI.
% 6) Compute ROI summary metrics:
%    - ROI size (voxels)
%    - Number of subjects overlapping each ROI
%    - Percent of subjects overlapping each ROI
% 7) Export ROI metrics to an Excel sheet and save a QC scatter plot of
%    ROI size vs. % subject overlap.
% 8) Select ROIs that meet a subject-consistency criterion
%    (>= `subject_proportion_thresh` of subjects).
% 9) Write each selected ROI as an individual binary NIfTI mask (one file
%    per ROI) for downstream fROI extraction/analysis.
%
% Key parameters:
% - p_thresh: p-value used to pick the overlap map and individual activation maps
% - smooth_FWHM: smoothing kernel (mm) used for the overlap map filename
% - overlap_thr_vox: minimum number of subjects required at a voxel to include
%   it in the watershed parcellation
% - subject_proportion_thresh: minimum proportion of subjects that must show
%   overlap with an ROI for it to be kept
% - nVoxels_thresh: minimum number of overlapping voxels within a subject to
%   count that subject as overlapping the ROI
% - apply_mask / mask_path: optionally restrict ROI definition to a binary mask
%
% Inputs:
% - Smoothed group overlap map:
%     <overlap_dir>/<contrast>_<p>_smooth<k>mm.nii
% - Individual subject thresholded activation maps:
%     <firstLevel_dir>/<contrast>_<p>.nii
% - (Optional) explicit binary mask NIfTI (mask_path)
%
% Outputs:
% - Labeled watershed parcellation NIfTI (all ROIs in one file)
% - ROI metrics table (.xlsx): size, #subjects overlapping, %overlap
% - Scatter plot PNG: ROI size vs. % subject overlap
% - Individual binary ROI masks for selected parcels (one NIfTI per ROI)
%
% Notes:
% - All images must be in the same space/resolution and aligned voxel-wise.
% - The overlap map is assumed to encode subject counts (higher = more subjects).
% - Update all directory paths (SPM/spm_ss paths, overlap_root, first-level path,
%   and mask_path) to match your environment before running.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% Add paths
addpath '/MATLAB/spm12';

addpath '/MATLAB/spm_ss';


%% Settings

path = '/results/fMRI_analysis/first_level_simple';
subj_folders = dir([path '/sub-*']);
subj_folders = subj_folders(cellfun(@(x) ~endsWith(x,'04'), {subj_folders.name}));

firstLevel_dir = 'combine_ses';
firstLevel_name = 's4fmap_brain_first_level_concatenated';

con_parcel = '123vs0back'; % 123_vs_0backcontrast with which to parcelate (run the watershed algorithm)
con_extract = '123vs0back'; %123_vs_0back contrast with which to extract activated voxels -> used here to compute number of subjects having activation in a parcel

p_thresh = 0.0001; % p-value threshold (uncorrected) for selecting overlap map & computing how many subjects have activation in a parcel

smooth_FWHM = 6; % overlap map smoothing kernel

overlap_thr_vox = 2; % how many subjects should at least have activations in a given voxel for that voxel to be considered for the watershed parcelation

subject_proportion_thresh = 0.8; % proportion of subjects that have overlap with a ROI for it to be included in the parcelation

nVoxels_thresh = 1; % how many voxels have to overlap with a ROI within each subject to count as "some overlap"?

overlap_root = '/results/fROI_Analysis/fROI_p_00001_s4to6fmap_brain_mask/';
overlap_dir = [overlap_root 'individual_activations_overlap'];


out_dir = overlap_dir;

subject_counter = numel(subj_folders);

% Define the mask path
mask_path = '/results/fROI_Analysis/masks/r_brain_s4to6fmap_BINARY_Mask_Cereb_all_SMGantpost_ANG_IFGpo_HOAtlas_cerebAtllas_FLIRT.nii';  % Update this path to your binary mask

apply_mask = 1;

%% Get smoothed overlap map
curr_smoothed_overlap_map_path = [overlap_dir '/' con_parcel '_' num2str(p_thresh) ...
    '_smooth' num2str(smooth_FWHM) 'mm.nii'];

curr_overlap_map = dir(curr_smoothed_overlap_map_path);

if numel(curr_overlap_map) == 1
    curr_overlap_map_path = [curr_overlap_map.folder '/' curr_overlap_map.name];
else
    error('Error: Could not find exactly 1 overlap map!');
end

if apply_mask == 1

    %% Load the mask file
    V_mask = spm_vol(mask_path);
    [Y_mask, XYZ_mask] = spm_read_vols(V_mask);
    
    %% Run Watershed algorithm
    V = spm_vol(curr_overlap_map_path);
    [Y, XYZ] = spm_read_vols(V);
    
    % Apply the binary mask to the overlap map
    Y(~Y_mask) = 0;  % Zero out areas outside the mask

    mask_setting = 'mask_CER_PAR_IFG';
else
    V = spm_vol(curr_overlap_map_path);
    [Y, XYZ] = spm_read_vols(V);

    mask_setting = 'whole_brain';
end

b3 = spm_ss_watershed(-Y, find(Y >= overlap_thr_vox));
n_ROIs = max(b3(:));
fprintf('Done. Defined %d regions\n', n_ROIs);

%% Output all watershed ROIs in 1 nifti file
out_path = [overlap_dir '/Watershed_min' num2str(overlap_thr_vox) 'subj_' ...
    con_parcel '_' num2str(p_thresh)  '_' mask_setting '_smooth' num2str(smooth_FWHM) 'mm.nii'];
a3 = struct('fname', out_path, 'mat', V.mat, 'dim', V.dim, ...
    'dt', [spm_type('float32') spm_platform('bigend')], 'pinfo', [1;0;0]);
spm_write_vol(a3, b3);

%% Compute number of subjects having voxels with each ROI
ROIs_overlap = zeros(n_ROIs, numel(subj_folders));

for iSubject = 1:numel(subj_folders)
    curr_folder = subj_folders(iSubject).name;
    curr_folder_path = [path '/' curr_folder];

    % Get individual activation map (thresholded at p_thresh)
    curr_firstLevel_dir = [curr_folder_path '/' firstLevel_dir];

    indiv_map = dir([curr_firstLevel_dir '/' con_extract '_' num2str(p_thresh) '.nii']);

    if numel(indiv_map) ~= 1
        error(['Error: Not exactly 1 individual thresholded map found at ' ...
            curr_firstLevel_dir '/' con_extract '_' num2str(p_thresh) '.nii']);
    end

    indiv_map_path = [indiv_map.folder '/' indiv_map.name];

    % Load map as matrix
    V_indiv_map = spm_vol(indiv_map_path);
    [Y_indiv_map, XYZ_indiv_map] = spm_read_vols(V_indiv_map);

    % Find voxels that are both part of a watershed-parcel & activated
    % in the individual activation map
    idx_voxelOverlap = intersect(find(b3 ~= 0), find(Y_indiv_map ~= 0));

    ROIs_voxelOverlap = b3(idx_voxelOverlap);

    % In each ROI, check whether number of activated voxels > nVoxels_thresh
    for iROI = 1:n_ROIs
        nVoxels_in_ROI = sum(iROI == ROIs_voxelOverlap);

        if nVoxels_in_ROI >= nVoxels_thresh
            ROIs_overlap(iROI, iSubject) = 1;
        end
    end
end

%% Compute ROI overlap across subjects & ROI sizes
ROIs_overlap_subjects = sum(ROIs_overlap, 2);
ROI_overlap_pct = (ROIs_overlap_subjects ./ iSubject) .* 100;

for iROI = 1:n_ROIs
    ROI_sizes(iROI) = sum(b3(:) == iROI);
end

%% Save ROI sizes as excel-sheet
ROI_sizes_table = table([1:n_ROIs]', ROI_sizes', ROIs_overlap_subjects, ROI_overlap_pct, ...
    'VariableNames', {'ROI', 'n_Voxels', 'Overlap_subjects', 'Overlap_pct'});

writetable(ROI_sizes_table, [overlap_dir '/ROIs_' con_extract ...
    '_Watershed_min' num2str(overlap_thr_vox) 'subj_' con_parcel '_' ...
    num2str(p_thresh) '_' mask_setting '_smooth' num2str(smooth_FWHM) 'mm.xlsx']);

%% Plot ROI size vs. percentage overlap with individual activation maps
scatter_plot = figure; hold on;
scatter(ROI_sizes, ROI_overlap_pct, 25, 'filled');
ylabel('% subjects');
xlabel('ROI size (voxels)');
set(gca, 'ylim', [0 100]);

% Add horizontal lines to plot
x_max = scatter_plot.CurrentAxes.XLim(2);
line([0 x_max], [50 50]);
line([0 x_max], [60 60]);
line([0 x_max], [65 65]);
line([0 x_max], [70 70]);
line([0 x_max], [75 75]);
line([0 x_max], [80 80]);
line([0 x_max], [85 85]);
line([0 x_max], [90 90]);
line([0 x_max], [95 95]);

%% Save the figure
saveas(scatter_plot, [overlap_dir '/ROI_size_VS_Subject_overlap_' con_extract ...
    '_Watershed_min' num2str(overlap_thr_vox) 'subj_' con_parcel '_' ...
    num2str(p_thresh) '_' mask_setting '_smooth' num2str(smooth_FWHM) 'mm.png']);

%% Select ROIs which have overlap with a certain percentage of subjects
subject_thresh = subject_proportion_thresh .* subject_counter;
ROIs_selected = find(ROIs_overlap_subjects >= subject_thresh);
numel(ROIs_selected)

%% Create an individual nifti-image for each selected ROI/parcel
% ROI output directory
outdir_ROIs = [overlap_root mask_setting '_ROI_' firstLevel_name ...
    '/Watershed_min' num2str(overlap_thr_vox) 'subj_' con_parcel '_' ...
    num2str(p_thresh) '_smooth' num2str(smooth_FWHM) 'mm'];

% Create ROI output directory
if ~exist(outdir_ROIs, 'dir')
    mkdir(outdir_ROIs)
end

% Output each ROI
for iROI = 1:numel(ROIs_selected)
    curr_ROI = ROIs_selected(iROI);

    ROI_img = b3;
    ROI_img(b3 == curr_ROI) = 1;
    ROI_img(b3 ~= curr_ROI) = 0;

    ROI_outpath = [outdir_ROIs '/Watershed_min' num2str(overlap_thr_vox)...
        'subj_' con_parcel '_' ...
        num2str(p_thresh) '_' mask_setting '_smooth' num2str(smooth_FWHM) 'mm' ...
        '_ROI_' num2str(curr_ROI) '.nii'];

    V_ROI = struct('fname', ROI_outpath, 'mat', V.mat, 'dim', V.dim, ...
        'dt', [spm_type('float32') spm_platform('bigend')], 'pinfo', [1;0;0]);
    spm_write_vol(V_ROI, ROI_img);
end
