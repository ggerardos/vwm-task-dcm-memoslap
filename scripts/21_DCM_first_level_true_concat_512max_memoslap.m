%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Specify and estimate subject-level DCMs from concatenated GLMs (GCM)
%
% Author: Gerasimos Gerardos
% Project: MeMoSLAP Project (P04–Leipzig)
% Paper: Task load modulates network interactions between bilateral fronto-
%        parietal and cerebellar areas during verbal working memory
%
% Description:
% This script builds and estimates a set of fMRI Dynamic Causal Models (DCMs)
% for a six-region fronto-parietal–cerebellar network using SPM12. It assumes
% that each subject has already been analyzed with a *concatenated* first-level
% GLM (nback1 + nback2 treated as a continuous run) and that VOI eigenvariates
% have been extracted for each region from that same concatenated model.
%
% Workflow:
% 1) Define the DCM model architecture:
%    - A-matrix: fully connected intrinsic connectivity (all-to-all).
%    - B-matrices: modulatory effects of task load (1-back, 2-back, 3-back)
%      on all connections (optionally excluding self-connections).
%    - C-matrix: driving input applied to all regions via the pooled task
%      regressor ("all_conditions").
%    - D-matrix: disabled (no nonlinear/modulatory gating terms).
%
% 2) For each subject:
%    - Load the subject’s concatenated first-level SPM.mat (DCM_folder).
%    - Load the subject’s ROI time series structures (xY) from VOI .mat files
%      (one per region), and label them consistently.
%    - Call spm_dcm_specify to create a DCM structure with the chosen options
%      (TR/TE, mean-centering choice, max iterations).
%    - Save the subject-specific DCM .mat file and store its path.
%
% 3) Group estimation:
%    - Collate all subject DCMs into a GCM (Group Connectivity Model) object.
%    - Estimate all DCMs using spm_dcm_fit (optionally with parallel workers).
%    - Save the estimated GCM to disk.
%
% 4) Diagnostics:
%    - Run spm_dcm_fmri_check to extract model diagnostics (e.g., variance
%      explained per subject).
%    - Save variance explained values to an Excel sheet for QC/reporting.
%
% Inputs:
% - Concatenated first-level GLM per subject:
%     data_path/sub-*/concat_for_DCM/SPM.mat
% - VOI eigenvariate files per subject and ROI:
%     ROI_dirs{iROI}/*<subj>*mat   (contains xY structure)
%
% Outputs:
% - Subject-level DCM files (unestimated specification):
%     data_path/sub-*/concat_for_DCM/DCM_<analysis_name>.mat
% - Estimated group model collection:
%     out_dir/GCM_<analysis_name>.mat
% - QC table:
%     out_dir/variance_explained_<analysis_name>.xlsx
%
% Key settings:
% - mean_center:
%     true  -> inputs mean-centered (A reflects mean connectivity across conditions)
%     false -> no mean-centering (A reflects implicit baseline connectivity)
% - DCM.options.maxit:
%     maximum number of iterations for model inversion (set to 512 here)
% - use_parfor:
%     estimate DCMs in parallel (set true/false depending on your setup)
%
% Notes:
% - ROI time series and SPM.mat must come from the same concatenated GLM.
% - Condition indices (0/1/2/3-back + all_conditions) must match the design
%   matrix ordering in the DCM GLM.
% - Update paths (SPM, data_path, ROI_dirs, out_dir) before running.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all

%% Set paths
addpath('/MATLAB/spm12');
spm fmri

data_path = '/results/fMRI_analysis/first_level_simple/';
subj_folders = dir([data_path '/sub-sham40*']);
subj_folders = subj_folders(cellfun(@(x) ~endsWith(x,'04'), {subj_folders.name}));

out_dir = '/results/DCM_models/DCM_true_concat_512max/';

DCM_folder = 'concat_for_DCM'; 

analysis_name = 'DCM_analysis_mask_CER_PAR_IFG_ROI_s4fmap_brain_true_concat_512max';

%% Settings

% MRI scanner settings
TR = 1;   % Repetition time (secs)
TE = 0.0308;  % Echo time (secs) -> mean of first (12) and second (33) echo
mean_center = true; % whether to mean center the input or not (no mean-centering: intrinsic connectivity = baseline; with mean-centering: mean across conds)

% Experiment settings
nsubjects   = 20;
nregions    = 6; 
nconditions = 5;

% Index of each condition in the DCM
zero_back = 1;
one_back = 2;
two_back = 3;
three_back = 4;
all_conditions = 5;


% Index of each region in the DCM
l_IFG = 1; 
r_IFG = 2; 
l_IPL = 3; 
r_IPL = 4; 
l_Cereb_sup = 5; 
r_Cereb_sup = 6;

% directories of ROI images (first eigenvariate of ROI timeseries)
ROI_dirs = {'/results/fROI_Analysis/VOI_timeseries/VOI_mask_CER_PAR_IFG_ROI_s4fmap_brain_first_level_true_concat/left_IFG';
            '/results/fROI_Analysis/VOI_timeseries/VOI_mask_CER_PAR_IFG_ROI_s4fmap_brain_first_level_true_concat/right_IFG';
            '/results/fROI_Analysis/VOI_timeseries/VOI_mask_CER_PAR_IFG_ROI_s4fmap_brain_first_level_true_concat/left_IPL';
            '/results/fROI_Analysis/VOI_timeseries/VOI_mask_CER_PAR_IFG_ROI_s4fmap_brain_first_level_true_concat/right_IPL';
            '/results/fROI_Analysis/VOI_timeseries/VOI_mask_CER_PAR_IFG_ROI_s4fmap_brain_first_level_true_concat/left_Cereb_Crus1';
            '/results/fROI_Analysis/VOI_timeseries/VOI_mask_CER_PAR_IFG_ROI_s4fmap_brain_first_level_true_concat/right_Cereb_Crus1'};

ROI_names = {'l_IFG';'r_IFG';'l_IPL'; 'r_IPL'; 'l_Cereb_sup'; 'r_Cereb_sup'};

%% Specify DCM architecture (full model)

% A-matrix
a = ones(nregions,nregions); % fully connected

% B-matrix
b = zeros(nregions, nregions, nconditions);

% modulatory inputs: 1back, 2back, 3back
b(:,:,one_back) = ones(nregions);
b(:,:,two_back) = ones(nregions);
b(:,:,three_back) = ones(nregions);

% turn off modulations of self-connections?
b(1,1,:) = 0;
b(2,2,:) = 0;
b(3,3,:) = 0;
b(4,4,:) = 0;
b(5,5,:) = 0;
b(6,6,:) = 0;
b

% C-matrix
c = zeros(nregions,nconditions);

% driving inputs: all trials
c(:,all_conditions) = 1;
c

% D-matrix
d = zeros(nregions,nregions,0); % disable

%% Loop over subjects & specify 1 DCM per subject
start_dir = pwd;
subj_counter = 1;
DCM_paths = {};

for iSubject = 1:numel(subj_folders)
    
    curr_folder = subj_folders(iSubject).name;
    curr_folder_path = [data_path '/' curr_folder]; 
    
    %%
    curr_DCM_folder = [curr_folder_path '/' DCM_folder];

    % Load SPM
    clear SPM
    SPM = load([curr_DCM_folder '/SPM.mat']);
    SPM = SPM.SPM;

   
 
% Load ROIs    
        for iROI = 1:numel(ROI_dirs)
            
            clear curr_ROI XY
            curr_ROI = dir([ROI_dirs{iROI} '/*' curr_folder '*mat']);
            
            if numel(curr_ROI) ~= 1
                error(['Error: Not exactly 1 ROI file found at ' ROI_dirs{iROI} '/*' curr_folder '*mat'])
            else
                curr_ROI_path = [curr_ROI.folder '/' curr_ROI.name];
            end
            
            XY = load(curr_ROI_path);
            %xY.name = ROI_names{iROI};
            xY(iROI) = XY.xY;
            xY(iROI).name = ROI_names{iROI};
        end

% Display loaded ROI names
disp('Loaded ROI names:');
for iROI = 1:numel(xY)
    disp(xY(iROI).name);
end 

    % Move to output directory
    cd(curr_DCM_folder)

    % Select whether to include each condition from the design matrix
    include = ones(nconditions, 1); % include all regressors    

    % Specify. Corresponds to the series of questions in the GUI.
    s = struct();
    s.name       = analysis_name;
    s.u          = include;                 % Conditions
    s.delays     = repmat(TR,1,nregions);   % Slice timing for each region
    s.TE         = TE;
    s.nonlinear  = false;
    s.two_state  = false;
    s.stochastic = false;
    s.centre     = mean_center; % mean-center the input? if yes, A-matrix = mean connectivity across conds; if no, A-matrix = connectivity of implicit baseline
    s.induced    = 0;
    s.a          = a;
    s.b          = b;
    s.c          = c;
    s.d          = d;
    DCM = spm_dcm_specify(SPM,xY,s);
    DCM.options.maxit = 512;

    save([curr_DCM_folder '/DCM_' s.name '.mat'])

    DCM_paths{subj_counter} = [curr_DCM_folder '/DCM_' s.name '.mat'];


    subj_counter = subj_counter + 1;

    % Return to script directory
    cd(start_dir);

end

DCM_paths = DCM_paths';

%% Collate into a GCM file and estimate

% Filenames -> DCM structures
GCM = spm_dcm_load(DCM_paths);

% Estimate DCMs (this won't affect original DCM files)
use_parfor = true;
GCM = spm_dcm_fit(GCM, use_parfor);

% Save estimated GCM
save([out_dir '/GCM_' analysis_name '.mat'],'GCM');

% diagnostics (e.g. variance explained)
DCM_check = spm_dcm_fmri_check(GCM)

varExp = NaN(numel(DCM_check),1);
for iSub = 1:numel(DCM_check)
    varExp(iSub) = DCM_check{iSub}.diagnostics(1);
end

varExp_table = table(varExp,'VariableNames',{'Variance_explained'});

writetable(varExp_table, [out_dir '/variance_explained_' analysis_name '.xlsx'])
