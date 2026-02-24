%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Specify and Estimate DCM with Session Modulation (Concatenated Data)
%
% Author: Gerasimos Gerardos
% Project: MeMoSLAP Project (P04–Leipzig)
% Paper: Task load modulates network interactions between bilateral fronto-
%        parietal and cerebellar areas during verbal working memory
%
% Purpose:
% This script specifies and estimates subject-level Dynamic Causal Models
% (DCMs) for fMRI data using temporally concatenated sessions. The model
% includes session effects as modulatory inputs and is designed for group-
% level PEB analysis.
%
% Main Steps:
% 1. Load first-level GLM results (SPM.mat) from concatenated sessions.
% 2. Load VOI time series (first eigenvariates) for predefined ROIs.
% 3. Specify a fully connected DCM architecture (A, B, C matrices).
% 4. Model session effects as modulatory inputs (B-matrix).
% 5. Specify driving inputs (C-matrix) for task conditions.
% 6. Create one DCM per subject using spm_dcm_specify.
% 7. Estimate all subject-level DCMs in parallel.
% 8. Combine estimated DCMs into a Group Connectivity Model (GCM).
% 9. Save variance explained and diagnostics.
%
% Model Characteristics:
% - A-matrix: Fully connected intrinsic connectivity.
% - B-matrix: Modulation by session regressors.
% - C-matrix: Driving inputs from task conditions.
% - D-matrix: Disabled (no nonlinear effects).
% - Data: Two sessions concatenated in time.
%
% Intended Use:
% This script is intended for preparing subject-level DCMs for subsequent
% Parametric Empirical Bayes (PEB) group analysis.
% It is NOT intended for VOI extraction or standard task GLM analysis.
%
% Output:
% - Subject-level DCM files in each subject folder.
% - Group-level GCM file.
% - Variance explained summary table.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% Set paths
addpath('/MATLAB/spm12');
spm fmri

data_path = '/results/fMRI_analysis/first_level_simple/';
subj_folders = dir([data_path '/sub-sham40*']);
subj_folders = subj_folders(cellfun(@(x) ~endsWith(x,'04'), {subj_folders.name}));

out_dir = '/results/DCM_models/DCM_sess_mod_restrict_512max/';

DCM_folder = 'concat_for_DCM_sess_mod'; 

analysis_name = 'DCM_analysis_mask_CER_PAR_IFG_ROI_s4fmap_brain_true_concat_sess_mod_restrict_512max';

%% Settings

% MRI scanner settings
TR = 1;   % Repetition time (secs)
TE = 0.0308;  % Echo time (secs) -> mean of first (12) and second (33) echo
mean_center = true; % whether to mean center the input or not (no mean-centering: intrinsic connectivity = baseline; with mean-centering: mean across conds)

% Experiment settings
nsubjects   = 20;
nregions    = 6; 
nconditions = 6;

% Index of each condition in the DCM
zero_back = 1;
one_back = 2;
two_back = 3;
three_back = 4;
session_1 = 5;
session_2 = 6;



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
b = zeros(nregions, nregions, nconditions-1);

% modulatory inputs: 1back, 2back, 3back
b(:,:,session_1) = ones(nregions);
b(:,:,session_2) = ones(nregions);

% turn off modulations of self-connections?
b(1,1,:) = 0;
b(2,2,:) = 0;
b(3,3,:) = 0;
b(4,4,:) = 0;
b(5,5,:) = 0;
b(6,6,:) = 0;
b

% C-matrix
c = zeros(nregions,nconditions-1);

% driving inputs: all trials
c(:,one_back) = 1;
c(:,two_back) = 1;
c(:,three_back) = 1;
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
