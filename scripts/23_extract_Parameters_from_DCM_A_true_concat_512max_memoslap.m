%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Export subject-level intrinsic (A-matrix) connection strengths from PEB/BMA
%
% Author: Gerasimos Gerardos
% Project: MeMoSLAP Project (P04–Leipzig)
% Paper: Task load modulates network interactions between bilateral fronto-
%        parietal and cerebellar areas during verbal working memory
%
% Description:
% This script extracts subject-specific intrinsic connectivity estimates
% (A-matrix parameters) for connections that are identified as reliable at the
% group level in a PEB/BMA analysis of DCMs.
%
% It loads:
% - the estimated GCM (one DCM per subject),
% - the PEB results (to access the reduced/empirical Bayes subject DCMs, RCM),
% - the BMA results (to obtain posterior probabilities for each A-parameter).
%
% For each PEB effect (regressor) in the group model (e.g., "Mean"):
% 1) Select all A-matrix parameters from the BMA parameter list.
% 2) Threshold connections by posterior probability (default Pp > 0.95).
% 3) Parse the surviving parameter labels (e.g., A(3,1)) into directed
%    connections (from-region -> to-region).
% 4) Extract the corresponding per-subject A-parameter estimates from each
%    subject’s reduced DCM (RCM{iSub}.Ep.A).
% 5) Export an Excel file containing subjects × significant connections for
%    that effect.
%
% Output:
% - One Excel file per PEB effect, saved in GCM_folder, containing per-subject
%   A-matrix values for all connections passing the Pp threshold, e.g.:
%     intrinsic_connectivity_A_<EffectName>_<analysis_name>.xlsx
%
% Notes:
% - The region ordering must match the DCM region ordering used during model
%   specification/inversion.
% - The posterior probability threshold (pp_thresh) controls how conservative
%   the selection of “significant” intrinsic connections is.
% - This script focuses on A-matrix parameters only (intrinsic connectivity).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
clear all
close all
clc

%%
addpath('/MATLAB/spm12');

%load GCM
GCM_folder = '/results/DCM_models/DCM_true_concat_512max' 
analysis_name = 'DCM_analysis_mask_CER_PAR_IFG_ROI_s4fmap_brain_true_concat_512max' 

GCM_path = [GCM_folder '/GCM_' analysis_name '.mat'];
GCM = load(GCM_path);
GCM = GCM.GCM;

n_subjects = size(GCM,1);

PEB_file = '/results/DCM_models/DCM_true_concat_512max/PEB_A_B_C_DCM_analysis_mask_CER_PAR_IFG_ROI_s4fmap_brain_true_concat_512max.mat'

dcm_results = load(PEB_file);


RCM = dcm_results.RCM_A_B_C;

BMA_file = '/results/DCM_models/DCM_true_concat_512max/BMA_A_B_C_DCM_analysis_mask_CER_PAR_IFG_ROI_s4fmap_brain_true_concat_512max.mat' 

% Load group-level BMA
load(BMA_file);

% Define region names
regions = {'lIFG', 'rIFG', 'lIPL', 'rIPL', 'lCereb', 'rCereb'};
n_regions = length(regions);



%  Get all A-matrix parameter indices 
A_all_idx = find(contains(BMA_A_B_C.Pnames, 'A('));
n_all_idx = length(BMA_A_B_C.Pnames);
nEffects = length(BMA_A_B_C.Xnames);
nParamsPerEffect = length(A_all_idx);

%  Significance threshold 
pp_thresh = 0.95;

%  Loop through regressors (effects) 
for iEffect = 1:nEffects
    effect_name = BMA_A_B_C.Xnames{iEffect};
    fprintf('\n--- Processing A-matrix effect: %s ---\n', effect_name);

    % Indices for current effect
    idx_start = 1 + (iEffect - 1) * n_all_idx;
    idx_end =  nParamsPerEffect + (iEffect - 1) * n_all_idx;
    A_idx = A_all_idx;
    A_pp = BMA_A_B_C.Pp(idx_start:idx_end);
    A_pnames = BMA_A_B_C.Pnames(A_idx);

    % Significant parameters
    sig_idx = find(A_pp > pp_thresh);
    if isempty(sig_idx)
        fprintf('No significant connections for effect: %s\n', effect_name);
        continue
    end

    % Parse parameter names
    sig_conn_names = {};
    row_extract = [];
    col_extract = [];

    for k = 1:length(sig_idx)
        pname = A_pnames{sig_idx(k)};  % e.g., 'A(3,1)'
        tokens = sscanf(pname, 'A(%d,%d)');
        to = tokens(1); from = tokens(2);
        row_extract(end+1) = to;
        col_extract(end+1) = from;
        conn_label = sprintf('%s_to_%s', regions{from}, regions{to});
        sig_conn_names{end+1} = conn_label;
    end

    % Extract individual values per participant
    nSubs = numel(RCM);
    nConns = length(sig_conn_names);
    subj_params = nan(nSubs, nConns);

    for iSub = 1:nSubs
        for j = 1:nConns
            subj_params(iSub, j) = RCM{iSub}.Ep.A(row_extract(j), col_extract(j));
        end
    end

    % Create table and export
    T = array2table(subj_params, 'VariableNames', sig_conn_names);
    fname = sprintf('%s/intrinsic_connectivity_A_%s_%s.xlsx', GCM_folder, effect_name, analysis_name);
    writetable(T, fname);
    fprintf('Saved: %s\n', fname);
end

disp('Done: One Excel file per effect with per-subject connection values.');
