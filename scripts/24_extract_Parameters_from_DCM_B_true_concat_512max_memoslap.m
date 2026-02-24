%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Export subject-level modulatory (B-matrix) connection strengths from PEB/BMA
%
% Author: Gerasimos Gerardos
% Project: MeMoSLAP Project (P04–Leipzig)
% Paper: Task load modulates network interactions between bilateral fronto-
%        parietal and cerebellar areas during verbal working memory
%
% Description:
% This script extracts subject-specific modulatory connectivity estimates
% (B-matrix parameters) for connections that are identified as reliable at the
% group level in a PEB/BMA analysis of DCMs.
%
% It loads:
% - the estimated GCM (one DCM per subject),
% - the PEB results (to access the reduced/empirical Bayes subject DCMs, RCM),
% - the BMA results (to obtain posterior probabilities for each B-parameter).
%
% For each PEB effect (regressor) and each modulatory condition (e.g., 1/2/3-back):
% 1) Select the subset of B-matrix parameters belonging to that condition.
% 2) Threshold parameters by posterior probability (default Pp > 0.95).
% 3) Parse surviving parameter labels into directed connections (from -> to),
%    and prepend the condition name to the exported column labels.
% 4) Extract per-subject B-parameter estimates from each subject’s reduced DCM
%    (RCM{iSub}.Ep.B(:,:,condition_index)).
% 5) Export an Excel file per condition × effect containing subjects ×
%    significant modulatory connections.
%
% Output:
% - Excel files saved in GCM_folder, one per (effect × condition), e.g.:
%     modulatory_connectivity_B_<cond>_<effect>_<analysis_name>.xlsx
%   Each file contains per-subject B-matrix values for all connections passing
%   the Pp threshold for that modulatory condition.
%
% Notes:
% - `conds` defines the modulatory conditions and must match the ordering of
%   modulators in the DCM (B(:,:,k)) and in the PEB/BMA parameterization.
% - `cond_adj` shifts the index used when extracting from Ep.B to account for
%   any preceding modulatory inputs in the model (e.g., if B(:,:,1) is unused).
% - The region ordering must match the DCM region ordering used during model
%   specification/inversion.
% - This script focuses on B-matrix parameters only (modulatory connectivity).
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
B_all_idx = find(contains(BMA_A_B_C.Pnames, 'B('));
n_all_idx = length(BMA_A_B_C.Pnames);
nEffects = length(BMA_A_B_C.Xnames);
conds = {'n1back','n2back','n3back'} %{'n1back','n2back','n3back'}; {'session_1','session_2'}
cond_adj = 1; % conditions perciding modulators
nParamsPerEffect = length(B_all_idx)/length(conds);

%  Significance threshold 
pp_thresh = 0.95;

%  Loop through regressors (effects) 
for iEffect = 1:nEffects
    effect_name = BMA_A_B_C.Xnames{iEffect};
    fprintf('\n--- Processing B-matrix effect: %s ---\n', effect_name);

    
    for iCond = 1:numel(conds)

        cond_name = conds{iCond};
        fprintf('\n--- Processing B-matrix effect: %s ---\n', cond_name);

        % Indices for current effect
        idx_start = B_all_idx(1) + (iEffect - 1) * n_all_idx  + nParamsPerEffect* (iCond - 1);
        idx_end =  B_all_idx(1) + nParamsPerEffect + (iEffect - 1) * n_all_idx + nParamsPerEffect* (iCond - 1)-1;
        cond_start = 1 + nParamsPerEffect * (iCond - 1);
        cond_end = nParamsPerEffect + nParamsPerEffect * (iCond - 1);
        B_idx = B_all_idx(cond_start:cond_end);
        B_pp = BMA_A_B_C.Pp(idx_start:idx_end);
        B_pnames = BMA_A_B_C.Pnames(B_idx);
    
        % Significant parameters
        sig_idx = find(B_pp > pp_thresh);
        if isempty(sig_idx)
            fprintf('No significant connections for effect: %s in %s\n', effect_name,cond_name);
            continue
        end
    
        % Parse parameter names
        sig_conn_names = {};
        row_extract = [];
        col_extract = [];
    
        for k = 1:length(sig_idx)
            pname = B_pnames{sig_idx(k)};  % e.g., 'A(3,1)'
            tokens = sscanf(pname, 'B(%d,%d)');
            to = tokens(1); from = tokens(2);
            row_extract(end+1) = to;
            col_extract(end+1) = from;
            conn_label = sprintf('%s_%s_to_%s', cond_name,regions{from}, regions{to});
            sig_conn_names{end+1} = conn_label;
        end
    
        % Extract individual values per participant
        nSubs = numel(RCM);
        nConns = length(sig_conn_names);
        subj_params = nan(nSubs, nConns);
    
        for iSub = 1:nSubs
            for j = 1:nConns
                subj_params(iSub, j) = RCM{iSub}.Ep.B(row_extract(j), col_extract(j),iCond+cond_adj);
            end
        end
    
        % Create table and export
        T = array2table(subj_params, 'VariableNames', sig_conn_names);
        fname = sprintf('%s/modulatory_connectivity_B_%s_%s_%s.xlsx', GCM_folder, cond_name, effect_name, analysis_name);
        writetable(T, fname);
        fprintf('Saved: %s\n', fname);
    end
end

disp(' Done: One Excel file per effect with per-subject connection values.');
