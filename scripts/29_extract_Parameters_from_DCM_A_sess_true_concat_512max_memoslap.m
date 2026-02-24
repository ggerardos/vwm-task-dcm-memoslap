%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Export per-subject intrinsic (A-matrix) connectivity values for
%         connections that are "significant" at the group level (BMA Pp > 0.95)
%
% Purpose:
% After you have estimated a group PEB model and performed BMR/BMA,
% this script:
%   1) Identifies which intrinsic (A) connections have high posterior
%      probability (Pp) in the group BMA (default threshold: 0.95),
%   2) Extracts the corresponding subject-specific A-parameter estimates
%      from the reduced/empirical Bayes DCMs (RCM),
%   3) Saves one Excel file per PEB effect (regressor) containing
%      per-subject values for those connections.
%
% Inputs:
% - GCM_<analysis_name>.mat
%     Cell array with estimated subject DCMs (used mainly for metadata).
% - PEB_A_B_C_<analysis_name>.mat
%     Contains RCM_A_B_C: subject-specific reduced/updated DCMs (RCM).
% - BMA_A_B_C_<analysis_name>.mat
%     Contains BMA_A_B_C: group-level BMA results including:
%       * Pnames: parameter labels (e.g., 'A(3,1)')
%       * Pp: posterior probabilities for each parameter (stacked by effect)
%       * Xnames: effect/regressor names (e.g., 'Mean', 'Group', 'Covariate')
%
% What it does (step-by-step):
% 1) Loads GCM, PEB results (to get RCM), and BMA results.
% 2) Finds the indices of all A-matrix parameters in the BMA parameter list
%    (Pnames containing 'A(').
% 3) For each PEB effect (each column of the design matrix M.X used in PEB):
%    a) Selects the slice of the BMA posterior probability vector (Pp)
%       corresponding to that effect.
%       - Pp is stored as a long vector where parameters are stacked by effect.
%    b) Thresholds A-parameters at pp_thresh (default 0.95).
%    c) Parses each significant parameter name 'A(to,from)' to determine:
%       - "from" region index (column)
%       - "to" region index (row)
%    d) For each subject, extracts the subject-specific intrinsic connectivity
%       estimate from the RCM:
%         RCM{iSub}.Ep.A(to, from)
%    e) Writes an Excel file with columns named like:
%         '<fromRegion>_to_<toRegion>'
%
% Outputs:
% - For each effect, an Excel file:
%     intrinsic_connectivity_A_<effect_name>_<analysis_name>.xlsx
%   saved to GCM_folder.
%
% Notes / Caveats:
% - This script assumes that within BMA_A_B_C, the parameter vector is
%   organised as [all params for effect1; all params for effect2; ...].
% - It also assumes that the first nParamsPerEffect entries in the effect
%   slice correspond to the A parameters *in the same order as A_all_idx*.
%   In most standard PEB outputs this works, but if you ever see mismatches,
%   the safe approach is to index Pp using A_all_idx offset by the effect.
% - The exported subject parameters come from the *RCM* (re-estimated/reduced
%   subject DCMs under the PEB), not necessarily the original GCM DCMs.
%   That is usually what you want for "empirical Bayes" subject estimates.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear all
close all
clc

%%
addpath('/MATLAB/spm12');

%load GCM
GCM_folder = '/results/DCM_models/DCM_sess_mod_restrict_512max'
analysis_name = 'DCM_analysis_mask_CER_PAR_IFG_ROI_s4fmap_brain_true_concat_sess_mod_restrict_512max'

GCM_path = [GCM_folder '/GCM_' analysis_name '.mat'];
GCM = load(GCM_path);
GCM = GCM.GCM;

n_subjects = size(GCM,1);

PEB_file = '/results/DCM_models/DCM_sess_mod_restrict_512max/PEB_A_B_C_DCM_analysis_mask_CER_PAR_IFG_ROI_s4fmap_brain_true_concat_sess_mod_restrict_512max.mat'

dcm_results = load(PEB_file);

RCM = dcm_results.RCM_A_B_C;

BMA_file = '/results/DCM_models/DCM_sess_mod_restrict_512max/BMA_A_B_C_DCM_analysis_mask_CER_PAR_IFG_ROI_s4fmap_brain_true_concat_sess_mod_restrict_512max.mat'

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
