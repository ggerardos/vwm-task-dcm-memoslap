%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Export per-subject modulatory (B-matrix) connectivity values for
%         connections that are "significant" at the group level (BMA Pp > 0.95)
%         separately for each modulatory condition (here: session_1, session_2)
%
% Purpose:
% After you estimated subject DCMs, built a group PEB on {A,B,C} and ran
% Bayesian Model Reduction / Averaging (BMR/BMA), this script:
%   1) Finds which B-matrix parameters (modulatory connections) have high
%      posterior probability (Pp) at the group level (default: Pp > 0.95),
%   2) Extracts the corresponding subject-specific B-parameter estimates from
%      the reduced/empirical Bayes subject DCMs (RCM),
%   3) Writes one Excel file per PEB effect *and per condition*, containing
%      subject-by-connection values.
%
% Inputs (files you already created earlier):
% - GCM_<analysis_name>.mat
%     Cell array with estimated subject DCMs (mostly used for metadata).
% - PEB_A_B_C_<analysis_name>.mat
%     Contains RCM_A_B_C, i.e., subject-specific reduced/updated DCMs (RCM).
% - BMA_A_B_C_<analysis_name>.mat
%     Contains BMA_A_B_C, i.e., group-level BMA results including:
%       * Pnames : parameter labels (e.g., 'B(3,1,5)')
%       * Pp     : posterior probabilities for each parameter (stacked by effect)
%       * Xnames : effect/regressor names in the PEB design matrix
%
% What it does (high level):
% 1) Loads GCM, RCM, and BMA.
% 2) Identifies the list of all B-parameters in the BMA parameter list
%    (Pnames entries that contain 'B(').
% 3) Assumes the B parameters are ordered by modulatory condition and splits
%    them into equal-sized blocks:
%       nParamsPerEffect = (#B parameters) / (#conditions)
%    For each PEB effect and for each condition:
%       a) Takes the corresponding slice of BMA_A_B_C.Pp for that effect+condition,
%          then thresholds at pp_thresh.
%       b) Parses each significant B name to get row/col indices:
%          B(to,from,cond)
%       c) For each subject, extracts the subject-specific value from:
%          RCM{iSub}.Ep.B(to,from, condIndex)
%          where condIndex is set via (iCond + cond_adj).
%       d) Saves an Excel file with columns named:
%          '<cond>_<fromRegion>_to_<toRegion>'
%
% Outputs:
% - One Excel file per effect and condition, saved in GCM_folder:
%     modulatory_connectivity_B_<cond_name>_<effect_name>_<analysis_name>.xlsx
%
% Key variables you set:
% - conds = {'session_1','session_2'}
%     These are the modulatory conditions you want to export separately.
%
% - pp_thresh = 0.95
%     Threshold for "significant" parameters in the BMA.
%
% - cond_adj = 4
%     IMPORTANT: This is used only when extracting subject-specific Ep.B:
%         Ep.B(:,:, iCond + cond_adj)
%     In your DCM, the third index of B corresponds to the *modulatory input
%     number* as stored in DCM.U (not necessarily the same as your "condition
%     list" length).
%
%     In your setup, you likely have:
%       U(1)=0back, U(2)=1back, U(3)=2back, U(4)=3back,
%       U(5)=session_1, U(6)=session_2
%     and you set:
%       cond_adj = 4  -> iCond=1 -> 5 (session_1), iCond=2 -> 6 (session_2)
%     So cond_adj maps your loop index (1..2) onto the correct U index (5..6).
%
% Notes / Caveats:
% - The indexing of BMA_A_B_C.Pp is the trickiest part: this script assumes
%   the B parameters for each condition are contiguous and equally sized.
%   That is often true if your model has the same B structure for each
%   modulatory input, but if your B structure differs across conditions or
%   SPM changes ordering, the safest way is to index Pp using B_all_idx
%   *within each effect slice* instead of using B_all_idx(1) + offsets.
%
% - The subject values are taken from RCM (empirical Bayes estimates), not
%   from the raw GCM DCMs. That is typically what you want for downstream
%   stats because it conditions subject estimates on the group model.
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
B_all_idx = find(contains(BMA_A_B_C.Pnames, 'B('));
n_all_idx = length(BMA_A_B_C.Pnames);
nEffects = length(BMA_A_B_C.Xnames);
conds = {'session_1','session_2'} %{'n1back','n2back','n3back'}; {'session_1','session_2'}
cond_adj = 4; % conditions perciding modulators
nParamsPerEffect = length(B_all_idx)/length(conds);

%  Significance threshold 
pp_thresh = 0.95;

%  Loop through regressors (effects) 
for iEffect = 1:nEffects
    effect_name = BMA_A_B_C.Xnames{iEffect};
    fprintf('\n Processing B-matrix effect: %s n', effect_name);

    
    for iCond = 1:numel(conds)

        cond_name = conds{iCond};
        fprintf('\n Processing B-matrix effect: %s n', cond_name);

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

disp('Done: One Excel file per effect with per-subject connection values.');
