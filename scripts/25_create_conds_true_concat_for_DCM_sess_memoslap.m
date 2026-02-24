%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Concatenate DCM-GLM condition files and add session modulators
%
% Author: Gerasimos Gerardos
% Project: MeMoSLAP Project (P04–Leipzig)
% Paper: Task load modulates network interactions between bilateral fronto-
%        parietal and cerebellar areas during verbal working memory
%
% Description:
% This script prepares subject-specific condition (.mat) files for a DCM-focused
% first-level GLM where two n-back sessions (nback1 and nback2) are concatenated
% into a single continuous time series. It loads the per-session condition files,
% shifts all session-2 onsets by the length of session 1 (t_concat), and merges
% session-1 and session-2 condition onsets/durations into one combined design.
%
% In addition to the standard task regressors (e.g., 0/1/2/3-back and a pooled
% "all_conditions" regressor), the script explicitly creates two extra regressors:
%   - session_1: events corresponding to task blocks occurring in session 1
%   - session_2: events corresponding to task blocks occurring in session 2
%
% These session regressors are intended to be used as modulatory inputs in a
% subsequent DCM analysis, allowing connectivity modulation (B-matrix) to differ
% between the first and second session.
%
% Workflow:
% 1) Identify all subjects (excluding those ending in '04').
% 2) For each subject:
%    - Load the session-1 and session-2 condition .mat files.
%    - Shift session-2 onsets forward by t_concat (= 1420*TR) to align them in
%      the concatenated timeline.
%    - Concatenate task-condition onsets/durations across sessions.
%    - Create additional "session_1" and "session_2" regressors by combining
%      selected task-condition onsets/durations separately for each session.
%    - Append these session regressors to the names/onsets/durations lists.
%    - Save a new combined condition .mat file for the concatenated GLM.
%
% Inputs:
% - Per-session DCM-GLM condition files:
%     root/*<subject_id>*1.mat   (session 1)
%     root/*<subject_id>*2.mat   (session 2)
%
% Output:
% - Concatenated condition file with added session regressors:
%     outdir/sub-<id>p4_task-all.mat
%   containing variables: names, onsets, durations
%
% Notes:
% - TR, n_Series, and t_concat must match the actual length of each session’s
%   fMRI time series used for concatenation.
% - The construction of session_1/session_2 currently pools specific task
%   regressors (here using indices {2,3,3}); adjust these indices if your
%   condition ordering differs.
% - Ensure outdir exists (or add a mkdir(outdir) block) before saving.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

root = '/results/fMRI_analysis/condition_matfiles_glm_for_dcm';
f_root = '/results/fMRIprep/s4fmap_gr1';
outdir = '/results/fMRI_analysis/concat_conditions/first_level_DCM_sess_mod';

subj_fMRIprep = dir([f_root '/sub-*']);
subj_fMRIprep = subj_fMRIprep(cellfun(@(x) ~endsWith(x,'04'), {subj_fMRIprep.name}));

subj_id = {subj_fMRIprep.name};
subj_id = cellfun(@(x) x(end-3:end), subj_id,'UniformOutput', false);

TR = 1; % repetition time
n_Series = 1420;
t_concat = n_Series * TR;

for iSub = 1:numel(subj_id)
    subj_conds = dir(fullfile(root, ['*' subj_id{iSub} '*']));
    
    conds_1 = subj_conds(cellfun(@(x) endsWith(x,'1.mat'), {subj_conds.name}));
    conds_2 = subj_conds(cellfun(@(x) endsWith(x,'2.mat'), {subj_conds.name}));
    
    c1_mat = load(fullfile(conds_1.folder,conds_1.name));
    c2_mat = load(fullfile(conds_2.folder,conds_2.name));
    
    c2_adj_onsets = cellfun(@(x) x + t_concat , c2_mat.onsets, 'UniformOutput' , false);
    
    c_concat_mat = c1_mat;
    
    for i=1:numel(c_concat_mat.onsets)
        if i ~= 5 
            c_concat_mat.onsets{i} = [c_concat_mat.onsets{i},c2_adj_onsets{i}];
            c_concat_mat.durations{i} = [c_concat_mat.durations{i},c2_mat.durations{i}];
        elseif i == 5
            c_concat_mat.onsets{i} = [c_concat_mat.onsets{i};c2_adj_onsets{i}];
            c_concat_mat.durations{i} = [c_concat_mat.durations{i};c2_mat.durations{i}];
        end
    end
    
    onsets = c_concat_mat.onsets;
   
    on_sess1 = [c1_mat.onsets{2}, c1_mat.onsets{3}, c1_mat.onsets{3}];
    on_sess2 = [c2_adj_onsets{2}, c2_adj_onsets{3}, c2_adj_onsets{3}];

    onsets = [onsets(1:4), on_sess1, on_sess2, onsets(5)];
    
    durations = c_concat_mat.durations;

    dur_sess1 =  [c1_mat.durations{2}, c1_mat.durations{3}, c1_mat.durations{3}];
    dur_sess2 =  [c2_mat.durations{2}, c2_mat.durations{3}, c2_mat.durations{3}];

    durations = [durations(1:4), dur_sess1, dur_sess2, durations(5)];
    
    names = c_concat_mat.names;

    names = [names(1:4), {'session_1'}, {'session_2'}, names(5)];

    ses_name = ['sub-' subj_id{iSub} 'p4_task-all'];
    
    save([outdir '/' ses_name '.mat'], ...
            'names','onsets','durations')
end

aa= load(fullfile(outdir,'sub-4001p4_task-all.mat') )
