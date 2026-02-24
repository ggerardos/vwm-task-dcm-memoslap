%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Concatenate session condition .mat files for combined first-level GLM
%
% Author: Gerasimos Gerardos
% Project: MeMoSLAP Project (P04–Leipzig)
% Paper: Task load modulates network interactions between bilateral fronto-
%        parietal and cerebellar areas during verbal working memory
%
% Description:
% This script concatenates two session-specific SPM multi-condition .mat files
% (e.g., nback1 and nback2) into a single “combined-session” condition file
% for each subject. It adjusts all onsets from the second session by adding
% the duration of the first session (in seconds), so that event timings are
% expressed in one continuous time axis (as if the two runs were back-to-back).
% The resulting .mat files (names, onsets, durations) can then be used for
% first-level GLM models that treat both sessions as a single concatenated run.
%
% Workflow (per subject):
% 1) Identify the subject’s two condition files in `root`:
%    - Session 1 condition file (*1.mat)
%    - Session 2 condition file (*2.mat)
% 2) Load both .mat files (expects variables: names, onsets, durations).
% 3) Compute the concatenation offset:
%    - t_concat = n_Series * TR  (total duration of session 1 in seconds)
% 4) Shift all session-2 onsets by +t_concat.
% 5) Merge session-1 and shifted session-2 onsets/durations into a single
%    structure:
%    - Most conditions are concatenated horizontally ([a, b]).
%    - One special-case condition (index i == 5) is concatenated vertically
%      ([a; b]) to preserve its expected shape/format in SPM.
% 6) Save a new combined condition file per subject in `outdir` with a
%    standardized name (e.g., sub-4001p4_task-all.mat).
%
% Inputs:
% - Session-specific multi-condition files:
%     root/*<subjectID>*1.mat
%     root/*<subjectID>*2.mat
% - TR and number of volumes per session (n_Series) to compute run duration
%
% Outputs:
% - One concatenated multi-condition file per subject in `outdir`:
%     sub-<ID>p4_task-all.mat
%   containing: names, onsets, durations
%
% Notes:
% - This assumes both sessions have the same duration: n_Series * TR.
% - Ensure the ordering and number of conditions match between session 1 and 2.
% - The special-case condition at index 5 likely corresponds to a condition
%   stored as a column vector/cell array that must be vertically concatenated.
% - Update `root`, `f_root`, and `outdir` to match your environment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all

root = '/results/fMRI_analysis/condition_matfiles_simple/';
f_root = '/results/fMRIprep/s4fmap_gr1';
outdir = '/results/fMRI_analysis/concat_conditions/first_level_complex';

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
    durations = c_concat_mat.durations;
    names = c_concat_mat.names;

    ses_name = ['sub-' subj_id{iSub} 'p4_task-all'];
    
    save([outdir '/' ses_name '.mat'], ...
            'names','onsets','durations')
end


aa= load('/data/p_02809/MeMoSLAP/fMRI_analysis/derivatives/first_level_fmap/concat_conditions/first_level_complex/sub-4001p4_task-all.mat')
