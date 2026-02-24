%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Concatenate DCM-GLM condition files across sessions (nback1 + nback2)
%
% Author: Gerasimos Gerardos
% Project: MeMoSLAP Project (P04–Leipzig)
% Paper: Task load modulates network interactions between bilateral fronto-
%        parietal and cerebellar areas during verbal working memory
%
% Description:
% This script creates *concatenated* SPM multi-condition files (names, onsets,
% durations) for the DCM GLM model by merging the two session-specific
% condition .mat files (session 1 and session 2) into a single continuous
% timing file per subject.
%
% For each subject, it:
% - Locates the two DCM-GLM condition files (*1.mat and *2.mat).
% - Loads session-1 and session-2 variables (names, onsets, durations).
% - Shifts all session-2 onsets forward by the duration of session 1:
%     t_concat = n_Series * TR
%   so that session-2 events occur after session-1 events in a single
%   continuous time axis.
% - Concatenates onsets and durations across sessions for each condition.
%   One condition (index i == 5) is concatenated vertically to preserve the
%   expected shape/format of that regressor in SPM.
% - Saves a single “task-all” multi-condition .mat file per subject, to be
%   used for the concatenated DCM GLM model.
%
% Inputs:
% - Session-specific DCM-GLM condition files:
%     root/*<subID>*1.mat
%     root/*<subID>*2.mat
% - TR and number of volumes per session (n_Series) to compute run duration
%
% Outputs:
% - One concatenated DCM-GLM condition file per subject:
%     outdir/sub-<ID>p4_task-all.mat
%   containing variables: names, onsets, durations
%
% Notes:
% - Assumes both sessions have equal length (n_Series volumes, TR seconds).
% - Ensure session 1 and session 2 condition files have identical condition
%   ordering and structure.
% - Update `root`, `f_root`, and `outdir` to match your environment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all

root = '/results/fMRI_analysis/condition_matfiles_glm_for_dcm';
f_root = '/results/fMRIprep/s4fmap_gr1';
outdir = '/results/fMRI_analysis/concat_conditions/first_level_DCM';


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


aa= load(fullfile(outdir,'sub-4001p4_task-all.mat') )
