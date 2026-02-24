%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Create SPM multi-condition files for the DCM GLM model
%
% Author: Gerasimos Gerardos
% Project: MeMoSLAP Project (P04–Leipzig)
% Paper: Task load modulates network interactions between bilateral fronto-
%        parietal and cerebellar areas during verbal working memory
%
% Description:
% This script converts block timing tables (onsets and durations) into SPM
% “multiple conditions” .mat files (names, onsets, durations) that are used
% to specify the first-level GLM for the DCM analysis pipeline.
%
% For each subject/session, it:
% - Loads block start times and block durations from CSV timing files.
% - Removes fixation blocks ("fix") and instruction blocks ("instr") so the
%   model contains task blocks only.
% - Keeps the separate n-back task conditions (e.g., 0/1/2/3-back).
% - Adds an additional pooled task regressor ("all_conditions") by duplicating
%   all "*back" blocks and relabeling them as one combined condition.
% - Extracts, sorts, and pairs onsets and durations per condition, and saves
%   the result as an SPM-compatible .mat file.
%
% Output:
% - One .mat file per subject/session in `outdir`, containing:
%     names, onsets, durations
%   (SPM multi-condition format for the DCM GLM model).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear all
close all

%%

outdir = 'results/fMRI_analysis/condition_matfiles_glm_for_dcm/';

onset_dir = '/results/fMRI_analysis/timings_SPM/block_starts/'; %change to your directory 
dur_dir = '/results/fMRI_analysis/timings_SPM/block_durations/'; %change to your directory 

start_dirs = dir([onset_dir '/sub-*']);
last_dirs = dir([dur_dir '/sub-*']);


for iSes = 1:numel(start_dirs)
    curr_ses = start_dirs(iSes).name(1:17);
    starts = readtable([onset_dir start_dirs(iSes).name], ReadVariableNames=true);
    lasts = readtable([dur_dir last_dirs(iSes).name], ReadVariableNames=true);
    conds = starts.block;

    %remove fixation 
    ifix = cell2mat(cellfun(@(x)  contains( x , 'fix' ) ,conds,'UniformOutput',false));
    starts = starts(~ifix,:);
    lasts = lasts(~ifix,:);
    conds = conds(~ifix,:);

    %remove instruction 
    iinstr = cell2mat(cellfun(@(x)  contains( x , 'instr' ) ,conds,'UniformOutput',false));
    starts = starts(~iinstr ,:);
    lasts = lasts(~iinstr ,:);
    conds = conds(~iinstr ,:);

    %unify conditions
    cstarts = starts;
    clasts = lasts;
    cconds = conds;
    iback = cell2mat(cellfun(@(x)  contains( x , 'back' ) ,cconds,'UniformOutput',false));
    cstarts(iback,'block') = {'all_conditions'};
    clasts(iback,'block') = {'all_conditions'};
    cconds(iback) = {'all_conditions'};

    starts = [starts;cstarts];
    lasts = [lasts;clasts];
    conds = [conds, cconds];

    uconds = unique(conds);
    
    names = {};
    onsets = {};
    durations = {};

    for iuCond = 1:numel(uconds)
        
        curr_cond = uconds{iuCond};
        iCond = cell2mat(cellfun(@(x)  contains( x , curr_cond ) ,conds,'UniformOutput',false));
    
        curr_onsets = table2array(starts(iCond,1:6));
        curr_onsets = curr_onsets(~isnan(curr_onsets));
        [curr_onsets, on_idx] = sort(curr_onsets);


        curr_durations = table2array(lasts(iCond, 1:6));
        curr_durations = curr_durations(~isnan(curr_durations));
        curr_durations = curr_durations(on_idx);

        names(iuCond) = {curr_cond};
        onsets(iuCond) = {curr_onsets};
        durations(iuCond) = {curr_durations};
    end

    

    ses_name = erase(curr_ses, ["_", "-"]);
    ses_name = [ses_name(1:3), '-', ses_name(4:9), '_', ses_name(10:end-1), '-', ses_name(end)];
    save([outdir '/' ses_name '.mat'], ...
            'names','onsets','durations')
end


