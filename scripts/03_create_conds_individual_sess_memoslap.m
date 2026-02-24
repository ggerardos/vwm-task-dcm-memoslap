%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Create SPM condition .mat files from block timing tables
%
% Author: Gerasimos Gerardos
% Project: MeMoSLAP Project (P04–Leipzig)
% Paper: Task load modulates network interactions between bilateral fronto-
%        parietal and cerebellar areas during verbal working memory
%
% Description:
% This script converts participant-level block onset and duration CSV tables
% (previously generated from behavioral logs) into SPM-compatible condition
% .mat files containing the variables: names, onsets, and durations.
%
% Workflow:
% - Reads block start times and durations for each subject/session.
% - Removes fixation blocks (any condition containing "fix").
% - Unifies all instruction blocks (any condition containing "instr") into a
%   single condition named "instr".
% - For each remaining unique condition, extracts all onsets and matching
%   durations (sorted by onset time).
% - Saves one .mat file per subject/session into the specified output folder.
%
% Output:
% - One .mat file per subject/session in outdir, containing:
%     names      : cell array of condition names
%     onsets     : cell array of onset vectors (seconds)
%     durations  : cell array of duration vectors (seconds)
%
% Notes:
% - Update outdir, onset_dir, and dur_dir to match your local directory
%   structure before running.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
clear all
close all

%%

outdir = '/results/fMRI_analysis/condition_matfiles_simple/'; %change to your directory 

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

    %unify instruction
    iinstr = cell2mat(cellfun(@(x)  contains( x , 'instr' ) ,conds,'UniformOutput',false));
    starts(iinstr,'block') = {'instr'};
    lasts(iinstr,'block') = {'instr'};
    conds(iinstr) = {'instr'};
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


