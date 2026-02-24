%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Second-level (group-level) one-sample t-tests in SPM12
%
% Author: Gerasimos Gerardos
% Project: MeMoSLAP Project (P04–Leipzig)
% Paper: Task load modulates network interactions between bilateral fronto-
%        parietal and cerebellar areas during verbal working memory
%
% Description:
% This script runs second-level (group-level) mass-univariate one-sample
% t-test analyses in SPM12 for a set of first-level contrast images.
% It assumes that subject-level GLMs have already been estimated and that
% corresponding contrast images (con_*.nii) exist for each subject.
%
% For each contrast specified in `con_images` / `con_names`, the script:
% 1) Collects the relevant first-level contrast image for every subject.
% 2) Specifies a one-sample t-test design (factorial_design.des.t1) to test
%    whether the group mean contrast differs from zero across subjects.
% 3) Estimates the second-level model (Classical / OLS estimation).
% 4) Defines two group-level contrasts:
%    - Positive:  +1  (A > B direction for the given first-level contrast)
%    - Negative:  -1  (B > A direction, i.e., sign-flipped)
% 5) Generates an initial results report thresholded at p < 0.001
%    uncorrected (extent threshold = 10 voxels) for quick inspection.
%
% Inputs:
% - First-level SPM directories per subject (e.g., first_level_dir/sub-*/...)
% - Subject-level contrast images in each subject’s model folder
%   (e.g., con_0005.nii, con_0006.nii, etc.)
%
% Outputs (per contrast):
% - A second-level results folder: outdir/<CONTRAST_NAME>/
% - Second-level SPM.mat, beta images, ResMS, etc.
% - Group-level contrast images (con_*.nii) for positive/negative effects
% - Thresholded results output from the SPM batch “Results” step
%
% Options / Notes:
% - You may optionally set `explicit_mask_path` to restrict inference
%   (e.g., to gray matter or a group brain mask). Leave empty ('') to skip.
% - The p < 0.001 threshold here is intended for first inspection only.
%   Final inference should be performed manually in SPM with multiple-
%   comparison correction (FWE/FDR at voxel or cluster level).
% - Update `spm_path`, `first_level_dir`, and `outdir` to match your setup.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Clear the workspace
clear all
close all

%% Add SPM
spm_path = '/MATLAB/spm12';
addpath(spm_path)
spm fmri
spm_get_defaults('cmdline', true);
set(0, 'DefaultFigureVisible', 'off');

%% Setup
first_level_dir = '/results/fMRI_analysis/first_level_simple/';

subj_dirs = dir([first_level_dir '/sub-*']);
subj_dirs = subj_dirs([1:3,5:end]);

% Define the contrast images & their corresponding names for which you
% want to run a one-sample t-test across subjects
con_images = {'con_0005.nii','con_0006.nii','con_0007.nii','con_0008.nii','con_0009.nii','con_0010.nii','con_0017.nii','con_0018.nii'};         
con_names = {'1BACK_VS_0BACK','2BACK_VS_0BACK','3BACK_VS_0BACK','2BACK_VS_1BACK','3BACK_VS_1BACK', '3BACK_VS_2BACK','123BACK_VS_0BACK','23BACK_VS_01BACK'};
         
outdir = '/results/fMRI_analysis/second_level_1T_simple';
if isempty(dir(outdir))
    mkdir(outdir)
end

% optional: add path to explicit mask (e.g. to restrict statistical tests to gray matter)
explicit_mask_path = ''; 
         
%% Run a second-level one-sample t-test for each contrast
for iCon = 1:numel(con_images)         

    curr_con_image = con_images{iCon}
    curr_con_name = con_names{iCon}
    
    curr_outdir = [outdir '/' curr_con_name]
    mkdir(curr_outdir)
    
    %% For the current contrast, get the con-images for each subject
    con_struct = {};
    for iSubject = 1:numel(subj_dirs)
        
        curr_subj = subj_dirs(iSubject);
        curr_subj_path = [curr_subj.folder '/' curr_subj.name '/combine_ses'];
        
        % get con-image for current subject
        subj_con_image = dir([curr_subj_path '/' curr_con_image]);
        
        if numel(subj_con_image) ~= 1
            error('Could not find exactly 1 con-image for current subject!')
        end
        
        subj_con_image_path = [subj_con_image.folder '/' subj_con_image.name];
        
        con_struct(iSubject) = {subj_con_image_path};
        
    end
    con_struct = con_struct'
    
    %% Specify and run the second-level batch
    clear matlabbatch
    
    % Specify second-level design
    matlabbatch{1}.spm.stats.factorial_design.dir = {curr_outdir};
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = con_struct;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {explicit_mask_path};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    % Estimate
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    %%% Contrasts
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    
    % Positive activations (A > B)
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = [curr_con_name '_positive'];
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
    % Negative activations (B > A)
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = [curr_con_name '_negative'];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = -1;
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.delete = 1;

    % Results report
    matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
    matlabbatch{4}.spm.stats.results.conspec.contrasts = Inf;
    matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'none';
    matlabbatch{4}.spm.stats.results.conspec.thresh = 0.001;
    matlabbatch{4}.spm.stats.results.conspec.extent = 10;
    matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
    matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
    matlabbatch{4}.spm.stats.results.units = 1;
    matlabbatch{4}.spm.stats.results.export = cell(1, 0);
    
    %%% Run the batch
    spm_jobman('run',matlabbatch)
    
end
