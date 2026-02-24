%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: First-level DCM-GLM specification and estimation with concatenated runs
%
% Author: Gerasimos Gerardos
% Project: MeMoSLAP Project (P04–Leipzig)
% Paper: Task load modulates network interactions between bilateral fronto-
%        parietal and cerebellar areas during verbal working memory
%
% Description:
% This script specifies and estimates a subject-level first-level GLM in SPM12
% that is tailored for the DCM analysis pipeline. It models the fMRI data as a
% *single concatenated run* by stacking two sessions (nback1 followed by nback2)
% into one continuous scan list, and pairing that concatenated time series with:
% - concatenated multi-condition timing files (names/onsets/durations; *all.mat)
% - concatenated nuisance regressors (R matrix; *all.mat)
%
% The script then calls `spm_fmri_concatenate` so SPM correctly handles the
% two-run structure internally (e.g., session boundary / filtering / whitening)
% while the design is specified as one session for downstream DCM use.
%
% Workflow (per subject):
% 1) Collect smoothed functional images for nback1 and nback2 (fMRIPrep outputs),
%    expand them into volume-wise scan entries, and stack them into `func_struct`
%    (assumes 1420 volumes per session; total 2840 volumes).
% 2) Load the subject’s concatenated DCM condition file (*all.mat) containing
%    names, onsets, and durations on a single continuous time axis.
% 3) Load the subject’s concatenated nuisance regressor file (*all.mat) containing
%    R (e.g., motion, FD spikes, aCompCor, block order, behavioral outlier blocks).
% 4) Specify the first-level GLM in SPM12 using canonical HRF (configurable),
%    high-pass filtering, AR(1) autocorrelation model, and an explicit brain mask.
% 5) Run model specification, then apply `spm_fmri_concatenate(SPM.mat,[1420,1420])`.
% 6) Estimate the model (Classical / OLS), producing an SPM.mat ready for DCM steps.
%
% Inputs:
% - Smoothed functional images (per subject, per session):
%     root_dir/sub-*/ses-nback*/func/smooth*MNI..._bold.nii(.gz)
% - Concatenated DCM condition files:
%     condition_matfiles_dir/*<subID>*all.mat
% - Concatenated nuisance regressor files:
%     nuisanceRegs_dir/*<subID>*all.mat
% - Subject brain mask (fMRIPrep):
%     mask_root/sub-*/anat/*space-MNI..._desc-brain_mask.nii(.gz)
%
% Outputs:
% - Subject-specific DCM-GLM directory:
%     outdir/sub-*/concat_for_DCM/
%   containing SPM.mat and estimated GLM outputs (betas, residual variance, etc.).
%
% Notes:
% - Assumes exactly 1420 volumes in each session; update if different.
% - The concatenated condition and nuisance files must match the scan order
%   (nback1 first, then nback2) and total length (2840 volumes).
% - This script stops after model estimation setup; contrasts can be added
%   later if needed, but are not required for DCM model estimation itself.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Clear the workspace
clear all
close all

%% SPM setup
spm_path = '/MATLAB/spm12'; % save SPM12 path for later
addpath(spm_path)
spm fmri
spm_get_defaults('cmdline', true);
set(0, 'DefaultFigureVisible', 'off');

%% Setup
root_dir = '/results/fMRIprep/s4fmap_gr1/'; %output_gr1/';
mask_root = '/results/fMRIprep/fmap_gr1/';

subj_dirs = dir([root_dir '/sub-*']);
subj_dirs = subj_dirs(cellfun(@(x) ~endsWith(x,'04'), {subj_dirs.name}));

% Directory where condition .mat files are stored
% (can be created using create_condition_matfiles.m)
condition_matfiles_dir = '/results/fMRI_analysis/concat_conditions/first_level_DCM';

% Standard nuisance regressors are the 6 realignment parameters. 
% Set nuisanceRegs_SPM12_motionRegressors = 1 if you want to use only them.
% If you need more/other nuisance regressors, create a .mat file containing 
% an 'R' matrix of nuisance regressors as columns, save it in
% nuisanceRegs_dir, and set nuisanceRegs_SPM12_motionRegressors = 0;
nuisanceRegs_SPM12_motionRegressors = 0;
if nuisanceRegs_SPM12_motionRegressors == 0
    nuisanceRegs_dir = '/results/fMRI_analysis/nuisance_regressors/motionRegs24_FD09_aCompCor10_order_withoutsess_concat';
end

TR = 1; % repetition time
microtime_resolution = 72; % set to number of slices
microtime_onset = 36; % set to middle slice (or reference slice when slice timing was performed)

func_folder = 'func';
func_prefix = 'swau'; % prefix of functional images to use for analysis
tasks = {'nback1', 'nback2'};

% Define your HRF model:
% [0 0] = canonical HRF only; [1 0] = canonical HRF + temporal derivative; 
% [1 1] = canonical HRF + temporal and dispersion derivatives
HRF_model = [0 0]; 

run_estimation = 1; % run model estimation? (set to 0 if done already, e.g. if you only want to define new contrasts)
run_results_report = 1; % run results report? (set to 0 to not show a results plot yet -> speeds up running first-level analyses across a lot of subjects)

outdir = '/results/fMRI_analysis/first_level_simple'; %first_level_analysis';
if isempty(dir(outdir))
    mkdir(outdir)
end

errors = {};

%% For each subject, run their first-level analysis
for iSubject = 1:numel(subj_dirs)
    if subj_dirs(iSubject).isdir == 1 %&& any(cellfun(@(x) all(join([subj_dirs(iSubject).name ".html"], "",2)==x), {subj_dirs.name}))
        curr_subj = subj_dirs(iSubject);
        curr_subj_name = curr_subj.name
        curr_subj_dir = [curr_subj.folder '/' curr_subj.name];
        
        curr_output_dir = [outdir '/' curr_subj_name '/concat_for_DCM']
        % Check and create the directory
        if ~exist(curr_output_dir, 'dir')
            mkdir(curr_output_dir);
        end


        
        co = dir(curr_output_dir);

        if any(cell2mat({co.isdir})==0)
            ierror = ['Error: Directory already processed! Skipping subject ' curr_subj_name]
            errors = [errors; ierror];
            continue;
        end
        %%
        clear matlabbatch
        
        %% fMRI model specification
        matlabbatch{1}.spm.stats.fmri_spec.dir = {curr_output_dir};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = microtime_resolution;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = microtime_onset;
        
        func_struct = {};
        for iTask = 1:numel(tasks)
            
            curr_task = tasks{iTask};

            %curr_output_dir = [curr_output_dir '/ses-' curr_task '/func'];
            %mkdir(curr_output_dir)
        
            

            
            %%% get fMRI data
            curr_func_dir = [curr_subj_dir '/ses-' curr_task '/' func_folder];
            func_images = dir([curr_func_dir '/' 'smooth*MNI152NLin2009cAsym_desc-preproc_bold.nii'])
            if isempty(func_images)
                zipped_func = dir([curr_func_dir '/' 'smooth*MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'])
                if isempty(zipped_func)
                    ierror = ['Error: Not exactly 1 zipped.gz nifti found! Skipping session ' curr_subj_name '_' curr_task]
                    errors = [errors; ierror];
                    continue;
                end
                gzipfilenames = [zipped_func.folder '/' zipped_func.name];
                gunzip(gzipfilenames, curr_func_dir)
                func_images = dir([curr_func_dir '/' 'smooth*MNI152NLin2009cAsym_desc-preproc_bold.nii'])
            end

        
            %% For each run, specify the design
      
            
            %%% get functional volumes
            curr_func_image = func_images;
            
            if numel(curr_func_image) ~= 1
               ierror = ['Error: Not exactly 1 functional nifti found! Skipping session ' curr_subj_name '_' curr_task]
               errors = [errors; ierror];
               continue;
            end
            
            curr_func_image_path = [curr_func_image.folder '/' curr_func_image.name];
            
            % load nifti & get number of volumes
            func_vols = spm_vol(curr_func_image_path);
            
            
            for iVol = 1:numel(func_vols)
                func_struct(iVol+((iTask-1)*1420)) = {[curr_func_image_path ',' num2str(iVol)]};
            end
        end
        func_struct = func_struct'

            


        %% get multiple conditions matfile
        cond_matfile = dir([condition_matfiles_dir '/*' cell2mat(regexp(curr_subj_name,'\d*','Match')) '*all.mat']);
        
        if numel(cond_matfile) ~= 1
            ierror = ['Error: Not exactly 1 condition matfile found! Skipping session ' curr_subj_name '_' curr_task]
            errors = [errors; ierror];
            continue;
        end
        
        cond_matfile_path = [cond_matfile.folder '/' cond_matfile.name];
        
        %% get nuisance regressors file

        if nuisanceRegs_SPM12_motionRegressors == 1
            nuisance_file = dir([curr_func_dir '/rp*run-0' num2str(iRun) '_bold.txt']);
        else

            nuisance_file = dir([nuisanceRegs_dir '/*' cell2mat(regexp(curr_subj_name,'\d*','Match')) '*all.mat']);
        end
        
        if numel(nuisance_file) ~= 1
            ierror = ['Error: Not exactly 1 nuisance file found! Skipping session ' curr_subj_name '_' curr_task]
            errors = [errors; ierror];
            continue;
        end
        
        nuisance_file_path = [nuisance_file.folder '/' nuisance_file.name];
        
        %% Define matlabbatch for current run
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = func_struct;
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {cond_matfile_path};
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {nuisance_file_path};
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128;
        
    
        
        if isempty(func_images)
            if isempty(zipped_func)
                ['Confirming Error: zip: ' num2str(isempty(zipped_func))]
                continue;
                if numel(curr_func_image) ~= 1 
                   ['Confirming Error: func img ~=1: ' num2str(numel(curr_func_image))]
                   continue;
                    if numel(cond_matfile) ~= 1 
                        ['Confirming Error: cond files ~=1: ' num2str(numel(cond_matfile))]
                        continue;
                        if numel(nuisance_file) ~= 1
                            ['Confirming Error: no-regress files ~=1: ' num2str(numel(nuisance_file))]
                            continue;
                        end
                    end
                end   
            end
        end

        %%% get mask data
        curr_mask_dir = [mask_root curr_subj_name '/anat'];
        mask_images = dir([curr_mask_dir '/' '*space-MNI152NLin2009cAsym_desc-brain_mask.nii'])
        if isempty(mask_images)
            zipped_mask = dir([curr_mask_dir '/' '*MNI152NLin2009cAsym_desc-brain_mask.nii.gz'])
            if isempty(zipped_mask)
                ierror = ['Error: Not exactly 1 zipped.gz nifti mask found! Skipping session ' curr_subj_name '_' curr_task]
                errors = [errors; ierror];
                continue;
            end
            gzipfilenames = [zipped_mask.folder '/' zipped_mask.name];
            gunzip(gzipfilenames, curr_mask_dir)
            mask_images = dir([curr_mask_dir '/' '*MNI152NLin2009cAsym_desc-brain_mask.nii'])
        end

        mask_dir = fullfile(mask_images.folder, [mask_images.name ',1'])

    
        %% Define the rest of the batch
        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = HRF_model;
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.1;
        matlabbatch{1}.spm.stats.fmri_spec.mask = {mask_dir};
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
        
        spm_jobman('run',matlabbatch)
        
        % get SPM.mat file
        SPM_mat = dir([curr_output_dir '/SPM.mat']);
    
        if numel(SPM_mat) ~= 1
            ierror = ['Error: Not exactly 1 SPM.mat file found! Skipping participant '  curr_subj_name]
            errors = [errors; ierror];
            continue;
        end
    
        SPM_mat_path = [SPM_mat.folder '/' SPM_mat.name];

        spm_fmri_concatenate(SPM_mat_path, [1420,1420]); 
        
        clear matlabbatch


    
        %%% Model estimation
        matlabbatch{1}.spm.stats.fmri_est.spmmat(1) = {SPM_mat_path};
        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    
        if run_estimation == 1
            spm_jobman('run',matlabbatch) % Run the batch
        end
    
        %% Contrast manager
    
    
        clear matlabbatch
    
        
    end
end

