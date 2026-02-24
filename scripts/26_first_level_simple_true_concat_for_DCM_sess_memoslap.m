%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: First-level GLM for DCM (concatenated sessions + session modulators)
%
% Author: Gerasimos Gerardos
% Project: MeMoSLAP Project (P04–Leipzig)
% Paper: Task load modulates network interactions between bilateral fronto-
%        parietal and cerebellar areas during verbal working memory
%
% Description:
% This script specifies and estimates a subject-level SPM12 GLM that is
% tailored for subsequent DCM analysis using a *concatenated* fMRI time series.
% Two sessions (nback1 and nback2) are concatenated into a single continuous
% scan list (e.g., 1420 + 1420 volumes). The corresponding condition file
% (multi-condition .mat) is assumed to already contain onsets that are aligned
% to the concatenated timeline, and it additionally contains explicit
% "session_1" and "session_2" regressors intended to be used as modulatory
% inputs in the DCM (B-matrix) to model session-specific connectivity changes.
%
% Workflow:
% 1) Loop over subjects (excluding IDs ending in '04').
% 2) Build a concatenated scan list (func_struct) across nback1 and nback2.
% 3) Load the concatenated multi-condition file (*all.mat) that includes:
%       - task regressors (e.g., 0/1/2/3-back, all_conditions)
%       - session_1 and session_2 regressors (for DCM modulation by session)
% 4) Load concatenated nuisance regressors (multi_reg; same temporal length as
%    the concatenated scans).
% 5) Specify the model (single SPM session with concatenated scans).
% 6) Run SPM model specification to create SPM.mat.
% 7) Call spm_fmri_concatenate(SPM.mat, [1420, 1420]) so SPM treats the design
%    correctly as concatenated sessions.
% 8) Estimate the model (spm.stats.fmri_est).
%
% Inputs:
% - Preprocessed functional images per session:
%     .../ses-nback1/func/smooth*MNI...bold.nii(.gz)
%     .../ses-nback2/func/smooth*MNI...bold.nii(.gz)
% - Concatenated condition .mat file with session modulators:
%     condition_matfiles_dir/*<subid>*all.mat
% - Concatenated nuisance regressor .mat file (variable 'R'):
%     nuisanceRegs_dir/*<subid>*all.mat
% - Brain mask from fmriprep:
%     mask_root/sub-*/anat/*brain_mask.nii(.gz)
%
% Output:
% - Subject-level SPM GLM directory:
%     outdir/sub-*/concat_for_DCM_sess_mod/
%   containing SPM.mat, beta images, residuals (optional), etc.
%
% Notes:
% - This script does NOT define contrasts (left empty on purpose for DCM).
% - Ensure that the concatenated condition and nuisance files match the
%   concatenated scan length exactly (e.g., 2840 rows for 2*1420 volumes).
% - If you change per-session length (1420), update both:
%     - scan concatenation indexing (iVol+((iTask-1)*1420))
%     - spm_fmri_concatenate call [1420,1420]
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
condition_matfiles_dir = '/results/fMRI_analysis/concat_conditions/first_level_DCM_sess_mod';

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
        
        curr_output_dir = [outdir '/' curr_subj_name '/concat_for_DCM_sess_mod']
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

