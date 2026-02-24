%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: First-level (subject-level) GLM analysis in SPM12 (n-back task)
%
% Author: Gerasimos Gerardos
% Project: MeMoSLAP Project (P04–Leipzig)
% Paper: Task load modulates network interactions between bilateral fronto-
%        parietal and cerebellar areas during verbal working memory
%
% Description:
% This script runs a first-level (subject-level) fMRI GLM analysis in SPM12
% for an n-back working memory paradigm (sessions: nback1, nback2). It is
% designed for BIDS/fMRIPrep-derived derivatives, but can be adapted to other
% folder structures.
%
% For each subject, the script:
% 1) Model specification (spm.stats.fmri_spec)
%    - Loads smoothed preprocessed BOLD time series (prefix defined by
%      func_prefix, e.g., 'smooth*' in MNI space).
%    - Sets timing parameters (TR, microtime resolution and onset).
%    - Loads task condition regressors from a precomputed SPM multi-condition
%      .mat file (names, onsets, durations), created separately (e.g., by
%      create_condition_matfiles.m).
%    - Adds nuisance regressors either as:
%        a) standard SPM motion regressors (6 realignment params), or
%        b) a custom nuisance matrix R stored in a .mat file (e.g., 24 motion
%           parameters + FD censoring + aCompCor + block order + session reg).
%    - Applies a high-pass filter (HPF = 128s), AR(1) serial correlations,
%      and an explicit brain mask from fMRIPrep.
%
% 2) Model estimation (spm.stats.fmri_est)
%    - Estimates the specified GLM (optional via run_estimation flag).
%
% 3) Contrast specification (spm.stats.con)
%    - Defines multiple t-contrasts for each load level (0/1/2/3-back) and
%      pairwise comparisons between loads (e.g., 3back > 0back), as well as
%      combined contrasts (e.g., 123back > 0back).
%    - Defines F-contrasts testing effects across multiple load conditions.
%    - Optionally deletes previously defined contrasts before saving new ones.
%
% 4) Optional results reporting (spm.stats.results)
%    - Generates a results report for all contrasts using FWE correction
%      (p < 0.05) and a minimum cluster extent threshold.
%    - Uses a subject-specific brain mask for masking in results.
%    - Runs in command-line mode with figures suppressed for batch execution.
%
% Inputs:
% - Smoothed functional images from fMRIPrep derivatives (MNI space NIfTI)
% - Multi-condition condition .mat files (names/onsets/durations) per session
% - Nuisance regressor files (either SPM rp*.txt or custom R .mat files)
% - Brain mask from fMRIPrep (space-MNI152NLin2009cAsym_desc-brain_mask.nii)
%
% Outputs:
% - First-level SPM model directory per subject (SPM.mat, beta images, etc.)
% - Contrast images (con_*.nii) and statistical maps (spmT_*.nii / spmF_*.nii)
% - (Optional) results report objects generated via SPM batch
%
% Notes:
% - Update paths (spm_path, root_dir, mask_root, condition_matfiles_dir,
%   nuisanceRegs_dir, outdir) before running.
% - Ensure the condition .mat and nuisance .mat files exist per subject/task.
% - This script currently combines sessions into a single output directory
%   (combine_ses) while using sess(iTask) for session-wise specification.
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

% Directory where condition .mat files are stored
% (can be created using create_condition_matfiles.m)
condition_matfiles_dir = '/results/fMRI_analysis/condition_matfiles_simple';

% Standard nuisance regressors are the 6 realignment parameters. 
% Set nuisanceRegs_SPM12_motionRegressors = 1 if you want to use only them.
% If you need more/other nuisance regressors, create a .mat file containing 
% an 'R' matrix of nuisance regressors as columns, save it in
% nuisanceRegs_dir, and set nuisanceRegs_SPM12_motionRegressors = 0;
nuisanceRegs_SPM12_motionRegressors = 0;
if nuisanceRegs_SPM12_motionRegressors == 0
    nuisanceRegs_dir = '/results/fMRI_analysis/nuisance_regressors/motionRegs24_FD09_aCompCor10_order_withsess';
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
        
        curr_output_dir = [outdir '/' curr_subj_name '/combine_ses']
        mkdir(curr_output_dir)
        
        co = dir(curr_output_dir);

        if any(cell2mat({co.isdir})==0)
            ierror = ['Error: Directory already processed! Skipping subject ' curr_subj_name]
            errors = [errors; ierror];
            continue;
        end
        %%
        clear matlabbatch

        for iTask = 1:numel(tasks)
            
            curr_task = tasks{iTask};

            %curr_output_dir = [curr_output_dir '/ses-' curr_task '/func'];
            %mkdir(curr_output_dir)
        
            
            
            %% fMRI model specification
            matlabbatch{1}.spm.stats.fmri_spec.dir = {curr_output_dir};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = microtime_resolution;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = microtime_onset;
            
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
            
            func_struct = {};
            for iVol = 1:numel(func_vols)
                func_struct(iVol) = {[curr_func_image_path ',' num2str(iVol)]};
            end
            func_struct = func_struct'
            


            %% get multiple conditions matfile
            cond_matfile = dir([condition_matfiles_dir '/*' cell2mat(regexp(curr_subj_name,'\d*','Match')) '*' num2str(iTask) '.mat']);
            
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

                nuisance_file = dir([nuisanceRegs_dir '/*' cell2mat(regexp(curr_subj_name,'\d*','Match')) '*' num2str(iTask) '.mat']);
            end
            
            if numel(nuisance_file) ~= 1
                ierror = ['Error: Not exactly 1 nuisance file found! Skipping session ' curr_subj_name '_' curr_task]
                errors = [errors; ierror];
                continue;
            end
            
            nuisance_file_path = [nuisance_file.folder '/' nuisance_file.name];
            
            %% Define matlabbatch for current run
            matlabbatch{1}.spm.stats.fmri_spec.sess(iTask).scans = func_struct;
            matlabbatch{1}.spm.stats.fmri_spec.sess(iTask).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(iTask).multi = {cond_matfile_path};
            matlabbatch{1}.spm.stats.fmri_spec.sess(iTask).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(iTask).multi_reg = {nuisance_file_path};
            matlabbatch{1}.spm.stats.fmri_spec.sess(iTask).hpf = 128;
            
        end
        
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
        
        %%% Model estimation
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
        if run_estimation == 1
            spm_jobman('run',matlabbatch) % Run the batch
        end
    
        %% Contrast manager
    
        % get SPM.mat file
        SPM_mat = dir([curr_output_dir '/SPM.mat']);
    
        if numel(SPM_mat) ~= 1
            ierror = ['Error: Not exactly 1 SPM.mat file found! Skipping participant '  curr_subj_name]
            errors = [errors; ierror];
            continue;
        end
    
        SPM_mat_path = [SPM_mat.folder '/' SPM_mat.name];
    
        clear matlabbatch
    
        matlabbatch{1}.spm.stats.con.spmmat(1) = {SPM_mat_path};
        
        %%% Define your contrasts here
        con_counter = 1;
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '0back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [1 0 0 0 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'repl';
        con_counter = con_counter + 1;
    
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '1back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [0 1 0 0 0 ];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'repl';
        con_counter = con_counter + 1;
    
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '2back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [0 0 1 0 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'repl';
        con_counter = con_counter + 1;

        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '3back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [0 0 0 1 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'repl';
        con_counter = con_counter + 1;

        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '1back > 0back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [-1 1 0 0 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'repl';
        con_counter = con_counter + 1;

        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '2back > 0back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [-1 0 1 0 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'repl';
        con_counter = con_counter + 1;
        
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '3back > 0back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [-1 0 0 1 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'repl';
        con_counter = con_counter + 1;

        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '2back > 1back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [0 -1 1 0 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'repl';
        con_counter = con_counter + 1;

        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '3back > 1back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [0 -1 0 1 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'repl';
        con_counter = con_counter + 1;
        
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '3back > 2back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [0 0 -1 1 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'repl';
        con_counter = con_counter + 1;

        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '0back > 1back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [1 -1 0 0 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'repl';
        con_counter = con_counter + 1;

        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '0back > 2back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [1 0 -1 0 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'repl';
        con_counter = con_counter + 1;
        
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '0back > 3back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [1 0 0 -1 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'repl';
        con_counter = con_counter + 1;

        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '1back > 2back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [0 1 -1 0 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'repl';
        con_counter = con_counter + 1;

        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '1back > 3back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [0 1 0 -1 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'repl';
        con_counter = con_counter + 1;
        
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '2back > 3back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [0 0 1 -1 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'repl';
        con_counter = con_counter + 1;

        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '123back > 0back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [-1 1/3 1/3 1/3 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'repl';
        con_counter = con_counter + 1;
        
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '23back > 01back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [-1/2 -1/2 1/2 1/2 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'repl';
        con_counter = con_counter + 1;

        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '0back > 123back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [1 -1/3 -1/3 -1/3 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'repl';
        con_counter = con_counter + 1;
        
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '01back > 23back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [1/2 1/2 -1/2 -1/2 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'repl';
        con_counter = con_counter + 1;

        % Define the F contrast matrix
        F_contrast_matrix = [1 0 0 0 ;
                             0 1 0 0 ;
                             0 0 1 0 ;
                             0 0 0 1 ;];
        
        % Add the F contrast to matlabbatch
        matlabbatch{1}.spm.stats.con.consess{con_counter}.fcon.name = 'F 0123back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.fcon.weights = F_contrast_matrix;
        matlabbatch{1}.spm.stats.con.consess{con_counter}.fcon.sessrep = 'repl';
        con_counter = con_counter + 1;
        
        % Define the F contrast matrix
        F_contrast_matrix = [0 0 0 0 ;
                             0 1 0 0 ;
                             0 0 1 0 ;
                             0 0 0 1 ;];
        
        % Add the F contrast to matlabbatch
        matlabbatch{1}.spm.stats.con.consess{con_counter}.fcon.name = 'F 123back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.fcon.weights = F_contrast_matrix;
        matlabbatch{1}.spm.stats.con.consess{con_counter}.fcon.sessrep = 'repl';

        matlabbatch{1}.spm.stats.con.delete = 1; % delete previously defined contrasts (if there are any)
        
        %%% Results report
        if run_results_report == 1

             %%% get fMRI data
            curr_anat_dir = ['/results/fMRI_analysis/fMRIprep/fmap_gr1/' curr_subj.name '/anat'];
            mask_images = dir([curr_anat_dir '/' '*space-MNI152NLin2009cAsym_desc-brain_mask.nii']);
            if isempty(mask_images)
                zipped_mask = dir([curr_anat_dir '/' '*space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz']);
                if isempty(zipped_mask)
                    ierror = ['Error: Not exactly 1 zipped.gz nifti mask found! Skipping session ' curr_subj_name];
                    errors = [errors; ierror];
                    continue;
                end
                mgzipfilenames = [zipped_mask.folder '/' zipped_mask.name];
                gunzip(mgzipfilenames, curr_anat_dir)
                mask_images = dir([curr_anat_dir '/' '*space-MNI152NLin2009cAsym_desc-brain_mask.nii']);
            end
        
            matlabbatch{2}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{2}.spm.stats.results.conspec.titlestr = '';
            matlabbatch{2}.spm.stats.results.conspec.contrasts = Inf; % Inf = show all contrasts
            matlabbatch{2}.spm.stats.results.conspec.threshdesc = 'FWE';
            matlabbatch{2}.spm.stats.results.conspec.thresh = 0.05;
            matlabbatch{2}.spm.stats.results.conspec.extent = 20;
            matlabbatch{2}.spm.stats.results.conspec.conjunction = 1;
            matlabbatch{2}.spm.stats.results.conspec.mask.image.name = {[mask_images.folder '/' mask_images.name ',1']};
            matlabbatch{2}.spm.stats.results.conspec.mask.image.mtype = '';
            matlabbatch{2}.spm.stats.results.units = 1;
            matlabbatch{2}.spm.stats.results.export = cell(1, 0);
    
        end
    
        spm_jobman('run',matlabbatch) % Run the batch
        
    end
end
