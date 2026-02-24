%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: First-level GLM in SPM12 using concatenated sessions (nback1+nback2)
%         for VOI extraction / DCM-ready time series
%
% Author: Gerasimos Gerardos
% Project: MeMoSLAP Project (P04–Leipzig)
% Paper: Task load modulates network interactions between bilateral fronto-
%        parietal and cerebellar areas during verbal working memory
%
% Purpose (DCM / VOI extraction):
% This script specifies and estimates a first-level SPM GLM using a single
% *concatenated* time series constructed by stacking two sessions (nback1 and
% nback2) into one continuous run. It uses concatenated condition timing
% files and concatenated nuisance regressors so that the estimated residuals
% and adjusted data are aligned with the concatenated series—this is
% especially useful for VOI extraction prior to DCM.
%
% Description:
% For each subject, the script:
% 1) Builds one scan list (`func_struct`) that contains all volumes from
%    nback1 followed by all volumes from nback2 (assumes 1420 volumes each).
% 2) Loads a *concatenated* multi-condition .mat file (names, onsets,
%    durations) where session-2 onsets have been shifted by the duration of
%    session 1 (created beforehand in `condition_matfiles_dir`).
% 3) Loads a *concatenated* nuisance regressor matrix R (motion, FD spikes,
%    aCompCor, block-order, and optional d′ block outliers), saved in
%    `nuisanceRegs_dir` (created beforehand).
% 4) Specifies the GLM in SPM12 (timing parameters, HRF model, HPF, AR(1),
%    and explicit fMRIPrep brain mask).
% 5) Runs `spm_fmri_concatenate` on the created SPM.mat with run lengths
%    [1420, 1420] so SPM treats the dataset as concatenated sessions while
%    handling session boundaries appropriately.
% 6) Estimates the model and defines a set of t- and F-contrasts for load
%    conditions (0/1/2/3-back) and load comparisons (e.g., 3>0, 23>01).
% 7) Optionally generates an automated results report (FWE p<0.05) using an
%    explicit subject brain mask.
%
% Inputs:
% - Smoothed functional images from fMRIPrep (per session):
%     root_dir/sub-*/ses-nback*/func/smooth*MNI..._bold.nii(.gz)
% - Concatenated multi-condition timing files:
%     condition_matfiles_dir/*<subID>*all.mat
% - Concatenated nuisance regressor files (variable R):
%     nuisanceRegs_dir/*<subID>*all.mat
% - Subject-specific fMRIPrep brain mask:
%     mask_root/sub-*/anat/*space-MNI..._desc-brain_mask.nii(.gz)
%
% Outputs:
% - First-level SPM model folder per subject:
%     outdir/sub-*/concat_ses/
%   including SPM.mat, beta images, residuals (if enabled), contrast images,
%   and statistical maps (spmT/spmF).
%
% Assumptions / Notes:
% - Assumes exactly 1420 volumes per session and that the concatenation order
%   is nback1 then nback2.
% - The “concat” condition and nuisance files must match this exact ordering
%   and onset-shifting scheme.
% - Set `run_results_report` to 0 to speed batch processing.
% - Update all directory paths (SPM path, root_dir, condition/nuisance dirs,
%   output dir) to match your environment.
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
condition_matfiles_dir = '/results/fMRI_analysis/concat_conditions/first_level_complex';

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
        
        curr_output_dir = [outdir '/' curr_subj_name '/concat_ses']
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
    
        matlabbatch{1}.spm.stats.con.spmmat(1) = {SPM_mat_path};
        
        %%% Define your contrasts here
        con_counter = 1;
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '0back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [1 0 0 0 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'none';
        con_counter = con_counter + 1;
    
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '1back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [0 1 0 0 0 ];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'none';
        con_counter = con_counter + 1;
    
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '2back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [0 0 1 0 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'none';
        con_counter = con_counter + 1;

        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '3back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [0 0 0 1 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'none';
        con_counter = con_counter + 1;

        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '1back > 0back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [-1 1 0 0 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'none';
        con_counter = con_counter + 1;

        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '2back > 0back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [-1 0 1 0 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'none';
        con_counter = con_counter + 1;
        
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '3back > 0back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [-1 0 0 1 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'none';
        con_counter = con_counter + 1;

        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '2back > 1back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [0 -1 1 0 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'none';
        con_counter = con_counter + 1;

        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '3back > 1back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [0 -1 0 1 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'none';
        con_counter = con_counter + 1;
        
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '3back > 2back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [0 0 -1 1 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'none';
        con_counter = con_counter + 1;

        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '0back > 1back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [1 -1 0 0 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'none';
        con_counter = con_counter + 1;

        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '0back > 2back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [1 0 -1 0 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'none';
        con_counter = con_counter + 1;
        
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '0back > 3back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [1 0 0 -1 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'none';
        con_counter = con_counter + 1;

        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '1back > 2back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [0 1 -1 0 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'none';
        con_counter = con_counter + 1;

        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '1back > 3back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [0 1 0 -1 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'none';
        con_counter = con_counter + 1;
        
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '2back > 3back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [0 0 1 -1 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'none';
        con_counter = con_counter + 1;

        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '123back > 0back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [-1 1/3 1/3 1/3 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'none';
        con_counter = con_counter + 1;
        
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '23back > 01back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [-1/2 -1/2 1/2 1/2 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'none';
        con_counter = con_counter + 1;

        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '0back > 123back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [1 -1/3 -1/3 -1/3 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'none';
        con_counter = con_counter + 1;
        
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.name = '01back > 23back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.weights = [1/2 1/2 -1/2 -1/2 0];
        matlabbatch{1}.spm.stats.con.consess{con_counter}.tcon.sessrep = 'none';
        con_counter = con_counter + 1;

        % Define the F contrast matrix
        F_contrast_matrix = [1 0 0 0 ;
                             0 1 0 0 ;
                             0 0 1 0 ;
                             0 0 0 1 ;];
        
        % Add the F contrast to matlabbatch
        matlabbatch{1}.spm.stats.con.consess{con_counter}.fcon.name = 'F 0123back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.fcon.weights = F_contrast_matrix;
        matlabbatch{1}.spm.stats.con.consess{con_counter}.fcon.sessrep = 'none';
        con_counter = con_counter + 1;
        
        % Define the F contrast matrix
        F_contrast_matrix = [0 0 0 0 ;
                             0 1 0 0 ;
                             0 0 1 0 ;
                             0 0 0 1 ;];
        
        % Add the F contrast to matlabbatch
        matlabbatch{1}.spm.stats.con.consess{con_counter}.fcon.name = 'F 123back';
        matlabbatch{1}.spm.stats.con.consess{con_counter}.fcon.weights = F_contrast_matrix;
        matlabbatch{1}.spm.stats.con.consess{con_counter}.fcon.sessrep = 'none';

        matlabbatch{1}.spm.stats.con.delete = 1; % delete previously defined contrasts (if there are any)
        
        %%% Results report
        if run_results_report == 1

             %%% get fMRI data
            curr_anat_dir = ['/data/p_02809/MeMoSLAP/fMRI_analysis/fMRIprep/fmap_gr1/' curr_subj.name '/anat'];
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
