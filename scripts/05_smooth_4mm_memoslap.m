%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Spatial smoothing of fMRIPrep preprocessed BOLD images (SPM12)
%
% Author: Gerasimos Gerardos
% Project: MeMoSLAP Project (P04–Leipzig)
% Paper: Task load modulates network interactions between bilateral fronto-
%        parietal and cerebellar areas during verbal working memory
%
% Description:
% This script performs spatial smoothing of fMRIPrep-preprocessed functional
% BOLD images using SPM12. For each subject and task session, it locates the
% preprocessed BOLD NIfTI files in MNI space, ensures they are available in
% uncompressed .nii format (unzipping .nii.gz files when needed), copies the
% selected run to a subject/session-specific output directory, and applies
% SPM’s Gaussian smoothing with a user-defined kernel (FWHM).
%
% Workflow:
% - Adds SPM12 to the MATLAB path and initializes SPM for fMRI.
% - Iterates over all subjects in the fMRIPrep output directory (sub-*).
% - For each task session (base, rest1, rest2, nback1, nback2) and run:
%     1) Checks whether the run has already been smoothed in the output folder
%        to avoid duplicate processing (logs skips in an error list).
%     2) Locates the preprocessed BOLD image
%        (*MNI152NLin2009cAsym_desc-preproc_bold.nii or .nii.gz).
%     3) Unzips .nii.gz files if necessary and selects the correct run file.
%     4) Copies the NIfTI file into the output directory for that session/run.
%     5) Builds an SPM volume list (one entry per timepoint/volume).
%     6) Runs SPM12 smoothing (spm.spatial.smooth) with the specified FWHM.
%
% Parameters:
% - smooth_FWHM: smoothing kernel in mm (e.g., [4 4 4])
% - tasks: list of sessions to process
%
% Outputs:
% - Smoothed functional images saved in:
%     outdir/sub-*/ses-*/func/
%   with prefix: 'smooth' (e.g., smooth*.nii)
% - A log (cell array) of skipped runs/sessions stored in `errors`
%
% Notes:
% - Update `spm_path`, `root_dir`, and `outdir` to match your local setup.
% - The script assumes fMRIPrep naming conventions and MNI152NLin2009cAsym
%   output space.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Clear the workspace
clear all
close all

%% SPM setup
spm_path = '/MATLAB/spm12'; % save SPM12 path for later
addpath(spm_path)
spm fmri

%% Setup
root_dir = '/results/fMRIprep/fmap_gr1/';

subj_dirs = dir([root_dir '/sub-*']); %sub-*4021*

outdir ='/results/fMRIprep/s4fmap_gr1'; %first_level_analysis';

func_folder = 'func';
tasks = {'base','rest1','rest2','nback1','nback2'}; %'base','rest1','rest2','nback1', 

errors = {};

smooth_FWHM = [4 4 4];

for iSubject = 1:numel(subj_dirs)
    if subj_dirs(iSubject).isdir == 1 && any(contains({subj_dirs.name}, subj_dirs(iSubject).name(end-3:end)) & endsWith({subj_dirs.name}, '.html'))
        curr_subj = subj_dirs(iSubject);
        curr_subj_name = curr_subj.name
        curr_subj_dir = [curr_subj.folder '/' curr_subj.name];
        
        curr_output_dir = [outdir '/' curr_subj_name '/']
        mkdir(curr_output_dir)
        
        %%
        clear matlabbatch

        for iTask = 1:numel(tasks)
            
            curr_task = tasks{iTask};
            
            numRun = length(dir([curr_subj_dir '/ses-' curr_task '/' func_folder '/' '*MNI152NLin2009cAsym_desc-preproc_bold.nii.gz']));

            for iRun = 1:numRun

                ses_output = [curr_output_dir 'ses-' curr_task '/' func_folder '/'];
                mkdir(ses_output)
                so = dir(ses_output);
    
                if sum(contains({so.name}, 'smooth')) == numRun 
                    ierror = ['Error: Directory already smoothed! Skipping subject ' curr_subj_name]
                    errors = [errors; ierror];
                    continue;
                elseif sum(contains({so.name}, ['run-',mat2str(iRun)])) == 2
                    ierror = ['Error: Directory already smoothed! Skipping run for subject ' curr_subj_name ' run-' mat2str(iRun)]
                    errors = [errors; ierror];
                    continue;
                end
                
                 

                %%% get fMRI data
                curr_func_dir = [curr_subj_dir '/ses-' curr_task '/' func_folder];
                func_images = dir([curr_func_dir '/' '*MNI152NLin2009cAsym_desc-preproc_bold.nii'])
                
                if isempty(func_images) | length(func_images) ~= numRun
                    zipped_func = dir([curr_func_dir '/' '*MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'])
                    zipped_func = zipped_func(iRun)
                    if isempty(zipped_func)
                        ierror = ['Error: Not exactly 1 zipped.gz nifti found! Skipping session ' curr_subj_name '_' curr_task]
                        errors = [errors; ierror];
                        continue;
                    end
                    gzipfilenames = [zipped_func.folder '/' zipped_func.name];
                    gunzip(gzipfilenames, curr_func_dir)
                    func_images = dir([curr_func_dir '/' zipped_func.name(1:end-3)])
                else
                    func_images = func_images(iRun)
                end
                
                func_img_path = [func_images.folder '/' func_images.name]
                copyfile(func_img_path, ses_output)
    
                 %%% get functional volumes
    
                curr_func_image_path = [ses_output func_images.name];
                
                % load nifti & get number of volumes
                func_vols = spm_vol(curr_func_image_path);
                
                func_struct = {};
                for iVol = 1:numel(func_vols)
                    func_struct(iVol) = {[curr_func_image_path ',' num2str(iVol)]};
                end
                func_struct = func_struct'
    
                clear matlabbatch
            
                matlabbatch{1}.spm.spatial.smooth.data = func_struct;
                matlabbatch{1}.spm.spatial.smooth.fwhm = smooth_FWHM; % define in header
                matlabbatch{1}.spm.spatial.smooth.dtype = 0;
                matlabbatch{1}.spm.spatial.smooth.im = 0;
                matlabbatch{1}.spm.spatial.smooth.prefix = 'smooth';
            
                spm_jobman('run',matlabbatch)
            end
        end
    end
end

