%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Build nuisance regressor matrices from fMRIPrep confounds + QC
%
% Author: Gerasimos Gerardos
% Project: MeMoSLAP Project (P04–Leipzig)
% Paper: Task load modulates network interactions between bilateral fronto-
%        parietal and cerebellar areas during verbal working memory
%
% Description:
% This script creates session-wise nuisance regressor matrices (R) for SPM
% analyses using fMRIPrep confounds, motion QC, and task-related regressors.
% It processes each participant and session, generates motion/FD QA figures,
% builds censoring regressors for high-motion volumes, optionally adds
% aCompCor components, encodes block order information across volumes, and
% includes task-performance outlier blocks (from a d' outlier table) as
% additional nuisance regressors.
%
% Inputs:
% - fMRIPrep confounds files (*confounds_timeseries.tsv) per subject/session
% - Block start and duration CSVs (timings_SPM/block_starts, block_durations)
% - Outlier block table (outliar_blocks.csv) containing poor-performance blocks
%
% Key steps:
% 1) Load fMRIPrep confounds and build motion regressors (6 or 24 params).
% 2) Generate motion parameter plots (translation/rotation) for QA.
% 3) Compute framewise displacement (FD), flag volumes with FD > threshold,
%    and create one-hot censoring regressors (spike regressors).
% 4) Plot FD time series for QA and store summary FD statistics per run.
% 5) Optionally append aCompCor regressors (top N components).
% 6) Use block start/duration tables to compute per-block volume ranges and
%    create an “orders” regressor that encodes the within-run block sequence.
% 7) If behavioral outlier blocks exist for a subject/session, add one-hot
%    regressors marking the volumes of those blocks.
% 8) Add a session regressor to center sessions across runs.
% 9) Replace NaNs with zeros and save the nuisance matrix R as a .mat file.
% 10) Save/update a CSV summary of FD outliers and descriptive FD statistics.
%
% Outputs:
% - Nuisance regressor .mat per subject/session: [outdir '/sub*_ses-*.mat']
%     variable: R (volumes x regressors)
% - QA figures (PNG): motion_*.png and FD_*.png in outdir_motionFigures
% - FD summary table: FD_outliers.csv in outdir_motionFigures
%
% Notes:
% - Update directory paths (fmriprep_dir, onset_dir, dur_dir, outdir, etc.)
%   before running.
% - Set options at the top (n_motion_regs, FD_threshold, aCompCor, n_aCompCor).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear all
close all

%%
fmriprep_dir = '/results/fMRIprep/fmap_gr1/';

sub_dirs = dir([fmriprep_dir '/sub-*ham*']);

onset_dir = '/results/fMRI_analysis/timings_SPM/block_starts/'; %change to your directory 
dur_dir = '/results/fMRI_analysis/timings_SPM/block_durations/'; %change to your directory 
outblock_dir = '/results/task_analysis/outliar_blocks.csv';

start_dirs = dir([onset_dir '/sub-*']);
last_dirs = dir([dur_dir '/sub-*']);
outblocks = readtable(outblock_dir, ReadVariableNames=true);

spm_path = '/data/u_gerados_software/MATLAB/spm12'; % save SPM12 path for later
addpath(spm_path)
%spm fmri

tasks = {'nback1', 'nback2'};

n_motion_regs = 24; % 6 or 24 motion regressors?

FD_threshold = 0.9; % choose FD threshold

aCompCor = 1; % include aCompCor regressors?
n_aCompCor = 10; % how many of the top aCompCor regressors?

outdir = '/results/fMRI_analysis/nuisance_regressors/motionRegs24_FD09_aCompCor10_order_withsess';
if isempty(dir(outdir))
    mkdir(outdir)
end

outdir_motionFigures = '/results/fMRI_analysis/QA_motion';
if isempty(dir(outdir_motionFigures))
    mkdir(outdir_motionFigures)
end

%%
FD_table = table;
FD_row = 1;

for iSub = 1:numel(sub_dirs)
    sub_name = sub_dirs(iSub).name;
    sub_code = arrayfun(@(x) str2double(x), sub_name);
    sub_code = arrayfun(@(x)  num2str(x) ,sub_code(~isnan(sub_code)));
   
    %%
    if sub_dirs(iSub).isdir == 1 && any(cellfun(@(x) all(join(["sub-ham",sub_code "_anat.html"], "",2)==x | join([sub_dirs(iSub).name ".html"], "",2)==x), {sub_dirs.name})) % fix second part
        curr_sub = sub_dirs(iSub).name
        curr_folder = [fmriprep_dir curr_sub];
        
        
        for iTask = 1:numel(tasks)
            
            curr_task = tasks{iTask};

            curr_func_folder = [curr_folder '/ses-' curr_task '/func'];
        
            confounds_file = dir([curr_func_folder '/*confounds_timeseries.tsv']);

            if numel(confounds_file) ~= 1
                [curr_func_folder ': Not exactly 1 confounds_timeseries.tsv file found! Continuing to the next session']
                continue; 
            end

            if ~isempty(dir([outdir '/' curr_sub '_ses-' curr_task '.mat']))
                [curr_func_folder ' has alredy been processed. Continuing to next loop']
                continue;
            end


            confounds_file_path = [confounds_file.folder '/' confounds_file.name];

            %% Import fmriprep confound file 
            confounds = readtable(confounds_file_path,'FileType','text');
            
            %% motion regressors
            
            % 6 motion regs
            if n_motion_regs == 6
                motion_regs = [confounds.trans_x confounds.trans_y confounds.trans_z ...
                    confounds.rot_x confounds.rot_y confounds.rot_z];
                
                idx_trans = 1:3;
                idx_rot = 4:6;
                
            % 24 motion regs
            elseif n_motion_regs == 24
                motion_regs = [confounds.trans_x confounds.trans_x_derivative1 ...
                    confounds.trans_x_power2 confounds.trans_x_derivative1_power2 ...
                    confounds.trans_y confounds.trans_y_derivative1 ...
                    confounds.trans_y_derivative1_power2 confounds.trans_y_power2 ...
                    confounds.trans_z confounds.trans_z_derivative1 ...
                    confounds.trans_z_derivative1_power2 confounds.trans_z_power2 ...
                    confounds.rot_x confounds.rot_x_derivative1 ...
                    confounds.rot_x_power2 confounds.rot_x_derivative1_power2 ...
                    confounds.rot_y confounds.rot_y_derivative1 ...
                    confounds.rot_y_derivative1_power2 confounds.rot_y_power2 ...
                    confounds.rot_z confounds.rot_z_derivative1 ...
                    confounds.rot_z_derivative1_power2 confounds.rot_z_power2];
                
                idx_trans = [1 5 9];
                idx_rot = [13 17 21];
            end
            
            %% plot motion regressors
            yscale_trans = [-6 6]
            yscale_rot = [-6 6];
            
             printfig = figure;
             set(printfig, 'Name', ['Motion parameters: ' curr_sub ' - Ses:' curr_task], 'Visible', 'on');
             %loadmot = load(deblank(b(i,:)));
             subplot(2,1,1);
             plot(motion_regs(:,idx_trans));
             line([1:size(motion_regs,1)], 5.*ones(size(motion_regs,1)), 'Color','k')
             line([1:size(motion_regs,1)], -5.*ones(size(motion_regs,1)), 'Color','k')
             grid on;
             ylim(yscale_trans);  % enable to always scale between fixed values as set above
             xlim([0 size(motion_regs,1)])
             title('Translation') 
             xlabel('Image')
             ylabel('mm')

             subplot(2,1,2);
             plot(motion_regs(:,idx_rot)*180/pi);
             line([1:size(motion_regs,1)], 5.*ones(size(motion_regs,1)), 'Color','k')
             line([1:size(motion_regs,1)], -5.*ones(size(motion_regs,1)), 'Color','k')
             grid on;
             ylim(yscale_rot);   % enable to always scale between fixed values as set above
             xlim([0 size(motion_regs,1)])
             title('Rotation');
             xlabel('Image')
             ylabel('degrees')

             motname = [outdir_motionFigures '/motion_' curr_sub '_ses-' curr_task '.png']
             print(printfig, '-dpng', '-noui', '-r300', motname);  % enable to print to file
             close all

            %% FD censoring
            FD = confounds.framewise_displacement;

            idx_motOutlier = find(FD > FD_threshold);
            num_outliers = numel(idx_motOutlier);
            
            % save outlier info in matrix
            FD_table.Sub(FD_row) = {[curr_sub '_' curr_task]};
            FD_table.Outliers(FD_row) = num_outliers;
            FD_table.max_FD(FD_row) = max(FD);
            FD_table.mean_FD(FD_row) = nanmean(FD);
            FD_table.SD_FD(FD_row) = nanstd(FD);

            FD_row = FD_row + 1;

            % create censoring regressors
            motion_outliers = [];
            for iOutlier = 1:numel(idx_motOutlier)

                col_outlier = zeros(numel(FD),1);
                col_outlier(idx_motOutlier(iOutlier)) = 1;
                motion_outliers = [motion_outliers col_outlier];

            end
            
            %% plot FD
            figure_FD = figure;
            plot(FD)
            line([1:size(FD,1)], 0.9 .* ones(size(FD,1)), 'Color','r')
            ylims_FD = [0 3];

            
            ylim(ylims_FD)
            xlim([[0 size(FD,1)]])
            grid on
            title(['FD: ' curr_sub ' - Task:' curr_task])
            xlabel('Volume')
            ylabel('Framewise displacement (FD)')
            
            outname = [outdir_motionFigures '/FD_' curr_sub '_ses-' curr_task '.png']
            print(figure_FD, '-dpng', '-noui', '-r300', outname);  % enable to print to file
            close all
            
            %% aCompCor
            
            aCompCor_regs = [];
            if aCompCor == 1
                
                for iComp = 1:n_aCompCor
                
                    aCompCor_regs = [aCompCor_regs confounds.(['a_comp_cor_' num2str(iComp-1,'%02.f')])];
                                        
                end
            
            end

            %% order of blocks
            
            start_dir = start_dirs(cellfun(@(x) contains(x, char(extract(curr_sub, digitsPattern))) & (contains(x, ['task-' num2str(iTask)]) | contains(x, ['task' num2str(iTask)])), {start_dirs.name}))
            last_dir = last_dirs(cellfun(@(x) contains(x, char(extract(curr_sub, digitsPattern))) & (contains(x, ['task-' num2str(iTask)]) | contains(x, ['task' num2str(iTask)])), {last_dirs.name}))

            starts = readtable([start_dir.folder '/' start_dir.name], ReadVariableNames=true);
            startn = starts(1:4,1:6);
            startn = table2array(startn);
            startn = startn(:);
            [startn, istart] = sort(startn(~isnan(startn)), 'ascend');

            lasts = readtable([last_dir.folder '/' last_dir.name], ReadVariableNames=true);
            lastn = lasts(1:4,1:6);
            lastn = table2array(lastn);
            lastn = lastn(:);
            lastn = lastn(~isnan(lastn));
            lastn = lastn(istart);

            endn = startn + lastn;

            volstart = ceil(startn);
            volend = floor(endn);
            
            orders = zeros(1420,1);
            for v=1:length(volstart)
                orders(volstart(v):volend(v))= v-mean(1:length(volstart))+0.5;
            end
            
            
          
            %% put stuff together
            R = [motion_regs motion_outliers aCompCor_regs orders];

            %% add outliars if any
            part_code = extract(curr_sub, digitsPattern);
            if any(outblocks.code == str2double(part_code{:})) && any(outblocks.session == iTask)
                curr_out = outblocks(outblocks.code == str2double(part_code{:}) & outblocks.session == iTask,:);
                outliars_reg = zeros(1420,size(curr_out,1));
                for o = 1:size(curr_out,1)
                    iout = outblocks.order(o)
                    outliars_reg(volstart(iout):volend(iout),o) = 1;
                end
                R = [R outliars_reg];
            end

            %% session regressor
           session_reg = ones(1420,1)*(iTask-mean(1:length(tasks)));
           R = [R session_reg];



            %% turn NaNs into 0s
            R(isnan(R)) = 0;
            
            %% save
            save([outdir '/' curr_sub '_ses-' curr_task '.mat'], 'R')
            
        end
        
    end
end

%% save FD matrix
if ~isempty(dir([outdir_motionFigures '/FD_outliers.csv']))
    FD_old = readtable([outdir_motionFigures '/FD_outliers.csv']);
    FD_update = vertcat(FD_old, FD_table)
    FD_update = sortrows(FD_update, 1)
    writetable(FD_update, [outdir_motionFigures '/FD_outliers.csv'])
else
    writetable(FD_table, [outdir_motionFigures '/FD_outliers.csv'])
end

