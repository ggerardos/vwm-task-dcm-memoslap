%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Plot fROI activation magnitudes across task load + nonparametric stats
%
% Author: Gerasimos Gerardos
% Project: MeMoSLAP Project (P04–Leipzig)
% Paper: Task load modulates network interactions between bilateral fronto-
%        parietal and cerebellar areas during verbal working memory
%
% Description:
% This script visualizes and statistically compares fROI activation
% magnitudes across n-back load contrasts using subject-level ROI summary
% measures (stored as .mat files). For each ROI, it creates a bar + scatter
% plot showing mean activation per condition with 95% confidence intervals,
% overlays individual subject data points (jittered), runs a Friedman test
% (nonparametric repeated-measures ANOVA), performs post-hoc pairwise
% comparisons (Dunn–Sidak), and annotates significant pairwise differences.
%
% Inputs:
% - Per-ROI activation magnitude files in `root_dir`:
%     <ROI>.mat  (expected variable: `betas`, subjects × conditions)
%   where conditions correspond to contrasts (e.g., 1>0, 2>0, 3>0).
%
% Workflow (per ROI):
% 1) Load ROI data matrix (`betas`) containing subject-wise activation values
%    for each condition/contrast.
% 2) Compute:
%    - Mean across subjects for each condition
%    - 95% CI using 1.96 * (SD / sqrt(N))
% 3) Plot:
%    - Bar chart of condition means
%    - Scatter of individual subjects (with x-jitter)
%    - Error bars for 95% CI
% 4) Statistics:
%    - Friedman test across conditions (within-subject)
%    - Post-hoc pairwise comparisons using multcompare with Dunn–Sidak
%    - Annotate significant comparisons with *, **, *** (and optionally n.s.)
% 5) Format axes and save/export figures (export_fig optional).
%
% Outputs:
% - One plot per ROI saved in `plot_outdir` (if export is enabled), typically:
%     <ROI>_indivConds.png
% - Console output for Friedman test tables and post-hoc results.
%
% Notes:
% - `export_fig` is optional; figures can also be saved manually from the GUI.
% - The script currently defines custom bar colors; adjust as needed.
% - Ensure `betas` in each ROI .mat file matches the expected dimensions:
%   rows = subjects, columns = conditions.
% - Update `root_dir` and `plot_outdir` paths to match your environment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Clear the workspace
clear all
close all

%% Setup
addpath /Toolboxes/export_fig

root_dir = '/results/fROI_Analysis/activation_magnitudes/top10pct'; % Main directory for your data

ROIs = {'l_IFG';'r_IFG';'l_IPL'; 'r_IPL'; 'l_Cereb_sup'; 'r_Cereb_sup'}; %{'ROI_1'; 'ROI_2'; 'ROI_7'; 
    %'ROI_8'; 'ROI_9'; 'ROI_11'};
ROI_names = {'L IFG'; 'R IFG'; 'L IPL'; 'R IPL'; 'L Cer'; 'R Cer'};
     
cond_names = {'1 > 0', '2 > 0', '3 > 0'};  % Simplified condition names

plot_outdir = 'results/fROI_Analysis/activation_magnitudes/top10pct/Plots';

% Check if output directory exists
if ~exist(plot_outdir, 'dir')
    mkdir(plot_outdir);
end

%% Plot CG
for iROI = 1:numel(ROIs)

    curr_ROI = ROIs{iROI};
    curr_ROI_name = ROI_names{iROI};
    
    ROI_file = [root_dir '/' curr_ROI '.mat']; % Adjusted to load directly from root_dir
    load(ROI_file) % Load the ROI data
    data_CG = betas; % Assuming 'betas' is in your .mat file
      
    %% Plot: all conditions
    plot_name = 'indivConds';

    close all
    means_plot = mean(data_CG, 1); % Calculate mean across subjects
    errors_plot = (std(data_CG) ./ sqrt(size(data_CG, 1))) .* 1.96; % 95% CI

    % Create the plot
    figure('Position', [50 50 1000 800]); hold on;  % Create new figure
    bar_handles = bar(1:size(means_plot, 2), means_plot, 'LineWidth', 2);  % Bar plot for mean values

    % Define custom colors for each condition
    bar_colors = [
        0,0.447,0.698;  % Very light blue for 1back
        0,0.620,0.451;  % Light blue for 2back
        0.835,0.369,0   % Dark blue for 3back
    ]; 

    % Assign colors to bars directly without a loop
    % set(bar_handles(1), 'FaceColor', bar_colors(1, :));
    % set(bar_handles(2), 'FaceColor', bar_colors(2, :));
    % set(bar_handles(3), 'FaceColor', bar_colors(3, :));

    % Set each bar's color
    bar_handles.FaceColor = 'flat';           % Enable per-bar color
    bar_handles.CData = bar_colors;  

    % Plot scatter points and error bars for each condition
    for iCond = 1:size(data_CG, 2)  % Loop over each condition
        % Scatter points for individual subjects
        xData = iCond + randn(size(data_CG, 1), 1) * 0.1;  % Add jitter for visibility
        yData = data_CG(:, iCond);  % Get data for the current condition

        % Plot scatter points
        scatter(xData, yData, 40, ...
            'MarkerFaceColor', 'k', ... #bar_colors(iCond, :)
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceAlpha', 0.5, ...
            'MarkerEdgeAlpha', 0.5);

        % Add error bars for mean values
        errorbar(iCond, means_plot(iCond), errors_plot(iCond), 'k.', 'LineWidth', 2);
    end

    % Set axes properties
    set(gca, 'xtick', 1:numel(cond_names), 'xticklabel', cond_names);
    set(gca, 'FontSize', 12);
    ylabel('Activation magnitude (± 95% CI)');
    title(curr_ROI_name);
    set(gcf, 'Color', 'w');
    set(gca, 'Color', 'w');

    % Set y-axis limits
    ylim([-0.5,14])%([-0.5, max(max(max(data_CG) + errors_plot)) + 3]);  % Adjust the upper limit as needed  max(means_plot + errors_plot) + 1
    

    [p_friedman, tbl_friedman, stats_friedman] = friedman(data_CG, 1, 'off');
    disp(tbl_friedman);
    % Save current figure handle so it’s not affected
    main_fig = gcf;
    
    % Create a dummy figure for multcompare, but hide it
    tmp_fig = figure('Visible', 'off');
    
    % Run multcompare silently
    [pairs, means, pvals, labels] = multcompare(stats_friedman, 'CType', 'dunn-sidak');
    
    % Close dummy figure so your main figure stays untouched
    close(tmp_fig);
    
    % Reactivate the main figure in case needed
    figure(main_fig);
   
    
    C1 = pairs(:,1);
    C2 = pairs(:,2);
    Pv = pairs(:,6);
    % 
    sig_pairs = [C1,C2,Pv];
    
    y_max = max(ylim) + 0.3;
    y_spacing = 0.5;
    
    for i = 1:size(sig_pairs, 1)
        c1 = sig_pairs(i, 1);
        c2 = sig_pairs(i, 2);
        pval = sig_pairs(i, 3);
    
        y =  max(max(data_CG)) + 0.5  + (i - 1) * y_spacing;

        hold on
        
        
        ns_adj = 0;
        % Add asterisk or p-value
        if pval < 0.001
            txt = '***';
            % Draw horizontal line
            plot(sort([c1, c2]), [y y], 'k-', 'LineWidth', 1.5);
            hold on
            text(mean([c1 c2]), y + 0.12 + ns_adj, txt, 'HorizontalAlignment', 'center', 'FontSize', 28);
        elseif pval < 0.01
            txt = '**';
            % Draw horizontal line
            plot(sort([c1, c2]), [y y], 'k-', 'LineWidth', 1.5);
            hold on
            text(mean([c1 c2]), y + 0.12 + ns_adj, txt, 'HorizontalAlignment', 'center', 'FontSize', 28);
        elseif pval < 0.05
            txt = '*';
            % Draw horizontal line
            plot(sort([c1, c2]), [y y], 'k-', 'LineWidth', 1.5);
            hold on
            text(mean([c1 c2]), y + 0.12 + ns_adj, txt, 'HorizontalAlignment', 'center', 'FontSize', 28);
        else
            txt = 'n.s.';
            ns_adj = 0.1;
        end
    
        
    end
    
    %anova_p = ranovatbl.pValue(1);  % Get p-value from ANOVA
    anova_p = p_friedman;
    posthoc_method = 'Post hoc pairwise comparisons: Dunn-Sidak'; %Tukey–Kramer
    
        % Format p-value (e.g., 0.00012 → "< 0.001")
    if anova_p < 0.001
        anova_str = 'p < 0.001';
    else
        anova_str = sprintf('p = %.3f', anova_p);
    end
    
    % Combine text
    stats_text = sprintf('Friedman test %s\n%s', anova_str, posthoc_method); %One-way repeated-measures ANOVA

    % Place text on the plot — adjust (x, y) as needed
    %text(2, max(ylim)-0.2, stats_text, ...
   % 'FontSize', 12, ...
   % 'VerticalAlignment', 'top', ...
   % 'BackgroundColor', 'white', ...
   % 'EdgeColor', 'none');

    yticks(0:2:16)

    ylim([min(ylim), y + 0.3]);  % Extend y-limits
    
    %xtickangle(20)

    ax = gca; % get current axes
    ax.FontSize = 24;

    %% Export figure or save from gui
    %export_fig(gcf, [plot_outdir '/' curr_ROI '_' plot_name '.png'],'-png', '-r1000','-p0.05');
    close all;  % Close the figure after saving
end

%% Export data for stats programs
for iROI = 1:numel(ROIs)

    curr_ROI = ROIs{iROI};
    curr_ROI_name = ROI_names{iROI};
    
    ROI_file = [root_dir '/' curr_ROI '.mat']; % Load each ROI's data from root_dir
    load(ROI_file);
    data_CG = betas;

    % You can add code here to save or process data_CG as needed
end



