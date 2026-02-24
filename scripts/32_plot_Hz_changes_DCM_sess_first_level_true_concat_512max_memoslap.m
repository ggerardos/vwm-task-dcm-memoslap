%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Summarize and plot intrinsic and session-modulated connectivity
%         from PEB/BMA DCM group results
%
% Purpose:
% This script loads the exported DCM group results table
% (results_table.xlsx) and summarizes intrinsic (A-matrix) and
% session-related modulatory (B-matrix) connectivity effects.
% It combines baseline and session-specific effects and visualizes
% reliable connections.
%
% After estimating subject DCMs, building a PEB model, and running
% BMR/BMA, this script:
%   1) Identifies significant intrinsic connections (A-matrix) based on
%      posterior probability,
%   2) Extracts corresponding Session 1 and Session 2 modulatory effects
%      (B-matrix, Hz_change),
%   3) Combines intrinsic and modulatory effects into summary tables,
%   4) Selects connections showing reliable session modulation,
%   5) Plots intrinsic, Session 1, and Session 2 connectivity values.
%
% Inputs:
% - results_table.xlsx
%     Exported PEB/BMA summary table containing:
%       * Field ('A','B','C')
%       * Connection_name
%       * ParameterEstimate
%       * PosteriorProbability
%       * Hz_change
%
% What it does (high level):
% 1) Loads the DCM results table.
% 2) Separates intrinsic (A) and modulatory (B) parameters.
% 3) Selects reliable intrinsic connections using posterior probability.
% 4) Matches modulatory parameters to intrinsic connections.
% 5) Builds tables of:
%      - Intrinsic connectivity
%      - Session 1 modulation
%      - Session 2 modulation
% 6) Identifies connections with significant session effects.
% 7) Produces stem/scatter plots of connectivity changes.
%
% Outputs:
% - Figure(s) showing intrinsic, Session 1, and Session 2 connectivity (Hz)
%   for selected connections.
%
% Key variables:
% - iSig
%     Logical index for significant intrinsic connections.
%
% - conns
%     Table containing intrinsic, sess1, and sess2 connectivity values.
%
% - sig_conns
%     Matrix containing posterior probabilities for each effect.
%
% Notes / Caveats:
% - Significance is currently defined as:
%       PosteriorProbability == 1
%   (can be changed to >= 0.95 if needed).
%
% - String indexing of Connection_name assumes a fixed naming format.
%   If naming conventions change, this part must be updated.
%
% - Only the first selected connection is plotted by default.
%   Modify the plotting loop to visualize all connections.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dcm_results = '/results/DCM_models/DCM_sess_mod_restrict_512max/results_table.xlsx';
dcm_data = readtable(dcm_results);

iA = cell2mat(dcm_data.Field) == 'A';
iB = cell2mat(dcm_data.Field) == 'B';
iC = cell2mat(dcm_data.Field) == 'C';
iSig = dcm_data.PosteriorProbability == 1;

conn_names = cellfun(@(x) x(16:end), dcm_data.Connection_name(iA & iSig), UniformOutput=false);

conns = table(zeros(length(conn_names), 1),zeros(length(conn_names), 1),zeros(length(conn_names), 1),...
    'RowNames',conn_names, 'VariableNames',{' intrinsic', 'sess1', 'sess2'})
sig_conns = zeros(26,3)

connB_names = cellfun(@(x) x(16:end-12), dcm_data.Connection_name(iB), UniformOutput=false);
connB_unique = unique(connB_names)

dcm_data.Connection_name(iB)=connB_names;


for i=1:numel(conn_names)
    name = conn_names{i};
    iname = cellfun(@(x) strcmp(x, name), dcm_data.Connection_name, UniformOutput=true);
    if sum(iname)==0
        continue
    else
        conn = dcm_data.Hz_change(iname & iB);
        sig = dcm_data.PosteriorProbability(iname & iB);
        conns(i,2:3)= num2cell(conn');
        sig_conns(i,2:3)= sig';
    end
end
        


conns(:,1) = num2cell(dcm_data.ParameterEstimate(iA & iSig));
sig_conns(:,1) = ones(26,1);


conn2plot = conns(any(sig_conns(:,2:3),2),:);
sig2plot = sig_conns(any(sig_conns(:,2:3),2),:);
sig2plot = sig2plot >= 0.95;

conn_mat = table2array(conn2plot);

names2plot = strrep(conn2plot.Properties.RowNames, '_', '-');


figure
nlevels =0:2;

colors= [0,0,0;
        0.7569,0.0980,0.5843;  % Very light blue for 1back
        1,0.7529,0;  % Light blue for 2back
        0.835,0.369,0   % Dark blue for 3back
    ]; 

for i = 1

    %subplot(3,3,i)
    stem(nlevels, conn_mat(i,:),'k', 'LineWidth', 1.75)
     % Add title to each subplot
    title(names2plot(i))

    hold on

    for n = nlevels
        scatter(n, conn_mat(i,n+1), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(n+1,:), 'SizeData', 81)
        hold on 
    end
       

    % Add x and y labels
    xlabel('Condition')
    ylabel('Connectivity (Hz)')

    hold on 

    scatter(nlevels(sig2plot(i,:)), [conn_mat(i,sig2plot(i,:)).*ones(1, length(nlevels(sig2plot(i,:))))]+0.15*sign(conn_mat(i,sig2plot(i,:))), '*','k')
    %max(conn_mat(i,:))*ones(length(nlevels(sig2plot(i,:))))+0.03
    xticks(nlevels)
    xticklabels({'Intrinsic', 'Session 1', 'Session 2'})
    yticks(-0.4:0.2:0.8)
    
    hold on 
    yline(0, 'k', 'LineWidth', 1, Alpha=0.3);
    
    hold on
    xlim([-0.1, 2.1])
    ylim([-0.1, 1]) %([min(conn_mat(i,:))-0.03,max(conn_mat(i,:))+0.06])
end

set(gcf, 'Color', 'w');
set(gca, 'Color', 'w');
