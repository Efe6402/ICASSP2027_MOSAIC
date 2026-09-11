%% PLOT_VIEWWISE_KAPPA_SNR_CURVES.m
% Re-plot the completed VIEW-WISE heteroscedastic kappa x SNR experiment.
%
% This script DOES NOT rerun any solver and DOES NOT change any hyperparameter,
% threshold, tau, graph, or noise realization.  It only reads TEST_SUMMARY.csv
% from the completed run and creates SNR-sweep curves for each kappa.
%
% Output:
%   FIGURES_SNR_BY_KAPPA/
%       VIEWWISE_KAPPA_SNR_TOP1_F1.png
%       VIEWWISE_KAPPA_SNR_TOP1_F1.pdf
%       VIEWWISE_KAPPA_SNR_TOP1_F1.fig
%
% The figure layout matches the intended interpretation:
%   top row    : Supra F1 versus SNR
%   bottom row : Mean(within,cross,copy) F1 versus SNR
%   columns    : MOSAIC, PGL2021, ZW2024
%   lines      : one line for each kappa
%
% Efe's completed result directory:
clear; clc; close all;

%% ------------------------------------------------------------------------
% 1. Completed result folder
% -------------------------------------------------------------------------
run_dir = ...
'/Users/efekarakoca/Desktop/university/terms/senior/selin_hoca_research/experiments/MOSAIC-experiments/mosaic-synthetic-12-august/MOSAIC_viewwise_heteroscedastic_kappa_snr_sweep_fixed_p5000_v2fixed/results/viewwise_kappa_snr_sweep_p05000/run_full_20260902_193816';

assert(isfolder(run_dir), ...
    'Result directory does not exist:\n%s', run_dir);

summary_file = fullfile(run_dir,'TEST_SUMMARY.csv');

% If the final table was not written because the run stopped at the very end,
% fall back to the checkpoint table.
if ~isfile(summary_file)
    checkpoint_file = fullfile(run_dir,'TEST_SUMMARY_CHECKPOINT.csv');
    assert(isfile(checkpoint_file), ...
        ['Could not find TEST_SUMMARY.csv or TEST_SUMMARY_CHECKPOINT.csv in:\n' ...
         '%s'], run_dir);
    summary_file = checkpoint_file;
end

fprintf('Reading:\n%s\n\n', summary_file);

T = readtable(summary_file,'TextType','string');

%% ------------------------------------------------------------------------
% 2. Validate the columns needed for plotting
% -------------------------------------------------------------------------
required_vars = { ...
    'method','selection_rank','candidate_id','kappa','target_snr_db', ...
    'f1_supra_mean','f1_component_mean'};

for q = 1:numel(required_vars)
    assert(ismember(required_vars{q},T.Properties.VariableNames), ...
        'Required column "%s" is missing from %s.', ...
        required_vars{q}, summary_file);
end

% We want the frozen rank-one configuration for each method, exactly as in
% the example figure.  No model re-selection is performed here.
T = T(T.selection_rank == 1,:);

methods = ["MOSAIC","PGL2021","ZW2024"];
for m = 1:numel(methods)
    assert(any(T.method == methods(m)), ...
        'No rank-one rows found for method %s.', methods(m));
end

%% ------------------------------------------------------------------------
% 3. Recover kappa and SNR grids directly from the saved results
% -------------------------------------------------------------------------
kappa_grid = sort(unique(double(T.kappa(:))))';

% Clean is stored as +Inf.  Put it first, followed by finite SNR values in
% descending order: clean, 20, 15, 12, ..., 0.
snr_all = unique(double(T.target_snr_db(:)));
has_clean = any(isinf(snr_all) & snr_all > 0);
snr_finite = sort(snr_all(isfinite(snr_all)),'descend')';

if has_clean
    snr_grid = [Inf, snr_finite];
else
    snr_grid = snr_finite;
end

x = 1:numel(snr_grid);
xlabels = strings(size(snr_grid));
for s = 1:numel(snr_grid)
    if isinf(snr_grid(s))
        xlabels(s) = "clean";
    else
        xlabels(s) = string(snr_grid(s));
    end
end

fprintf('Kappa grid found: %s\n', mat2str(kappa_grid));
fprintf('SNR grid found  : ');
disp(xlabels);

%% ------------------------------------------------------------------------
% 4. Output folder
% -------------------------------------------------------------------------
fig_dir = fullfile(run_dir,'FIGURES_SNR_BY_KAPPA');
if ~isfolder(fig_dir)
    mkdir(fig_dir);
end

%% ------------------------------------------------------------------------
% 5. Plot: one SNR curve for each kappa
% -------------------------------------------------------------------------
fig = figure( ...
    'Color','w', ...
    'Position',[40 40 1800 980], ...
    'Name','View-wise kappa x SNR test F1');

tl = tiledlayout(fig,2,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');

title(tl, ...
    sprintf(['p=5000 | fixed p=5000 rank-one configurations | ' ...
             'VIEW-WISE \\kappa and SNR test F1']), ...
    'FontWeight','bold', ...
    'FontSize',15);

% Let MATLAB choose distinct default colors.  We keep the same color assigned
% to a kappa in every panel.
color_order = lines(numel(kappa_grid));

metric_names = ["f1_supra_mean","f1_component_mean"];
row_titles   = ["Supra F1","Mean(within,cross,copy) F1"];

legend_handles = gobjects(numel(kappa_grid),1);

for row = 1:2
    metric = metric_names(row);

    for m = 1:numel(methods)
        ax = nexttile(tl,(row-1)*3+m);
        hold(ax,'on');
        grid(ax,'on');
        box(ax,'on');

        method_rows = T(T.method == methods(m),:);

        for ki = 1:numel(kappa_grid)
            kap = kappa_grid(ki);
            y = nan(size(snr_grid));

            for si = 1:numel(snr_grid)
                if isinf(snr_grid(si))
                    Q = method_rows( ...
                        abs(double(method_rows.kappa)-kap) < 1e-12 & ...
                        isinf(double(method_rows.target_snr_db)), :);
                else
                    Q = method_rows( ...
                        abs(double(method_rows.kappa)-kap) < 1e-12 & ...
                        abs(double(method_rows.target_snr_db)-snr_grid(si)) < 1e-12, :);
                end

                if ~isempty(Q)
                    vals = double(Q.(metric));
                    vals = vals(isfinite(vals));
                    if ~isempty(vals)
                        % There should normally be exactly one aggregate row.
                        % mean() keeps this plotting script robust if duplicate
                        % checkpoint rows exist.
                        y(si) = mean(vals,'omitnan');
                    end
                end
            end

            h = plot(ax,x,y, ...
                '-o', ...
                'LineWidth',1.8, ...
                'MarkerSize',5.5, ...
                'Color',color_order(ki,:), ...
                'MarkerFaceColor','none');

            if row == 1 && m == 1
                legend_handles(ki) = h;
            end
        end

        title(ax,sprintf('%s | %s',methods(m),row_titles(row)), ...
            'FontWeight','bold');

        xlabel(ax,'prescribed SNR (dB)');
        ylabel(ax,'F1');

        xticks(ax,x);
        xticklabels(ax,xlabels);
        xlim(ax,[0.85 numel(x)+0.15]);
        ylim(ax,[0 1.02]);

        ax.FontSize = 10;
        ax.Layer = 'top';
    end
end

% One common legend for all panels.
legend_labels = "kappa=" + string(kappa_grid);
lgd = legend(legend_handles,legend_labels, ...
    'Orientation','horizontal', ...
    'Location','southoutside');
lgd.Layout.Tile = 'south';

%% ------------------------------------------------------------------------
% 6. Save publication-quality versions
% -------------------------------------------------------------------------
png_file = fullfile(fig_dir,'VIEWWISE_KAPPA_SNR_TOP1_F1.png');
pdf_file = fullfile(fig_dir,'VIEWWISE_KAPPA_SNR_TOP1_F1.pdf');
fig_file = fullfile(fig_dir,'VIEWWISE_KAPPA_SNR_TOP1_F1.fig');

exportgraphics(fig,png_file,'Resolution',300);
exportgraphics(fig,pdf_file,'ContentType','vector');
savefig(fig,fig_file);

fprintf('\nSaved:\n');
fprintf('  %s\n',png_file);
fprintf('  %s\n',pdf_file);
fprintf('  %s\n',fig_file);

%% ------------------------------------------------------------------------
% 7. Optional audit table: exact plotted values
% -------------------------------------------------------------------------
% Save the rank-one rows used by this plotting script so the figure can be
% audited without rerunning anything.
plot_table = sortrows(T,{'method','kappa','target_snr_db'});
audit_file = fullfile(fig_dir,'VIEWWISE_KAPPA_SNR_TOP1_F1_PLOTTED_ROWS.csv');
writetable(plot_table,audit_file);

fprintf('  %s\n',audit_file);
fprintf('\nDone. No solver or hyperparameter file was modified.\n');
