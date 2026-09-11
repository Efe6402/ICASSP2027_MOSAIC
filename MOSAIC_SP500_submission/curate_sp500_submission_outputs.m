function curate_sp500_submission_outputs(runRoot,outRoot)
%CURATE_SP500_SUBMISSION_OUTPUTS Collect the selected model and analyses.

selectionRoot = fullfile(runRoot,'01_OFFSET_VALIDATION_RANK_HYPERPARAMETER_SELECTION');
targetRoot = fullfile(runRoot,'02_FINAL_OFFSET0_TARGET_TOP1');
summaryRoot = fullfile(runRoot,'03_MASTER_SUMMARY');

candidate = dir(fullfile(targetRoot,'top01_*'));
candidate = candidate([candidate.isdir]);
assert(numel(candidate)==1,'Expected one selected target-fit directory.');
candidateRoot = fullfile(candidate(1).folder,candidate(1).name);

figRoot = fullfile(candidateRoot,'figures');
postRoot = fullfile(candidateRoot,'cross_within_mode_postprocess');
viewRoot = fullfile(candidateRoot,'view_graph_clustering');

make_dir(outRoot);
make_dir(fullfile(outRoot,'selection'));
make_dir(fullfile(outRoot,'model'));
make_dir(fullfile(outRoot,'figures'));
make_dir(fullfile(outRoot,'tables'));
make_dir(fullfile(outRoot,'matrices'));
make_dir(fullfile(outRoot,'view_graph'));

packageRoot = fileparts(mfilename('fullpath'));
copy_required(fullfile(packageRoot,'selected_top1','figures', ...
    '01_SP500_paper_figure.png'), ...
    fullfile(outRoot,'figures','00_SP500_paper_figure.png'));

copy_required(fullfile(runRoot,'MASTER_RUN_LOG.txt'), ...
    fullfile(outRoot,'RUN_LOG.txt'));
copy_required(fullfile(summaryRoot,'FINAL_RANK_SELECTION_SUMMARY.csv'), ...
    fullfile(outRoot,'selection','rank_selection_summary.csv'));
copy_required(fullfile(summaryRoot,'SELECTED_TOP1_CONFIGURATION.csv'), ...
    fullfile(outRoot,'selection','SELECTED_TOP1_CONFIGURATION.csv'));
copy_required(fullfile(summaryRoot,'FINAL_TOP1_TARGET_FITS.csv'), ...
    fullfile(outRoot,'selection','SELECTED_TOP1_FIT_SUMMARY.csv'));
copy_required(fullfile(selectionRoot, ...
    'HYPERPARAMETER_BANK_REFERENCE64_PLUS_HP59_NEIGHBORS.csv'), ...
    fullfile(outRoot,'selection','HYPERPARAMETER_BANK.csv'));
copy_required(fullfile(selectionRoot,'selection_profile_rank_aggregate.csv'), ...
    fullfile(outRoot,'selection','PROFILE_RANK_VALIDATION_SCORES.csv'));

copy_required(fullfile(candidateRoot,'selected_representative_fit.mat'), ...
    fullfile(outRoot,'model','selected_top1_model.mat'));
copy_required(fullfile(targetRoot, ...
    'PRIMARY_FINAL_MODEL_TOP1_AUXILIARY_VALIDATION.mat'), ...
    fullfile(outRoot,'model','selected_top1_compact_model.mat'));

figureMap = {
    fullfile(figRoot,'02_modal_profiles_sector_graphs_copy.png'), ...
        '01_modal_profiles_stock_graphs_and_copy.png';
    fullfile(figRoot,'08_sector_diagnostics.png'), ...
        '02_within_view_sector_connectivity.png';
    fullfile(postRoot,'figures','01_top_cross_vs_within_edge_weights.png'), ...
        '03_cross_view_relations.png';
    fullfile(viewRoot,'01_view_graph_components.png'), ...
        '04_induced_view_graph.png';
    fullfile(viewRoot,'05_eigengap_clustering.png'), ...
        '05_view_graph_clustering.png';
    fullfile(viewRoot,'02_spectral_cluster_diagnostics.png'), ...
        '06_view_graph_eigengap.png';
    fullfile(viewRoot,'06_modal_view_graph_decomposition.png'), ...
        '07_modal_view_graph_contributions.png'};
copy_map(figureMap,fullfile(outRoot,'figures'));

tableMap = {
    fullfile(candidateRoot,'tables','within_horizon_sector_statistics.csv'), ...
        'within_horizon_sector_statistics.csv';
    fullfile(candidateRoot,'tables','mode_sector_statistics.csv'), ...
        'mode_sector_statistics.csv';
    fullfile(postRoot,'tables','top_cross_relations_cross_vs_within.csv'), ...
        'top_cross_relations_cross_vs_within.csv';
    fullfile(postRoot,'tables','top_cross_relations_cross_vs_within_ranking.csv'), ...
        'top_cross_relations_cross_vs_within_ranking.csv'};
copy_map(tableMap,fullfile(outRoot,'tables'));

matrixNames = {'B_view_mode_profiles.csv','Gamma_stock_copy_coefficients.csv', ...
    'A_mode_01.csv','A_mode_02.csv','A_mode_03.csv', ...
    'stock_index_metadata.csv','view_index_metadata.csv'};
for q=1:numel(matrixNames)
    copy_required(fullfile(postRoot,'matrices',matrixNames{q}), ...
        fullfile(outRoot,'matrices',matrixNames{q}));
end

viewNames = {'G_view_total.csv','G_view_crossnode_contribution.csv', ...
    'G_view_copy_contribution.csv','view_graph_edge_list.csv', ...
    'view_cluster_assignments.csv','spectral_cluster_metrics.csv', ...
    'modal_view_masses.csv'};
for q=1:numel(viewNames)
    copy_required(fullfile(viewRoot,viewNames{q}), ...
        fullfile(outRoot,'view_graph',viewNames{q}));
end
end

function copy_map(map,destination)
for q=1:size(map,1)
    copy_required(map{q,1},fullfile(destination,map{q,2}));
end
end

function copy_required(source,destination)
assert(isfile(source),'Missing expected output: %s',source);
[ok,msg] = copyfile(source,destination,'f');
assert(ok,'Could not copy %s: %s',source,msg);
end

function make_dir(path)
if ~isfolder(path), mkdir(path); end
end
