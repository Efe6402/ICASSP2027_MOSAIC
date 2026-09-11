%% Verify the S&P 500 submission package and prepared offset-0 dataset
clearvars; clc;
root = fileparts(mfilename('fullpath'));

required = {
    'README.md'
    'RUN_SP500_SUBMISSION.m'
    'RUN_DATA_PREPARATION.m'
    'run_sp500_submission_core.m'
    'curate_sp500_submission_outputs.m'
    'convert_rds_to_csv_with_dates.R'
    fullfile('data','build_mosaic_finance_multihorizon_dataset.m')
    fullfile('data','OFFSET_DATASET_DESIGN.csv')
    fullfile('data','stocks','sp500-data-2016-2020.rds')
    fullfile('data','stocks','sp500-data-2016-2020_converted.csv')
    fullfile('data','stocks','SP500-sectors.csv')
    fullfile('data','stocks','mosaic_selected_tickers.csv')
    fullfile('data','prepared_offset0','mosaic_finance_multihorizon_dataset.mat')
    fullfile('data','prepared_offset0','common_ending_dates.csv')
    fullfile('data','prepared_offset0','selected_stocks_and_sectors.csv')
    fullfile('selected_top1','SELECTED_TOP1_CONFIGURATION.csv')
    fullfile('selected_top1','selected_top1_model.mat')
    fullfile('selected_top1','figures','01_SP500_paper_figure.png')};

for q=1:numel(required)
    assert(isfile(fullfile(root,required{q})),'Missing: %s',required{q});
end

core = fileread(fullfile(root,'run_sp500_submission_core.m'));
main = fileread(fullfile(root,'RUN_SP500_SUBMISSION.m'));

assert(contains(core,'cfg.seed = 20260813;'),'Unexpected experiment seed.');
assert(contains(core,'cfg.rank_candidates = 2:5;'),'Unexpected rank search.');
assert(contains(core,'cfg.top_k_profiles = 1;'),'Final selection is not top-1.');
assert(contains(core,'cfg.transfer_pairs = [2 3; 4 5];'), ...
    'Unexpected offset-validation configuration.');
assert(contains(main,'prepared_offset0'), ...
    'Main experiment does not reference the prepared target.');
datasetPath = fullfile(root,'data','prepared_offset0', ...
    'mosaic_finance_multihorizon_dataset.mat');
loaded = load(datasetPath,'data');
assert(isfield(loaded,'data') && isstruct(loaded.data), ...
    'Prepared MAT file does not contain the data struct.');
data = loaded.data;
assert(iscell(data.X) && numel(data.X)==5,'Expected five signal views.');
assert(isequal(double(data.horizons(:).'),[1 20 45 70 120]), ...
    'Unexpected return horizons.');
for k=1:5
    assert(isequal(size(data.X{k}),[100 178]), ...
        'Unexpected signal dimensions in view %d.',k);
    assert(all(isfinite(data.X{k}(:))),'Non-finite values in view %d.',k);
end
assert(numel(data.tickers)==100,'Unexpected stock count.');
assert(numel(data.common_dates)==178,'Unexpected coordinate count.');

selected = readtable(fullfile(root,'selected_top1', ...
    'SELECTED_TOP1_CONFIGURATION.csv'),'VariableNamingRule','preserve');
assert(height(selected)==1,'Expected one reported configuration.');
assert(selected.r(1)==3 && selected.profile_id(1)==124, ...
    'Reported rank/profile does not match the selected model.');

fprintf('Package verification passed.\n');
