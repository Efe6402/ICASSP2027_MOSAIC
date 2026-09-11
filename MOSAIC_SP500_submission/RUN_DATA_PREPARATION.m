%% Construct the offset-0 S&P 500 multihorizon signal dataset
root = fileparts(mfilename('fullpath'));
run(fullfile(root,'data','build_mosaic_finance_multihorizon_dataset.m'));
