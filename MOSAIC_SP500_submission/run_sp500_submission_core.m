function results = run_sp500_submission_core(userCfg)
%RUN_SP500_SUBMISSION_CORE End-to-end MOSAIC finance experiment.
% -------------------------------------------------------------------------
% Single-file end-to-end MOSAIC finance experiment.
%
% WHAT THIS FILE DOES, SEQUENTIALLY
%   1) Builds five S&P multi-horizon datasets for 100 stocks, return
%      horizons [1 20 45 70 120], and the five endpoint residue classes
%      modulo a 5-trading-day endpoint step.
%
%      dataset 1 : endpoint offset 0  -> final estimation target
%      dataset 2 : endpoint offset 1  -> estimation configuration A
%      dataset 3 : endpoint offset 2  -> validation configuration A
%      dataset 4 : endpoint offset 3  -> estimation configuration B
%      dataset 5 : endpoint offset 4  -> validation configuration B
%
%      Thus the model-selection transfer design is
%             dataset 2  -> dataset 3
%             dataset 4  -> dataset 5
%      and dataset 1 is never used to choose hyperparameters or rank.
%
%   2) Builds a 128-profile hyperparameter bank comprising 64 finance
%      profiles and a focused LHS neighborhood around profile 59 (lhs_51).
%
%   3) For every candidate modal rank r and every hyperparameter profile:
%        - fit MOSAIC independently on dataset 2 and dataset 4;
%        - mass-normalize each learned supra graph;
%        - score its deterministic smoothness on its validation dataset
%          (2->3 and 4->5) WITHOUT refitting;
%        - average the two held-out transfer scores.
%
%   4) At every r, select the profile with the lowest mean validation-configuration
%      smoothness V*(r).  Select r* by a geometric endpoint-chord elbow and
%      retain the lowest-scoring hyperparameter profile at r*.
%
%   5) Open dataset 1 (offset 0) only after selection. Fit the selected
%      profile from scratch at r*. The target does not alter the selection.
%
%   6) Save finance, cross-vs-within, modal, and view-graph analyses for the
%      selected target fit, including:
%        - convergence, modal B profiles, A^(m), within-horizon adjacencies,
%          supra overview, horizon dependence, strongest cross/copy edges,
%          sector diagnostics, stock/sector strengths;
%        - induced view graph + spectral clustering diagnostics;
%        - cross-vs-within edge comparisons and rank percentiles;
%        - complete pair-specific KxK horizon-coupling matrices;
%        - complete within-view and mode adjacency CSVs;
%        - mode decomposition of strongest relations.
%
% IMPORTANT INTERPRETATION
%   The five offset datasets are distinct endpoint-sampling configurations,
%   not statistically independent samples.  This pipeline uses them only as
%   deterministic local offset-validation configurations for cross-configuration
%   smoothness selection.  No probabilistic iid/Gaussian claim is required.
%
% DEFAULT RANKS
%   r = 2:5.  K=5, so these provide four points for an elbow while covering
%   the previously interesting r=3 and r=4 regimes.
%
% OUTPUT
%   A timestamped folder is created beside this file unless cfg.output_root
%   is supplied.
%
% USAGE
%   results = run_sp500_submission_core;
%
% Optional overrides:
%   cfg = struct();
%   cfg.rank_candidates = 2:5;
%   cfg.local_neighbor_count = 64;     % 64 reference + 64 local = 128
%   cfg.selection_n_restarts = 2;
%   cfg.selection_max_iter = 400;
%   cfg.final_n_restarts = 5;
%   cfg.final_max_iter = 800;
%   results = run_sp500_submission_core(cfg);
%
% Resume an interrupted run by reusing the same output root:
%   cfg.output_root = '/absolute/path/MOSAIC_FINANCE_OFFSET_VALIDATION_RANK_SELECTION_...';
%   cfg.resume = true;
%   results = run_sp500_submission_core(cfg);
%
% August 2026
% -------------------------------------------------------------------------

    if nargin < 1 || isempty(userCfg), userCfg = struct(); end
    cfg = cf_defaults(userCfg);

    scriptDir = fileparts(mfilename('fullpath'));
    if isempty(scriptDir), scriptDir = pwd; end

    % Resolve raw-data paths only after scriptDir is known.
    if isempty(cfg.stocks_dir), cfg.stocks_dir = fullfile(scriptDir,'stocks'); end
    if isempty(cfg.price_csv), cfg.price_csv = fullfile(cfg.stocks_dir,'sp500-data-2016-2020_converted.csv'); end
    if isempty(cfg.sector_csv), cfg.sector_csv = fullfile(cfg.stocks_dir,'SP500-sectors.csv'); end
    if isempty(cfg.frozen_ticker_csv), cfg.frozen_ticker_csv = fullfile(cfg.stocks_dir,'mosaic_selected_tickers.csv'); end
    if isempty(cfg.prepared_target_mat)
        candidate = fullfile(scriptDir,'data','prepared_offset0', ...
            'mosaic_finance_multihorizon_dataset.mat');
        if isfile(candidate), cfg.prepared_target_mat = candidate; end
    end

    if isempty(cfg.output_root)
        stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
        outputRoot = fullfile(scriptDir,['MOSAIC_FINANCE_OFFSET_VALIDATION_RANK_SELECTION_',stamp]);
    else
        outputRoot = char(cfg.output_root);
    end
    cfg.output_root = outputRoot;

    cf_ensure_folder(outputRoot);
    datasetRoot = fullfile(outputRoot,'00_OFFSET_DATASETS');
    selectionRoot = fullfile(outputRoot,'01_OFFSET_VALIDATION_RANK_HYPERPARAMETER_SELECTION');
    finalRoot = fullfile(outputRoot,'02_FINAL_OFFSET0_TARGET_TOP1');
    summaryRoot = fullfile(outputRoot,'03_MASTER_SUMMARY');
    cf_ensure_folder(datasetRoot); cf_ensure_folder(selectionRoot);
    cf_ensure_folder(finalRoot); cf_ensure_folder(summaryRoot);

    oldVisibility = get(groot,'DefaultFigureVisible');
    set(groot,'DefaultFigureVisible','off');
    visibilityCleanup = onCleanup(@() set(groot,'DefaultFigureVisible',oldVisibility)); %#ok<NASGU>

    diary('off');
    logFile = fullfile(outputRoot,'MASTER_RUN_LOG.txt');
    diary(logFile);
    diaryCleanup = onCleanup(@() diary('off')); %#ok<NASGU>

    fprintf('\n============================================================\n');
    fprintf('MOSAIC FINANCE OFFSET-VALIDATION / RANK-SELECTION FULL PIPELINE\n');
    fprintf('============================================================\n');
    fprintf('Output root : %s\n',outputRoot);
    fprintf('Raw prices  : %s\n',cfg.price_csv);
    fprintf('Horizons    : %s\n',mat2str(cfg.horizons));
    fprintf('Ranks       : %s\n',mat2str(cfg.rank_candidates));
    fprintf('Offsets     : %s\n',mat2str(cfg.offsets));
    fprintf('Transfers   : dataset 2->3 and dataset 4->5\n');
    fprintf('Final target: dataset 1 (offset 0)\n\n');

    % ------------------------------------------------------------------
    % 1) DATASETS
    % ------------------------------------------------------------------
    allDataFile = fullfile(datasetRoot,'ALL_FIVE_OFFSET_DATASETS.mat');
    if cfg.resume && isfile(allDataFile)
        fprintf('[DATA] Loading previously generated five-offset datasets...\n');
        Z = load(allDataFile,'dataSets','marketInfo','datasetDesign');
        dataSets = Z.dataSets; marketInfo = Z.marketInfo; datasetDesign = Z.datasetDesign;
    else
        fprintf('[DATA] Building five endpoint-offset configurations from raw prices...\n');
        [dataSets,marketInfo,datasetDesign] = cf_build_five_offset_datasets(cfg,datasetRoot);
        save(allDataFile,'dataSets','marketInfo','datasetDesign','cfg','-v7.3');
    end
    if ~isempty(cfg.prepared_target_mat)
        dataSets{1} = cf_load_and_verify_prepared_target( ...
            cfg.prepared_target_mat,dataSets{1},cfg);
    end
    writetable(datasetDesign,fullfile(datasetRoot,'OFFSET_DATASET_DESIGN.csv'));
    cf_verify_dataset_alignment(dataSets,cfg);
    cf_save_offset_overlap_table(dataSets,datasetRoot);

    % Prepare normalized geometry and diagnostics independently by offset.
    prepared = cell(5,1);
    for dId = 1:5
        diagDir = fullfile(datasetRoot,sprintf('dataset_%d_offset_%d',dId,dataSets{dId}.endpoint_offset), ...
            'input_diagnostics');
        cf_ensure_folder(diagDir);
        [dataSets{dId},Yd,obsd,inputInfo] = finance_prepare_observed(dataSets{dId},cfg);
        prepared{dId} = struct('Y',Yd,'observed',obsd,'inputInfo',inputInfo);
        finance_save_input_diagnostics(dataSets{dId},Yd,inputInfo,diagDir,cfg);
        data = dataSets{dId}; %#ok<NASGU>
        save(fullfile(datasetRoot,sprintf('dataset_%d_offset_%d',dId,data.endpoint_offset), ...
            'mosaic_finance_multihorizon_dataset.mat'),'data','-v7.3');
    end
    save(allDataFile,'dataSets','marketInfo','datasetDesign','cfg','-v7.3');

    % ------------------------------------------------------------------
    % 2) HYPERPARAMETER BANK: REFERENCE 64 + HP59-LOCAL NEIGHBORHOOD
    % ------------------------------------------------------------------
    profiles = cf_build_profile_bank(cfg);
    profileTable = cf_profiles_to_table(profiles);
    writetable(profileTable,fullfile(selectionRoot,'HYPERPARAMETER_BANK_REFERENCE64_PLUS_HP59_NEIGHBORS.csv'));
    hp59 = profiles([profiles.profile_id]==59);
    writetable(cf_profiles_to_table(hp59),fullfile(selectionRoot,'HP59_LHS51_EXACT_ANCHOR.csv'));
    save(fullfile(selectionRoot,'HYPERPARAMETER_BANK.mat'),'profiles','profileTable','hp59','cfg','-v7.3');
    fprintf('[HP] Bank size: %d profiles per rank (%d reference + %d local).\n', ...
        numel(profiles),64,cfg.local_neighbor_count);

    % ------------------------------------------------------------------
    % 3-4) VALIDATION TRANSFER SELECTION + RANK ELBOW
    % ------------------------------------------------------------------
    selectionFile = fullfile(selectionRoot,'OFFSET_VALIDATION_SELECTION_RESULT.mat');
    if strcmpi(cfg.stage,'final')
        if ~isfile(selectionFile)
            error('cfg.stage=final requires %s',selectionFile);
        end
        Z = load(selectionFile,'selection'); selection = Z.selection;
    else
        selection = cf_run_offset_validation_selection(dataSets,prepared,profiles,cfg,selectionRoot);
        save(selectionFile,'selection','profiles','profileTable','cfg','-v7.3');
    end

    if strcmpi(cfg.stage,'selection')
        results = struct('output_root',outputRoot,'datasets',{dataSets}, ...
            'selection',selection,'configuration',cfg);
        save(fullfile(summaryRoot,'MASTER_RESULTS_SELECTION_ONLY.mat'),'results','-v7.3');
        cf_write_master_readme(fullfile(outputRoot,'README_PIPELINE.txt'),cfg,selection,[]);
        fprintf('\nSelection-only stage complete.\n');
        return;
    end

    % ------------------------------------------------------------------
    % 5-6) FINAL OFFSET-0 TARGET FITS + FULL ANALYSIS/POSTPROCESSING
    % ------------------------------------------------------------------
    final = cf_run_final_target(dataSets{1},prepared{1},profiles,selection,cfg,finalRoot);

    results = struct();
    results.output_root = outputRoot;
    results.configuration = cfg;
    results.dataset_design = datasetDesign;
    results.selection = selection;
    results.final = final;
    results.primary_profile_id = selection.top1_table.profile_id(1);
    results.primary_profile_name = selection.top1_table.profile_name{1};
    results.selected_rank = selection.selected_rank;
    results.target_used_for_reselection = false;
    save(fullfile(summaryRoot,'MASTER_RESULTS.mat'),'results','-v7.3');

    writetable(selection.rank_table,fullfile(summaryRoot,'FINAL_RANK_SELECTION_SUMMARY.csv'));
    writetable(selection.top1_table,fullfile(summaryRoot,'SELECTED_TOP1_CONFIGURATION.csv'));
    writetable(final.summary_table,fullfile(summaryRoot,'FINAL_TOP1_TARGET_FITS.csv'));
    cf_write_master_readme(fullfile(outputRoot,'README_PIPELINE.txt'),cfg,selection,final);

    fprintf('\n============================================================\n');
    fprintf('FULL PIPELINE COMPLETED\n');
    fprintf('============================================================\n');
    fprintf('Selected rank r* : %d\n',selection.selected_rank);
    fprintf('Primary profile  : %d %s\n',selection.top1_table.profile_id(1), ...
        selection.top1_table.profile_name{1});
    fprintf('Output root      : %s\n\n',outputRoot);
end


%% ========================================================================
%% PIPELINE DEFAULTS
%% ========================================================================
function cfg = cf_defaults(userCfg)
    cfg = struct();
    cfg.seed = 20260813;
    cfg.stage = 'all';                    % 'all', 'selection', 'final'
    cfg.resume = true;
    cfg.output_root = '';

    cfg.stocks_dir = '';
    cfg.price_csv = '';
    cfg.sector_csv = '';
    cfg.frozen_ticker_csv = '';
    cfg.prepared_target_mat = '';
    cfg.num_stocks = 100;
    cfg.horizons = [1 20 45 70 120];
    cfg.endpoint_step = 5;
    cfg.offsets = 0:4;

    % Dataset IDs are one-based: offset 0 is dataset 1.
    cfg.transfer_pairs = [2 3; 4 5];
    cfg.target_dataset_id = 1;

    % Rank search.  K=5, and r=2:5 provides enough points for an elbow.
    cfg.rank_candidates = 2:5;
    cfg.top_k_profiles = 1;

    % Profile bank: 64 finance profiles plus an hp59/lhs51 neighborhood.
    cfg.local_neighbor_count = 64;
    cfg.local_alpha_multiplier = [0.55 1.80];
    cfg.local_beta_multiplier  = [0.60 1.70];
    cfg.local_tau_multiplier   = [0.55 1.80];
    cfg.local_eta_multiplier   = [0.60 1.65];
    cfg.local_lambda_multiplier= [0.45 2.00];
    cfg.local_lambda_clip      = [0.001 0.20];

    % Selection solver settings.  Two restarts mirrors the focused COIL run.
    cfg.selection_n_restarts = 2;
    cfg.selection_max_iter = 400;
    cfg.final_n_restarts = 5;
    cfg.final_max_iter = 800;
    cfg.stationarity_tol = 1e-6;
    cfg.objective_tol = 1e-9;
    cfg.solver_verbose = false;

    % Prepared signals are already standardized separately per stock/horizon.
    cfg.center_each_signal = false;
    cfg.l2_normalize_each_signal = false;

    % Finance diagnostics.
    cfg.display_relative_threshold = 0.05;
    cfg.top_within_edges_per_view = 20;
    cfg.top_cross_edges = 30;
    cfg.top_copy_edges = 25;
    cfg.top_stock_strengths = 25;
    cfg.max_cross_pair_repetitions = 5;
    cfg.max_stock_pair_repetitions = 3;
    cfg.final_top_relations = 40;
    cfg.figure_resolution = 220;
    cfg.save_detailed_figures = true;

    % View-graph clustering inherited from the stability code.
    cfg.view_cluster_counts = 2:4;
    cfg.view_kmeans_restarts = 100;
    cfg.view_kmeans_max_iter = 500;
    cfg.max_eigengap_cluster_count = 4;

    % Cross-vs-within / mode postprocessing inherited from v3.
    cfg.top_cross_relations = 15;
    cfg.top_within_edges_per_horizon = 25;
    cfg.top_mode_edges = 25;
    cfg.save_individual_pair_figures = false;
    cfg.max_pair_figures = 15;
    cfg.matrix_csv_precision = '%.16g';

    cfg = cf_merge(cfg,userCfg);
    cfg.rank_candidates = unique(double(cfg.rank_candidates(:).'));
    cfg.offsets = double(cfg.offsets(:).');
    cfg.view_cluster_counts = unique(double(cfg.view_cluster_counts(:).'));

    if ~isequal(cfg.offsets,0:4)
        error('The offset-validation design requires cfg.offsets = 0:4.');
    end
    if size(cfg.transfer_pairs,2)~=2
        error('cfg.transfer_pairs must have two columns [source validation].');
    end
    if numel(cfg.rank_candidates)<3
        error('At least three candidate ranks are required for elbow selection.');
    end
    validStages = {'all','selection','final'};
    if ~any(strcmpi(cfg.stage,validStages))
        error('cfg.stage must be all, selection, or final.');
    end
end


%% ========================================================================
%% FIVE OFFSET DATASETS
%% ========================================================================
function [dataSets,marketInfo,designT] = cf_build_five_offset_datasets(cfg,outRoot)
    assert(isfolder(cfg.stocks_dir),'Cannot find stocks folder: %s',cfg.stocks_dir);
    assert(isfile(cfg.price_csv),'Cannot find price CSV: %s',cfg.price_csv);
    assert(isfile(cfg.sector_csv),'Cannot find sector CSV: %s',cfg.sector_csv);

    Traw = readtable(cfg.price_csv,'VariableNamingRule','preserve');
    rawDate = Traw{:,1};
    priceDates = parse_dates_strict(rawDate);
    priceTable = Traw(:,2:end);
    [Praw,numericKeep] = table_to_numeric_matrix_with_mask(priceTable);
    tickersAll = string(priceTable.Properties.VariableNames);
    tickersAll = tickersAll(numericKeep);

    goodDate = ~isnat(priceDates);
    priceDates = priceDates(goodDate); Praw = Praw(goodDate,:);
    [priceDates,ord] = sort(priceDates,'ascend'); Praw = Praw(ord,:);
    [priceDatesUnique,uidx] = unique(priceDates,'stable');
    priceDates = priceDatesUnique; Praw = Praw(uidx,:);

    validStock = all(isfinite(Praw),1) & all(Praw>0,1);
    if isfile(cfg.frozen_ticker_csv)
        F = readtable(cfg.frozen_ticker_csv,'VariableNamingRule','preserve');
        frozenTickers = string(F{:,1});
        frozenTickers = frozenTickers(strlength(strtrim(frozenTickers))>0);
        [found,selectedIdx] = match_tickers(frozenTickers,tickersAll);
        if any(~found), error('Frozen ticker list contains tickers absent from price CSV.'); end
        if any(~validStock(selectedIdx)), error('Frozen ticker list contains invalid/incomplete price series.'); end
        if numel(selectedIdx)~=cfg.num_stocks
            error('Frozen ticker count %d does not match cfg.num_stocks=%d.',numel(selectedIdx),cfg.num_stocks);
        end
    else
        validIdx = find(validStock);
        if numel(validIdx)<cfg.num_stocks
            error('Only %d complete positive-price stocks are available.',numel(validIdx));
        end
        selectedIdx = validIdx(1:cfg.num_stocks);
        frozenTickers = tickersAll(selectedIdx);
        writetable(table(frozenTickers(:),'VariableNames',{'Ticker'}),cfg.frozen_ticker_csv);
    end

    P = double(Praw(:,selectedIdx));
    tickers = tickersAll(selectedIdx);
    n = size(P,2);
    [sectorLabels,sectorKnown] = load_sector_labels(cfg.sector_csv,tickers);
    [sectorOrder,sectorNames,sectorBoundaries] = make_sector_order(sectorLabels);
    logP = log(P);

    sharedDir = fullfile(outRoot,'shared_market_input'); cf_ensure_folder(sharedDir);
    stockMeta = table((1:n)',tickers(:),sectorLabels(:),sectorKnown(:), ...
        'VariableNames',{'NodeIndex','Ticker','GICSSector','SectorKnown'});
    writetable(stockMeta,fullfile(sharedDir,'selected_stocks_and_sectors.csv'));
    save(fullfile(sharedDir,'shared_prices_and_metadata.mat'),'P','logP','priceDates', ...
        'tickers','sectorLabels','sectorKnown','sectorOrder','sectorNames', ...
        'sectorBoundaries','selectedIdx','-v7.3');

    K = numel(cfg.horizons); Hmax = max(cfg.horizons); Ndates = size(P,1);
    dataSets = cell(5,1);
    rows = repmat(struct('dataset_id',NaN,'offset',NaN,'role','', ...
        'source_for_dataset',NaN,'validation_for_dataset',NaN,'n',NaN,'K',NaN, ...
        'p',NaN,'first_ending_date',NaT,'last_ending_date',NaT),5,1);

    for offset = 0:4
        dId = offset+1;
        startIdx = Hmax+1+offset;
        endIdx = startIdx:cfg.endpoint_step:Ndates;
        if isempty(endIdx), error('Offset %d leaves no valid endpoints.',offset); end
        commonDates = priceDates(endIdx);
        p = numel(endIdx);

        X = cell(1,K); Xraw = cell(1,K);
        normMu = zeros(n,K); normSigma = zeros(n,K);
        viewNames = cell(1,K);
        for k=1:K
            h = cfg.horizons(k);
            Rh = (logP(endIdx,:) - logP(endIdx-h,:)).';
            mu = mean(Rh,2); sig = std(Rh,0,2);
            bad = sig<1e-12 | ~isfinite(sig); sig(bad)=1;
            Xh = (Rh-mu)./sig;
            if any(~isfinite(Xh(:))), error('Nonfinite standardized values at offset=%d h=%d.',offset,h); end
            X{k}=Xh; Xraw{k}=Rh; normMu(:,k)=mu; normSigma(:,k)=sig;
            viewNames{k}=sprintf('%d-day',h);
        end

        data = struct();
        data.description = sprintf(['S&P 500 MOSAIC endpoint-offset dataset. ', ...
            'Offset=%d modulo a five-trading-day endpoint grid.'],offset);
        data.X=X; data.X_raw=Xraw; data.horizons=cfg.horizons(:).';
        data.view_names=viewNames; data.view_labels=string(viewNames(:));
        data.K=K; data.n=n; data.p=p; data.num_views=K; data.num_stocks=n;
        data.num_objects=n; data.signal_dimension=p;
        data.endpoint_step=cfg.endpoint_step; data.endpoint_offset=offset;
        data.dataset_id=dId; data.common_dates=commonDates;
        data.end_price_indices=endIdx; data.tickers=tickers(:); data.node_names=tickers(:);
        data.sector_labels=sectorLabels(:); data.sector_known=sectorKnown(:);
        data.sector_order=sectorOrder; data.sector_names=sectorNames;
        data.sector_boundaries=sectorBoundaries;
        data.normalization='row_zscore_within_each_offset_and_horizon';
        data.norm_mu=normMu; data.norm_sigma=normSigma;
        data.preprocessing_notes = { ...
            'No prebuilt correlation/adjacency is supplied to MOSAIC.'; ...
            'All five offset datasets use the same stocks and horizon definitions.'; ...
            'Each offset uses a disjoint residue class of endpoint dates modulo 5.'; ...
            'Rolling-return windows can overlap across offsets; offsets are offset-validation configurations, not iid samples.'; ...
            'Each stock is z-scored separately within each offset and horizon.'};
        dataSets{dId}=data;

        dDir = fullfile(outRoot,sprintf('dataset_%d_offset_%d',dId,offset));
        cf_ensure_folder(dDir);
        writetable(table((1:p)',endIdx(:),commonDates(:), ...
            'VariableNames',{'CoordinateIndex','PriceRowIndex','EndingDate'}), ...
            fullfile(dDir,'common_ending_dates.csv'));

        rows(dId).dataset_id=dId; rows(dId).offset=offset; rows(dId).n=n;
        rows(dId).K=K; rows(dId).p=p; rows(dId).first_ending_date=commonDates(1);
        rows(dId).last_ending_date=commonDates(end);
        if dId==1
            rows(dId).role='FINAL_TARGET_UNTOUCHED_FOR_SELECTION';
        elseif dId==2
            rows(dId).role='ESTIMATION_CONFIGURATION_A'; rows(dId).validation_for_dataset=3;
        elseif dId==3
            rows(dId).role='VALIDATION_CONFIGURATION_A'; rows(dId).source_for_dataset=2;
        elseif dId==4
            rows(dId).role='ESTIMATION_CONFIGURATION_B'; rows(dId).validation_for_dataset=5;
        else
            rows(dId).role='VALIDATION_CONFIGURATION_B'; rows(dId).source_for_dataset=4;
        end
    end

    designT = struct2table(rows,'AsArray',true);
    marketInfo = struct('price_csv',cfg.price_csv,'sector_csv',cfg.sector_csv, ...
        'frozen_ticker_csv',cfg.frozen_ticker_csv,'price_dates',priceDates, ...
        'tickers',tickers,'n',n,'K',K,'horizons',cfg.horizons);
end


function cf_verify_dataset_alignment(dataSets,cfg)
    ref = dataSets{1};
    for q=1:numel(dataSets)
        d=dataSets{q};
        if d.num_stocks~=ref.num_stocks || d.num_views~=ref.num_views
            error('Offset dataset %d has incompatible n/K.',q);
        end
        if ~isequal(string(d.tickers(:)),string(ref.tickers(:)))
            error('Offset dataset %d does not preserve stock identity/order.',q);
        end
        if ~isequal(double(d.horizons(:).'),double(ref.horizons(:).'))
            error('Offset dataset %d does not preserve horizon views.',q);
        end
    end
    if ref.num_stocks~=cfg.num_stocks, error('Unexpected stock count.'); end
end


function data = cf_load_and_verify_prepared_target(matPath,generated,cfg)
    assert(isfile(matPath),'Cannot find prepared offset-0 dataset: %s',matPath);
    loaded = load(matPath,'data');
    assert(isfield(loaded,'data') && isstruct(loaded.data), ...
        'The prepared dataset must contain a struct named data.');
    data = loaded.data;

    assert(isfield(data,'X') && iscell(data.X), ...
        'The prepared dataset does not contain the required view signals.');
    assert(numel(data.X)==numel(generated.X), ...
        'Prepared and generated datasets have different view counts.');
    for k=1:numel(data.X)
        assert(isequal(size(data.X{k}),size(generated.X{k})), ...
            'Prepared and generated signals differ in size for view %d.',k);
        discrepancy=max(abs(double(data.X{k}(:))-double(generated.X{k}(:))));
        assert(discrepancy<=1e-12, ...
            'Prepared and generated signals differ for view %d (max error %.3g).', ...
            k,discrepancy);
    end

    assert(isequal(double(data.horizons(:).'),double(cfg.horizons(:).')), ...
        'Prepared dataset horizons do not match the experiment configuration.');
    assert(isequal(string(data.tickers(:)),string(generated.tickers(:))), ...
        'Prepared and generated datasets use different stocks or row ordering.');
    assert(isequal(data.common_dates(:),generated.common_dates(:)), ...
        'Prepared and generated datasets use different ending dates.');

    data.dataset_id=1;
    data.endpoint_offset=0;
    data.endpoint_step=cfg.endpoint_step;
    data.num_stocks=size(data.X{1},1);
    data.num_views=numel(data.X);
    data.num_samples=size(data.X{1},2);
    data.num_objects=data.num_stocks;
    data.signal_dimension=data.num_samples;
    fprintf('[DATA] Verified and loaded prepared offset-0 target: %s\n',matPath);
end


function cf_save_offset_overlap_table(dataSets,outRoot)
    m = numel(dataSets); M=zeros(m);
    for i=1:m
        for j=1:m
            M(i,j)=numel(intersect(dataSets{i}.end_price_indices,dataSets{j}.end_price_indices));
        end
    end
    names = matlab.lang.makeValidName("dataset"+string((1:m)')+"_offset"+string((0:m-1)'));
    T=array2table(M,'VariableNames',cellstr(names),'RowNames',cellstr(names));
    writetable(T,fullfile(outRoot,'ENDPOINT_INDEX_OVERLAP_COUNTS.csv'),'WriteRowNames',true);
end


%% ========================================================================
%% HYPERPARAMETER BANK
%% ========================================================================
function profiles = cf_build_profile_bank(cfg)
    % Generate the 64-profile finance bank; profile 59 is lhs_51.
    b = struct();
    b.seed = 20260801;
    b.n_profiles_per_rank = 64;
    b.alpha_base_range = [0.0015 0.0120];
    b.alpha_ratio_range = [0.35 2.85];
    b.beta_base_range = [0.035 0.180];
    b.beta_ratio_range = [0.50 2.00];
    b.tau_base_range = [0.040 0.220];
    b.tau_inter_intra_ratio_range = [0.35 2.85];
    b.eta_c_range = [0.40 2.50];
    b.lambda_B_range = [0.005 0.120];
    referenceProfiles = stability_generate_profiles(b);
    for q=1:numel(referenceProfiles), referenceProfiles(q).origin='reference64'; end

    anchor = referenceProfiles(59);
    nLocal = cfg.local_neighbor_count;
    local = repmat(anchor,nLocal,1);
    U = cf_latin_hypercube(nLocal,8,cfg.seed+590051);
    for q=1:nLocal
        p=anchor;
        p.profile_id = 64+q;
        p.name = sprintf('hp59_neighbor_%02d',q);
        p.slug = sprintf('hp59nbr%02d',q);
        p.origin = 'hp59_local_neighborhood';
        p.alpha_a = anchor.alpha_a*cf_log_multiplier(U(q,1),cfg.local_alpha_multiplier);
        p.alpha_c = anchor.alpha_c*cf_log_multiplier(U(q,2),cfg.local_alpha_multiplier);
        p.beta_a  = anchor.beta_a *cf_log_multiplier(U(q,3),cfg.local_beta_multiplier);
        p.beta_c  = anchor.beta_c *cf_log_multiplier(U(q,4),cfg.local_beta_multiplier);
        p.tau_intra = anchor.tau_intra*cf_log_multiplier(U(q,5),cfg.local_tau_multiplier);
        p.tau_inter = anchor.tau_inter*cf_log_multiplier(U(q,6),cfg.local_tau_multiplier);
        p.eta_c = anchor.eta_c*cf_log_multiplier(U(q,7),cfg.local_eta_multiplier);
        p.lambda_B = anchor.lambda_B*cf_log_multiplier(U(q,8),cfg.local_lambda_multiplier);
        p.lambda_B = min(max(p.lambda_B,cfg.local_lambda_clip(1)),cfg.local_lambda_clip(2));
        p.family = stability_classify_profile(p.alpha_a,p.alpha_c,p.beta_a,p.beta_c, ...
            p.tau_intra,p.tau_inter,p.eta_c,p.lambda_B);
        local(q)=p;
    end
    profiles=[referenceProfiles(:);local(:)];
end


function U=cf_latin_hypercube(n,d,seed)
    old=rng; cleanup=onCleanup(@()rng(old)); %#ok<NASGU>
    rng(seed,'twister'); U=zeros(n,d);
    for j=1:d
        perm=randperm(n); jitter=rand(n,1); U(:,j)=(perm(:)-jitter)/n;
    end
    U=min(max(U,eps),1-eps);
end


function m=cf_log_multiplier(u,bounds)
    m=exp(log(bounds(1))+u*(log(bounds(2))-log(bounds(1))));
end


function T=cf_profiles_to_table(profiles)
    if isempty(profiles), T=table(); return; end
    T=struct2table(profiles);
    T.alpha_a_over_alpha_c=T.alpha_a./T.alpha_c;
    T.beta_a_over_beta_c=T.beta_a./T.beta_c;
    T.tau_inter_over_tau_intra=T.tau_inter./T.tau_intra;
end


function opts=cf_profile_to_options(profile,cfg,seed,stage)
    opts=struct(); opts.seed=seed;
    if strcmpi(stage,'selection')
        opts.n_restarts=cfg.selection_n_restarts; opts.max_iter=cfg.selection_max_iter;
    elseif strcmpi(stage,'final')
        opts.n_restarts=cfg.final_n_restarts; opts.max_iter=cfg.final_max_iter;
    else
        error('Unknown fitting stage: %s',stage);
    end
    opts.stationarity_tol=cfg.stationarity_tol;
    opts.objective_tol=cfg.objective_tol;
    opts.verbose=cfg.solver_verbose;
    opts.alpha_a=profile.alpha_a; opts.alpha_c=profile.alpha_c;
    opts.beta_a=profile.beta_a; opts.beta_c=profile.beta_c;
    opts.tau_intra=profile.tau_intra; opts.tau_inter=profile.tau_inter;
    opts.eta_c=profile.eta_c; opts.lambda_B=profile.lambda_B;
    opts.epsilon_intra=profile.epsilon_intra;
    opts.epsilon_inter=profile.epsilon_inter;
end


%% ========================================================================
%% OFFSET CROSS-CONFIGURATION SELECTION
%% ========================================================================
function selection=cf_run_offset_validation_selection(dataSets,prepared,profiles,cfg,outRoot)
    ranks=cfg.rank_candidates; pairs=cfg.transfer_pairs;
    nPairs=size(pairs,1); nProfiles=numel(profiles);
    fitRoot=fullfile(outRoot,'candidate_records'); cf_ensure_folder(fitRoot);

    total=nPairs*numel(ranks)*nProfiles;
    rows=repmat(cf_empty_transfer_row(),total,1); cursor=0;

    for pp=1:nPairs
        srcId=pairs(pp,1); valId=pairs(pp,2);
        src=dataSets{srcId}; val=dataSets{valId};
        sourceDir=fullfile(fitRoot,sprintf('dataset_%d_offset_%d_to_dataset_%d_offset_%d', ...
            srcId,src.endpoint_offset,valId,val.endpoint_offset));
        cf_ensure_folder(sourceDir);
        fprintf('\n[TRANSFER %d/%d] dataset %d(offset %d) -> dataset %d(offset %d)\n', ...
            pp,nPairs,srcId,src.endpoint_offset,valId,val.endpoint_offset);

        for ir=1:numel(ranks)
            r=ranks(ir); rankDir=fullfile(sourceDir,sprintf('r_%02d',r)); cf_ensure_folder(rankDir);
            pairedSeed=cfg.seed+1000000*srcId+10000*r;
            observed=prepared{srcId}.observed; observed.r=r;
            Sval=prepared{valId}.observed.S;

            for ip=1:nProfiles
                profile=profiles(ip); cursor=cursor+1;
                recFile=fullfile(rankDir,sprintf('profile_%03d_%s.mat',profile.profile_id,cf_safe(profile.slug)));
                rec=cf_empty_transfer_row(); reuse=false;
                if cfg.resume && isfile(recFile)
                    Z=load(recFile,'record');
                    if isfield(Z,'record') && isfield(Z.record,'completed') && Z.record.completed
                        rec=Z.record; reuse=true;
                    end
                end
                if ~reuse
                    rec.source_dataset_id=srcId; rec.validation_dataset_id=valId;
                    rec.source_offset=src.endpoint_offset; rec.validation_offset=val.endpoint_offset;
                    rec.r=r; rec.profile_id=profile.profile_id; rec.profile_name=profile.name;
                    rec.profile_origin=profile.origin; rec.profile_family=profile.family;
                    rec.alpha_a=profile.alpha_a; rec.alpha_c=profile.alpha_c;
                    rec.beta_a=profile.beta_a; rec.beta_c=profile.beta_c;
                    rec.tau_intra=profile.tau_intra; rec.tau_inter=profile.tau_inter;
                    rec.eta_c=profile.eta_c; rec.lambda_B=profile.lambda_B;
                    rec.seed=pairedSeed;
                    opts=cf_profile_to_options(profile,cfg,pairedSeed,'selection');
                    fprintf('  r=%d profile %3d/%3d %-22s ',r,ip,nProfiles,profile.slug);
                    timer=tic;
                    try
                        fit=solve_mosaic_crossview_embedded(observed,opts);
                        rec.runtime_seconds=toc(timer); rec.status=fit.status;
                        rec.objective=fit.objective; rec.stationarity=fit.stationarity_combined;
                        rec.iterations=fit.iterations; rec.best_restart=fit.best_restart;
                        vm=cf_graph_validation_metrics(fit.A_supra,Sval,src.num_stocks,src.num_views);
                        gd=cf_graph_diagnostics(fit,src.num_stocks,src.num_views);
                        rec.validation_score=vm.total_score;
                        rec.validation_within_score=vm.within_score;
                        rec.validation_cross_score=vm.cross_score;
                        rec.validation_copy_score=vm.copy_score;
                        rec.within_mass_fraction=vm.within_mass_fraction;
                        rec.cross_mass_fraction=vm.cross_mass_fraction;
                        rec.copy_mass_fraction=vm.copy_mass_fraction;
                        rec.graph_total_mass=vm.total_mass;
                        rec.max_B_cosine=gd.max_B_cosine;
                        rec.min_modal_mass_fraction=gd.min_modal_mass_fraction;
                        rec.valid=isfinite(vm.total_score) && vm.total_mass>0;
                        rec.completed=true;
                        fprintf('V=%.6g residual=%.2e\n',rec.validation_score,rec.stationarity);
                    catch ME
                        rec.runtime_seconds=toc(timer); rec.status=['error:' ME.identifier];
                        rec.error_message=ME.message; rec.valid=false; rec.completed=true;
                        fprintf('FAILED: %s\n',ME.message);
                    end
                    record=rec; %#ok<NASGU>
                    save(recFile,'record','profile','opts');
                end
                rows(cursor)=rec;
            end
        end
        writetable(struct2table(rows(1:cursor),'AsArray',true), ...
            fullfile(outRoot,'selection_transfer_scores_checkpoint.csv'));
    end

    rows=rows(1:cursor); transferT=struct2table(rows,'AsArray',true);
    writetable(transferT,fullfile(outRoot,'selection_transfer_scores_long.csv'));

    aggRows=repmat(cf_empty_aggregate_row(),numel(ranks)*nProfiles,1); c=0;
    for ir=1:numel(ranks)
        r=ranks(ir);
        for ip=1:nProfiles
            profile=profiles(ip); c=c+1;
            sub=transferT(transferT.r==r & transferT.profile_id==profile.profile_id,:);
            scores=sub.validation_score(sub.valid & isfinite(sub.validation_score));
            a=cf_empty_aggregate_row(); a.r=r; a.profile_id=profile.profile_id;
            a.profile_name=profile.name; a.profile_origin=profile.origin; a.profile_family=profile.family;
            a.alpha_a=profile.alpha_a; a.alpha_c=profile.alpha_c; a.beta_a=profile.beta_a; a.beta_c=profile.beta_c;
            a.tau_intra=profile.tau_intra; a.tau_inter=profile.tau_inter; a.eta_c=profile.eta_c; a.lambda_B=profile.lambda_B;
            a.n_transfer_scores=numel(scores); a.expected_transfer_scores=nPairs;
            a.complete_transfer_set=(numel(scores)==nPairs);
            if a.complete_transfer_set
                a.valid=true; a.mean_validation_score=mean(scores); a.std_validation_score=std(scores,0);
                a.se_validation_score=a.std_validation_score/sqrt(max(numel(scores),1));
                a.min_validation_score=min(scores); a.max_validation_score=max(scores);
                a.mean_within_score=cf_mean_finite(sub.validation_within_score(sub.valid));
                a.mean_cross_score=cf_mean_finite(sub.validation_cross_score(sub.valid));
                a.mean_copy_score=cf_mean_finite(sub.validation_copy_score(sub.valid));
                a.mean_stationarity=cf_mean_finite(sub.stationarity(sub.valid));
            end
            aggRows(c)=a;
        end
    end
    aggregateT=struct2table(aggRows,'AsArray',true);
    aggregateT=sortrows(aggregateT,{'r','mean_validation_score'},{'ascend','ascend'});
    writetable(aggregateT,fullfile(outRoot,'selection_profile_rank_aggregate.csv'));

    rankRows=repmat(cf_empty_rank_row(),numel(ranks),1);
    for ir=1:numel(ranks)
        r=ranks(ir);
        sub=aggregateT(aggregateT.r==r & aggregateT.valid & isfinite(aggregateT.mean_validation_score),:);
        if isempty(sub), error('No complete valid candidate for r=%d.',r); end
        sub=sortrows(sub,'mean_validation_score','ascend'); b=sub(1,:);
        rankRows(ir).r=r; rankRows(ir).best_profile_id=b.profile_id;
        rankRows(ir).best_profile_name=b.profile_name{1};
        rankRows(ir).validation_score=b.mean_validation_score;
        rankRows(ir).validation_std=b.std_validation_score;
        rankRows(ir).validation_se=b.se_validation_score;
        if ir>1
            rankRows(ir).incremental_gain=rankRows(ir-1).validation_score-rankRows(ir).validation_score;
            rankRows(ir).relative_gain=rankRows(ir).incremental_gain/max(abs(rankRows(ir-1).validation_score),eps);
        end
    end
    rr=[rankRows.r]; vv=[rankRows.validation_score];
    [kneeScore,selectedIndex,kneeStatus]=cf_geometric_knee(rr,vv);
    for ir=1:numel(rankRows)
        rankRows(ir).knee_score=kneeScore(ir); rankRows(ir).selected=(ir==selectedIndex);
    end
    rankT=struct2table(rankRows,'AsArray',true);
    writetable(rankT,fullfile(outRoot,'rank_selection_summary.csv'));
    rstar=rankRows(selectedIndex).r;

    atR=aggregateT(aggregateT.r==rstar & aggregateT.valid & isfinite(aggregateT.mean_validation_score),:);
    atR=sortrows(atR,'mean_validation_score','ascend');
    k=min(cfg.top_k_profiles,height(atR)); top1T=atR(1:k,:);
    top1T.auxiliary_rank=(1:k)'; top1T.primary=(1:k)'==1;
    writetable(top1T,fullfile(outRoot,'TOP1_HYPERPARAMETERS_AT_SELECTED_R.csv'));

    cf_plot_selection_results(aggregateT,rankT,top1T,profiles,ranks,outRoot,cfg);
    cf_write_selection_readme(fullfile(outRoot,'SELECTION_README.txt'),cfg,rankT,top1T,kneeStatus);

    selection=struct(); selection.transfer_table=transferT; selection.aggregate_table=aggregateT;
    selection.rank_table=rankT; selection.selected_rank=rstar; selection.selected_rank_index=selectedIndex;
    selection.knee_status=kneeStatus; selection.top1_table=top1T;
end


function row=cf_empty_transfer_row()
    row=struct('source_dataset_id',NaN,'validation_dataset_id',NaN, ...
        'source_offset',NaN,'validation_offset',NaN,'r',NaN,'profile_id',NaN, ...
        'profile_name','','profile_origin','','profile_family','', ...
        'alpha_a',NaN,'alpha_c',NaN,'beta_a',NaN,'beta_c',NaN, ...
        'tau_intra',NaN,'tau_inter',NaN,'eta_c',NaN,'lambda_B',NaN,'seed',NaN, ...
        'completed',false,'valid',false,'status','not_run','error_message','', ...
        'runtime_seconds',NaN,'objective',NaN,'stationarity',NaN,'iterations',NaN, ...
        'best_restart',NaN,'validation_score',NaN,'validation_within_score',NaN, ...
        'validation_cross_score',NaN,'validation_copy_score',NaN, ...
        'graph_total_mass',NaN,'within_mass_fraction',NaN,'cross_mass_fraction',NaN, ...
        'copy_mass_fraction',NaN,'max_B_cosine',NaN,'min_modal_mass_fraction',NaN);
end

function row=cf_empty_aggregate_row()
    row=struct('r',NaN,'profile_id',NaN,'profile_name','','profile_origin','', ...
        'profile_family','','alpha_a',NaN,'alpha_c',NaN,'beta_a',NaN,'beta_c',NaN, ...
        'tau_intra',NaN,'tau_inter',NaN,'eta_c',NaN,'lambda_B',NaN, ...
        'valid',false,'n_transfer_scores',0,'expected_transfer_scores',0, ...
        'complete_transfer_set',false,'mean_validation_score',NaN,'std_validation_score',NaN, ...
        'se_validation_score',NaN,'min_validation_score',NaN,'max_validation_score',NaN, ...
        'mean_within_score',NaN,'mean_cross_score',NaN,'mean_copy_score',NaN, ...
        'mean_stationarity',NaN);
end

function row=cf_empty_rank_row()
    row=struct('r',NaN,'best_profile_id',NaN,'best_profile_name','', ...
        'validation_score',NaN,'validation_std',NaN,'validation_se',NaN, ...
        'incremental_gain',NaN,'relative_gain',NaN,'knee_score',NaN,'selected',false);
end

function [knee,idx,status]=cf_geometric_knee(r,V)
    r=double(r(:)); V=double(V(:)); n=numel(r);
    if n<3, error('At least three candidate ranks are required for elbow selection.'); end
    x=(r-r(1))/max(r(end)-r(1),eps);
    totalDrop=V(1)-V(end);
    if ~isfinite(totalDrop) || totalDrop<=1e-10*max(1,abs(V(1)))
        knee=zeros(n,1); [~,idx]=min(V); status='fallback_minimum_no_overall_decrease'; return;
    end
    y=(V(1)-V)/totalDrop; knee=y-x; knee(1)=0; knee(end)=0;
    interior=2:n-1; [mx,j]=max(knee(interior)); idx=interior(j);
    if ~isfinite(mx)||mx<=0
        [~,idx]=min(V); status='fallback_minimum_no_positive_knee';
    else
        status='geometric_endpoint_chord_knee';
    end
end

function m=cf_mean_finite(x)
    x=x(isfinite(x)); if isempty(x),m=NaN;else,m=mean(x);end
end

function [Sraw,scale,Snorm]=cf_distance_matrix(Y)
    p=size(Y,2); d=sum(Y.^2,2)/p;
    Sraw=max(d+d'-2*(Y*Y'/p),0); Sraw=0.5*(Sraw+Sraw');
    Sraw(1:size(Sraw,1)+1:end)=0;
    vals=Sraw(triu(true(size(Sraw)),1)); vals=vals(isfinite(vals)&vals>0);
    if isempty(vals), error('All positive pairwise distances vanished.'); end
    scale=median(vals); Snorm=Sraw/scale;
end

function vm=cf_graph_validation_metrics(A,S,n,K)
    N=n*K; A=max(0,0.5*(A+A')); A(1:N+1:end)=0;
    maskU=triu(true(N),1); mass=sum(A(maskU));
    if ~isfinite(mass)||mass<=0
        vm=struct('total_score',Inf,'within_score',Inf,'cross_score',Inf,'copy_score',Inf, ...
            'total_mass',mass,'within_mass_fraction',NaN,'cross_mass_fraction',NaN,'copy_mass_fraction',NaN); return;
    end
    Abar=A/mass; within=false(N); cross=false(N); copy=false(N);
    for k=1:K
        Ik=(k-1)*n+(1:n); within(Ik,Ik)=~eye(n);
        for ell=k+1:K
            Il=(ell-1)*n+(1:n); blk=true(n); blk(1:n+1:end)=false; cross(Ik,Il)=blk;
            cp=false(n); cp(1:n+1:end)=true; copy(Ik,Il)=cp;
        end
    end
    within=within&maskU; cross=cross&maskU; copy=copy&maskU;
    vm=struct('total_score',sum(Abar(maskU).*S(maskU)), ...
        'within_score',sum(Abar(within).*S(within)), ...
        'cross_score',sum(Abar(cross).*S(cross)), ...
        'copy_score',sum(Abar(copy).*S(copy)), ...
        'total_mass',mass,'within_mass_fraction',sum(A(within))/mass, ...
        'cross_mass_fraction',sum(A(cross))/mass,'copy_mass_fraction',sum(A(copy))/mass);
end

function gd=cf_graph_diagnostics(fit,n,K)
    r=size(fit.B,2); modeMass=zeros(r,1);
    for m=1:r
        bb=fit.B(:,m)*fit.B(:,m)'; M=bb-diag(diag(bb));
        Am=kron(bb,fit.A_state(:,:,m))+kron(M,diag(fit.Gamma(:,m)));
        modeMass(m)=sum(Am(triu(true(size(Am)),1)));
    end
    if sum(modeMass)>0,mf=modeMass/sum(modeMass);else,mf=zeros(r,1);end
    gd=struct('max_B_cosine',cf_max_column_cosine(fit.B), ...
        'min_modal_mass_fraction',min(mf),'max_modal_mass_fraction',max(mf));
end

function c=cf_max_column_cosine(B)
    r=size(B,2); c=NaN; if r<2,return;end; c=0;
    for i=1:r,for j=i+1:r,c=max(c,(B(:,i)'*B(:,j))/max(norm(B(:,i))*norm(B(:,j)),eps));end,end
end


%% ========================================================================
%% SELECTION FIGURES / REPORTS
%% ========================================================================
function cf_plot_selection_results(aggregateT,rankT,top1T,profiles,ranks,outRoot,cfg)
    % Heatmap: profiles x ranks.
    M=nan(numel(profiles),numel(ranks));
    for ip=1:numel(profiles)
        for ir=1:numel(ranks)
            z=aggregateT.profile_id==profiles(ip).profile_id & aggregateT.r==ranks(ir);
            if any(z), M(ip,ir)=aggregateT.mean_validation_score(find(z,1)); end
        end
    end
    f=cf_new_figure([1000 1500]); imagesc(M); colorbar; axis tight;
    xlabel('candidate rank r'); ylabel('hyperparameter profile ID');
    set(gca,'XTick',1:numel(ranks),'XTickLabel',string(ranks));
    title('Mean validation-configuration smoothness: lower is better');
    cf_export_figure(f,fullfile(outRoot,'01_profile_rank_validation_heatmap.png'),cfg);

    f=cf_new_figure([1050 700]);
    errorbar(rankT.r,rankT.validation_score,rankT.validation_se,'-o','LineWidth',1.5,'MarkerSize',7);
    grid on; box on; xlabel('modal rank r'); ylabel('V^*(r): best mean validation-configuration smoothness');
    title('Rank selection curve');
    xline(rankT.r(rankT.selected),'--','selected r');
    cf_export_figure(f,fullfile(outRoot,'02_rank_elbow_curve.png'),cfg);

    f=cf_new_figure([1050 700]);
    plot(rankT.r,rankT.incremental_gain,'-o','LineWidth',1.5,'MarkerSize',7); grid on; box on;
    xlabel('modal rank r'); ylabel('V^*(r-1)-V^*(r)'); title('Incremental validation gain by rank');
    cf_export_figure(f,fullfile(outRoot,'03_incremental_gain.png'),cfg);

    f=cf_new_figure([1050 700]);
    plot(rankT.r,rankT.knee_score,'-o','LineWidth',1.5,'MarkerSize',7); grid on; box on;
    xlabel('modal rank r'); ylabel('geometric knee score'); title('Endpoint-chord elbow score');
    xline(rankT.r(rankT.selected),'--','selected r');
    cf_export_figure(f,fullfile(outRoot,'04_knee_score.png'),cfg);

    f=cf_new_figure([1200 700]); bar(1:height(top1T),top1T.mean_validation_score); grid on; box on;
    xticks(1:height(top1T)); xticklabels(top1T.profile_name); xtickangle(25);
    ylabel('mean validation-configuration smoothness'); title(sprintf('Selected top-1 at r=%d',rankT.r(rankT.selected)));
    cf_export_figure(f,fullfile(outRoot,'05_top1_hyperparameter_scores.png'),cfg);
end

function cf_write_selection_readme(path,cfg,rankT,top1T,kneeStatus)
    fid=fopen(path,'w'); if fid<0,return;end; cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
    fprintf(fid,'MOSAIC FINANCE OFFSET CROSS-CONFIGURATION SELECTION\n');
    fprintf(fid,'===================================================\n\n');
    fprintf(fid,'Selection transfers: dataset 2(offset1)->dataset 3(offset2), dataset 4(offset3)->dataset 5(offset4).\n');
    fprintf(fid,'Dataset 1(offset0) is excluded from hyperparameter and rank selection.\n');
    fprintf(fid,'Each source fit uses all of its signal coordinates; no feature-coordinate holdout is used.\n');
    fprintf(fid,'Hyperparameter bank: 64 reference profiles + %d hp59-local profiles.\n',cfg.local_neighbor_count);
    fprintf(fid,'Ranks: %s\n',mat2str(cfg.rank_candidates));
    fprintf(fid,'Selection restarts per exact source/r/profile: %d\n',cfg.selection_n_restarts);
    fprintf(fid,'Rank rule: V*(r)=min_profile mean validation-configuration smoothness; geometric endpoint-chord elbow.\n');
    fprintf(fid,'Elbow status: %s\n',kneeStatus);
    fprintf(fid,'Selected r: %d\n\n',rankT.r(rankT.selected));
    fprintf(fid,'Selected top-1 configuration at the chosen rank (target not consulted):\n');
    for q=1:height(top1T)
        fprintf(fid,'  %d) profile %d %s | Vbar=%.12g\n',q,top1T.profile_id(q), ...
            top1T.profile_name{q},top1T.mean_validation_score(q));
    end
end


%% ========================================================================
%% FINAL TARGET FITS + FULL OUTPUTS
%% ========================================================================
function final=cf_run_final_target(targetData,targetPrepared,profiles,selection,cfg,outRoot)
    rstar=selection.selected_rank; top1=selection.top1_table;
    writetable(top1,fullfile(outRoot,'SELECTED_TOP1_FROM_OFFSET_VALIDATION.csv'));
    Y=targetPrepared.Y; observed=targetPrepared.observed; observed.r=rstar;
    inputInfo=targetPrepared.inputInfo;
    pData=pcm_prepare_data(targetData);

    nTop=height(top1); rows=repmat(cf_empty_final_row(),nTop,1);
    allRelations=cell(nTop,1); allModes=cell(nTop,1); allCouplings=cell(nTop,1);
    allWithinEdges=cell(nTop,1); allModeEdges=cell(nTop,1); allRankings=cell(nTop,1);
    contexts=cell(nTop,1);

    pairedSeed=cfg.seed+900000+10000*rstar;
    for q=1:nTop
        pid=top1.profile_id(q); ix=find([profiles.profile_id]==pid,1); profile=profiles(ix);
        candidateId=sprintf('top%02d_r%d_profile%03d_%s',q,rstar,pid,cf_safe(profile.slug));
        candDir=fullfile(outRoot,candidateId); tablesDir=fullfile(candDir,'tables');
        figDir=fullfile(candDir,'figures'); viewDir=fullfile(candDir,'view_graph_clustering');
        postDir=fullfile(candDir,'cross_within_mode_postprocess'); matrixDir=fullfile(postDir,'matrices');
        postTableDir=fullfile(postDir,'tables'); postFigDir=fullfile(postDir,'figures');
        pairDir=fullfile(postFigDir,'individual_pair_diagnostics');
        cf_ensure_folder(candDir); cf_ensure_folder(tablesDir); cf_ensure_folder(figDir);
        cf_ensure_folder(viewDir); cf_ensure_folder(postDir); cf_ensure_folder(matrixDir);
        cf_ensure_folder(postTableDir); cf_ensure_folder(postFigDir);
        if cfg.save_individual_pair_figures, cf_ensure_folder(pairDir); end

        fitFile=fullfile(candDir,'selected_representative_fit.mat'); rec=cf_empty_final_row();
        reuse=false; compactFit=[]; analysis=[]; finalRelations=[];
        if cfg.resume && isfile(fitFile)
            Z=load(fitFile);
            if isfield(Z,'record') && Z.record.valid && isfield(Z,'compactFit')
                rec=Z.record; compactFit=Z.compactFit; analysis=Z.analysis;
                if isfield(Z,'finalRelations'),finalRelations=Z.finalRelations;end
                reuse=true;
            end
        end
        if ~reuse
            fprintf('\n[FINAL TARGET] %s\n',candidateId);
            opts=cf_profile_to_options(profile,cfg,pairedSeed,'final'); timer=tic;
            try
                fit=solve_mosaic_crossview_embedded(observed,opts);
                rec.runtime_seconds=toc(timer); rec.valid=true; rec.status=fit.status;
                rec.r=rstar; rec.auxiliary_rank=q; rec.profile_id=pid; rec.profile_name=profile.name;
                rec.auxiliary_validation_score=top1.mean_validation_score(q);
                rec.objective=fit.objective; rec.stationarity=fit.stationarity_combined;
                rec.iterations=fit.iterations; rec.best_restart=fit.best_restart;
                rec.seed=pairedSeed;
                analysis=finance_analyze_fit(fit,targetData,inputInfo,cfg);
                compactFit=finance_compact_fit(fit);
                finalRelations=pcm_compute_top_cross_relations(compactFit.A_supra,pData,cfg.final_top_relations);
                record=rec; %#ok<NASGU>
                save(fitFile,'compactFit','analysis','profile','record','finalRelations','inputInfo','-v7.3');
            catch ME
                rec.runtime_seconds=toc(timer); rec.status=['error:' ME.identifier]; rec.error_message=ME.message;
                record=rec; save(fitFile,'record','profile');
                rows(q)=rec; warning('Final fit failed for %s: %s',candidateId,ME.message); continue;
            end
        end

        if isempty(finalRelations)
            finalRelations=pcm_compute_top_cross_relations(compactFit.A_supra,pData,cfg.final_top_relations);
        end
        if isempty(analysis), analysis=finance_analyze_fit(compactFit,targetData,inputInfo,cfg); end

        % Earlier stability-run style outputs.
        finance_save_candidate_tables(targetData,compactFit,analysis,tablesDir,candidateId,cfg);
        finance_save_candidate_figures(targetData,Y,compactFit,analysis,figDir,candidateId,cfg);
        viewResult=stability_view_graph_postprocess(compactFit,targetData,inputInfo,viewDir,candidateId,cfg); %#ok<NASGU>
        writetable(finalRelations,fullfile(tablesDir,'final_top_cross_relations_unconstrained.csv'));

        % Cross-vs-within/mode v3 outputs.
        within=pcm_reconstruct_within(compactFit,pData);
        saved=struct('finalRelations',finalRelations);
        topRel=pcm_get_top_cross_relations(saved,compactFit,pData,cfg.top_cross_relations);
        [relationAnalysis,modeLong,couplingLong,pairMatrices]= ...
            pcm_analyze_cross_relations(topRel,compactFit,within,pData,candidateId,rstar);
        rankingTable=pcm_extract_ranking_table(relationAnalysis);
        withinEdges=pcm_top_within_edges(within,pData,cfg.top_within_edges_per_horizon,candidateId,rstar);
        modeEdges=pcm_top_mode_edges(compactFit.A_state,pData,cfg.top_mode_edges,candidateId,rstar);

        writetable(relationAnalysis,fullfile(postTableDir,'top_cross_relations_cross_vs_within.csv'));
        writetable(rankingTable,fullfile(postTableDir,'top_cross_relations_cross_vs_within_ranking.csv'));
        writetable(modeLong,fullfile(postTableDir,'top_cross_relations_mode_decomposition_long.csv'));
        writetable(couplingLong,fullfile(postTableDir,'top_cross_stock_pair_horizon_coupling_long.csv'));
        writetable(withinEdges,fullfile(postTableDir,'top_within_view_edges_all_horizons.csv'));
        writetable(modeEdges,fullfile(postTableDir,'top_mode_A_edges_all_modes.csv'));
        pcm_save_complete_matrices(compactFit,within,pData,matrixDir);
        pcm_plot_cross_vs_within(relationAnalysis,postFigDir,candidateId,cfg);
        pcm_plot_cross_vs_within_ranks(rankingTable,postFigDir,candidateId,cfg);
        pcm_plot_pair_coupling_grid(pairMatrices,relationAnalysis,pData,postFigDir,candidateId,cfg);
        pcm_plot_within_graphs(within,pData,postFigDir,candidateId,cfg);
        pcm_plot_mode_graphs(compactFit,pData,postFigDir,candidateId,cfg);
        pcm_plot_mode_decomposition(modeLong,relationAnalysis,postFigDir,candidateId,rstar,cfg);
        if cfg.save_individual_pair_figures
            pcm_save_individual_pair_figures(pairMatrices,relationAnalysis,modeLong,pData,pairDir,candidateId,cfg);
        end

        context=struct('representative_id',candidateId,'rank',rstar,'profile',profile, ...
            'relation_analysis',relationAnalysis,'ranking_analysis',rankingTable, ...
            'mode_decomposition',modeLong,'horizon_coupling',couplingLong, ...
            'within_edges',withinEdges,'mode_edges',modeEdges,'within_adjacencies',within, ...
            'mode_adjacencies',compactFit.A_state,'B',compactFit.B,'Gamma',compactFit.Gamma, ...
            'pair_coupling_matrices',{pairMatrices});
        save(fullfile(postDir,'cross_within_mode_postprocess.mat'),'context','-v7.3');

        allRelations{q}=relationAnalysis; allModes{q}=modeLong; allCouplings{q}=couplingLong;
        allWithinEdges{q}=withinEdges; allModeEdges{q}=modeEdges; allRankings{q}=rankingTable;
        contexts{q}=context; rows(q)=rec;
    end

    summaryT=struct2table(rows,'AsArray',true);
    writetable(summaryT,fullfile(outRoot,'FINAL_TOP1_TARGET_FITS.csv'));

    summaryDir=fullfile(outRoot,'ACROSS_TOP1_POSTPROCESS_SUMMARY'); cf_ensure_folder(summaryDir);
    R=pcm_vertcat_nonempty(allRelations); M=pcm_vertcat_nonempty(allModes);
    C=pcm_vertcat_nonempty(allCouplings); W=pcm_vertcat_nonempty(allWithinEdges);
    A=pcm_vertcat_nonempty(allModeEdges); Q=pcm_vertcat_nonempty(allRankings);
    if ~isempty(R)
        relationConsensus=pcm_relation_consensus(R,nTop);
        stockPairConsensus=pcm_stock_pair_consensus(R,nTop);
        writetable(R,fullfile(summaryDir,'all_top1_cross_vs_within_long.csv'));
        writetable(M,fullfile(summaryDir,'all_top1_mode_decomposition_long.csv'));
        writetable(C,fullfile(summaryDir,'all_top1_pair_coupling_long.csv'));
        writetable(W,fullfile(summaryDir,'all_top1_top_within_edges.csv'));
        writetable(A,fullfile(summaryDir,'all_top1_top_mode_edges.csv'));
        writetable(Q,fullfile(summaryDir,'all_top1_cross_vs_within_ranking_long.csv'));
        writetable(relationConsensus,fullfile(summaryDir,'cross_relation_recurrence_and_within_comparison.csv'));
        writetable(stockPairConsensus,fullfile(summaryDir,'stock_pair_multiscale_recurrence_summary.csv'));
        pcm_plot_across_representative_consensus(relationConsensus,summaryDir,cfg);
        pcm_plot_stock_pair_multiscale_summary(stockPairConsensus,summaryDir,cfg);
    else
        relationConsensus=table(); stockPairConsensus=table();
    end

    if nTop>=1 && rows(1).valid
        Z=load(fullfile(outRoot,sprintf('top%02d_r%d_profile%03d_%s',1,rstar,top1.profile_id(1), ...
            cf_safe(profiles(find([profiles.profile_id]==top1.profile_id(1),1)).slug)), ...
            'selected_representative_fit.mat'));
        primaryFit=Z.compactFit; primaryAnalysis=Z.analysis; primaryProfile=Z.profile; primaryRecord=Z.record; %#ok<NASGU>
        save(fullfile(outRoot,'PRIMARY_FINAL_MODEL_TOP1_AUXILIARY_VALIDATION.mat'), ...
            'primaryFit','primaryAnalysis','primaryProfile','primaryRecord','rstar','-v7.3');
    end

    final=struct('summary_table',summaryT,'r',rstar,'selected_top1_table',top1, ...
        'contexts',{contexts},'relation_consensus',relationConsensus, ...
        'stock_pair_consensus',stockPairConsensus,'target_used_for_reselection',false);
end

function row=cf_empty_final_row()
    row=struct('r',NaN,'auxiliary_rank',NaN,'profile_id',NaN,'profile_name','', ...
        'auxiliary_validation_score',NaN,'seed',NaN,'valid',false,'status','not_run', ...
        'error_message','','runtime_seconds',NaN,'objective',NaN,'stationarity',NaN, ...
        'iterations',NaN,'best_restart',NaN);
end


%% ========================================================================
%% MASTER README / SMALL UTILITIES
%% ========================================================================
function cf_write_master_readme(path,cfg,selection,final)
    fid=fopen(path,'w'); if fid<0,return;end; cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
    fprintf(fid,'MOSAIC FINANCE FIVE-OFFSET VALIDATION RANK-SELECTION PIPELINE\n');
    fprintf(fid,'========================================================\n\n');
    fprintf(fid,'Five endpoint residue classes modulo endpointStep=%d are generated.\n',cfg.endpoint_step);
    fprintf(fid,'Dataset 1 offset0 = untouched final target for model selection.\n');
    fprintf(fid,'Selection transfers: dataset2(offset1)->dataset3(offset2), dataset4(offset3)->dataset5(offset4).\n');
    fprintf(fid,'These are local deterministic offset-validation configurations; no statistical independence claim is made.\n\n');
    fprintf(fid,'Horizons: %s\n',mat2str(cfg.horizons));
    fprintf(fid,'Ranks searched: %s\n',mat2str(cfg.rank_candidates));
    fprintf(fid,'Hyperparameters: 64 reference finance profiles plus %d LHS neighbors around profile59/lhs51.\n',cfg.local_neighbor_count);
    fprintf(fid,'Selection restarts: %d; final restarts: %d.\n\n',cfg.selection_n_restarts,cfg.final_n_restarts);
    if ~isempty(selection)
        fprintf(fid,'Selected r*: %d (%s).\n',selection.selected_rank,selection.knee_status);
        fprintf(fid,'Selected top-1 configuration:\n');
        for q=1:height(selection.top1_table)
            fprintf(fid,'  %d) profile %d %s | validation Vbar=%.12g\n',q, ...
                selection.top1_table.profile_id(q),selection.top1_table.profile_name{q}, ...
                selection.top1_table.mean_validation_score(q));
        end
    end
    if ~isempty(final)
        fprintf(fid,'\nThe offset-0 target is fit from scratch with the selected top-1 configuration.\n');
        fprintf(fid,'Primary remains auxiliary top-1; target data do not re-rank profiles or rank.\n');
    end
end

function s=cf_safe(s)
    s=regexprep(char(s),'[^A-Za-z0-9_\-]','_');
end

function f=cf_new_figure(sz)
    f=figure('Visible','off','Color','w','Position',[80 80 sz]);
end

function cf_export_figure(f,path,cfg)
    [folder,~,~]=fileparts(path); cf_ensure_folder(folder);
    try, exportgraphics(f,path,'Resolution',cfg.figure_resolution); catch, print(f,path,'-dpng',sprintf('-r%d',cfg.figure_resolution)); end
    close(f);
end

function cf_ensure_folder(path)
    if ~exist(path,'dir'), mkdir(path); end
end

function out=cf_merge(base,override)
    out=base; if isempty(override),return;end
    f=fieldnames(override); for q=1:numel(f),out.(f{q})=override.(f{q});end
end


%% ========================================================================
%% EXACT ORIGINAL FINANCE PROFILE-GENERATION HELPERS
%% ========================================================================
function profiles = stability_generate_profiles(cfg)
    template = struct( ...
        'profile_id',NaN,'name','','slug','','family','', ...
        'alpha_a',NaN,'alpha_c',NaN, ...
        'beta_a',NaN,'beta_c',NaN, ...
        'tau_intra',NaN,'tau_inter',NaN, ...
        'eta_c',NaN,'lambda_B',NaN, ...
        'epsilon_intra',1e-4,'epsilon_inter',1e-4);

    anchors = repmat(template,8,1);
    anchors(1) = stability_make_profile(template,1, ...
        'weak_regularization','weak','anchor_weak', ...
        .001,.001,.03,.03,.04,.04,1.00,.00);
    anchors(2) = stability_make_profile(template,2, ...
        'balanced_low','bal_low','anchor_balanced', ...
        .003,.003,.06,.06,.06,.06,1.00,.02);
    anchors(3) = stability_make_profile(template,3, ...
        'balanced','balanced','anchor_balanced', ...
        .010,.010,.10,.10,.10,.10,1.00,.03);
    anchors(4) = stability_make_profile(template,4, ...
        'strong_connectivity','connect','anchor_connectivity', ...
        .004,.004,.15,.15,.22,.22,1.00,.03);
    anchors(5) = stability_make_profile(template,5, ...
        'copy_persistence_emphasis','copy','anchor_copy', ...
        .006,.001,.08,.04,.08,.25,.50,.03);
    anchors(6) = stability_make_profile(template,6, ...
        'cross_stock_emphasis','cross','anchor_cross', ...
        .001,.008,.04,.12,.08,.25,1.00,.03);
    anchors(7) = stability_make_profile(template,7, ...
        'within_horizon_emphasis','within','anchor_within', ...
        .003,.012,.08,.15,.20,.05,1.00,.03);
    anchors(8) = stability_make_profile(template,8, ...
        'strong_mode_diversity','diverse','anchor_diverse', ...
        .006,.006,.08,.08,.10,.10,1.00,.12);

    nTotal = cfg.n_profiles_per_rank;
    if nTotal < numel(anchors)
        profiles = anchors(1:nTotal);
        return;
    end

    nLhs = nTotal-numel(anchors);
    profiles = repmat(template,nTotal,1);
    profiles(1:numel(anchors)) = anchors;

    if nLhs == 0
        return;
    end

    U = stability_latin_hypercube(nLhs,8,cfg.seed+271828);
    for q = 1:nLhs
        alphaBase = stability_log_map(U(q,1),cfg.alpha_base_range);
        alphaRatio = stability_log_map(U(q,2),cfg.alpha_ratio_range);
        betaBase = stability_log_map(U(q,3),cfg.beta_base_range);
        betaRatio = stability_log_map(U(q,4),cfg.beta_ratio_range);
        tauBase = stability_log_map(U(q,5),cfg.tau_base_range);
        tauRatio = stability_log_map( ...
            U(q,6),cfg.tau_inter_intra_ratio_range);
        etaC = stability_log_map(U(q,7),cfg.eta_c_range);
        lambdaB = stability_linear_map(U(q,8),cfg.lambda_B_range);

        aa = alphaBase*sqrt(alphaRatio);
        ac = alphaBase/sqrt(alphaRatio);
        ba = betaBase*sqrt(betaRatio);
        bc = betaBase/sqrt(betaRatio);
        ti = tauBase/sqrt(tauRatio);
        te = tauBase*sqrt(tauRatio);

        id = numel(anchors)+q;
        family = stability_classify_profile( ...
            aa,ac,ba,bc,ti,te,etaC,lambdaB);
        profiles(id) = stability_make_profile(template,id, ...
            sprintf('lhs_%02d',q),sprintf('lhs%02d',q),family, ...
            aa,ac,ba,bc,ti,te,etaC,lambdaB);
    end
end


function p = stability_make_profile(template,id,name,slug,family, ...
        aa,ac,ba,bc,ti,te,eta,lb)
    p = template;
    p.profile_id = id;
    p.name = name;
    p.slug = slug;
    p.family = family;
    p.alpha_a = aa;
    p.alpha_c = ac;
    p.beta_a = ba;
    p.beta_c = bc;
    p.tau_intra = ti;
    p.tau_inter = te;
    p.eta_c = eta;
    p.lambda_B = lb;
end


function family = stability_classify_profile(aa,ac,ba,bc,ti,te,eta,lb)
    if lb >= 0.085
        family = 'mode_diversity';
    elseif te/ti >= 1.75 && eta <= 0.80
        family = 'copy_persistence';
    elseif ac/aa >= 1.75 || bc/ba >= 1.60
        family = 'cross_stock';
    elseif ti/te >= 1.75
        family = 'within_horizon';
    elseif max([ba,bc,ti,te]) >= 0.18
        family = 'strong_connectivity';
    else
        family = 'balanced';
    end
end


function U = stability_latin_hypercube(n,d,seed)
% Toolbox-free randomized Latin hypercube on (0,1)^d.
    rng(seed,'twister');
    U = zeros(n,d);
    for j = 1:d
        perm = randperm(n);
        jitter = rand(n,1);
        U(:,j) = (perm(:)-jitter)/n;
    end
    U = min(max(U,eps),1-eps);
end


function value = stability_log_map(u,bounds)
    value = exp(log(bounds(1)) + u*(log(bounds(2))-log(bounds(1))));
end


function value = stability_linear_map(u,bounds)
    value = bounds(1) + u*(bounds(2)-bounds(1));
end


function T = stability_profiles_to_table(profiles)
    T = struct2table(profiles);
    T.alpha_a_over_alpha_c = T.alpha_a./T.alpha_c;
    T.beta_a_over_beta_c = T.beta_a./T.beta_c;
    T.tau_inter_over_tau_intra = T.tau_inter./T.tau_intra;
end



%% ========================================================================
%% ORIGINAL DATA-BUILDER INPUT/PARSING HELPERS
%% ========================================================================
function dates = parse_dates_strict(rawDate)
% Parse a date column. Fails rather than inventing artificial dates because
% calendar alignment is essential for this dataset.

    if isdatetime(rawDate)
        dates = rawDate(:);
        return;
    end

    if isnumeric(rawDate)
        v = double(rawDate(:));
        dates = NaT(size(v));

        % Unix seconds
        if median(v,'omitnan') > 1e8
            try
                dates = datetime(v, 'ConvertFrom', 'posixtime', ...
                    'TimeZone', 'UTC');
                dates.TimeZone = '';
                return;
            catch
            end
        end

        % MATLAB datenum / Excel-like serial date fallback
        try
            tmp = datetime(v, 'ConvertFrom', 'datenum');
            if nnz(~isnat(tmp)) >= 0.9*numel(v)
                dates = tmp;
                return;
            end
        catch
        end
    end

    s = string(rawDate(:));
    s = strtrim(s);
    dates = NaT(size(s));

    formats = { ...
        'yyyy-MM-dd', ...
        'yyyy/MM/dd', ...
        'MM/dd/yyyy', ...
        'dd/MM/yyyy', ...
        'dd-MMM-yyyy', ...
        'yyyy-MM-dd HH:mm:ss', ...
        'MM/dd/yyyy HH:mm:ss'};

    bestDates = dates;
    bestCount = 0;
    for i = 1:numel(formats)
        try
            tmp = datetime(s, 'InputFormat', formats{i});
            count = nnz(~isnat(tmp));
            if count > bestCount
                bestCount = count;
                bestDates = tmp;
            end
        catch
        end
    end

    if bestCount < 0.9*numel(s)
        try
            tmp = datetime(s);
            count = nnz(~isnat(tmp));
            if count > bestCount
                bestCount = count;
                bestDates = tmp;
            end
        catch
        end
    end

    if bestCount < 0.9*numel(s)
        error(['Could not reliably parse the date column. Parsed %d of %d ' ...
               'entries. Check the CSV date format.'], bestCount, numel(s));
    end

    dates = bestDates(:);
end

function [X, keep] = table_to_numeric_matrix_with_mask(T)
% Convert numeric or numeric-looking table columns to a matrix while
% retaining a mask that keeps the variable names aligned.

    nRows = height(T);
    nCols = width(T);
    Xall = nan(nRows,nCols);
    keep = false(1,nCols);

    for j = 1:nCols
        v = T{:,j};

        if isnumeric(v) || islogical(v)
            xj = double(v);
        elseif iscell(v) || isstring(v) || ischar(v) || iscategorical(v)
            xj = str2double(string(v));
        else
            try
                xj = double(v);
            catch
                xj = nan(nRows,1);
            end
        end

        xj = xj(:);
        if numel(xj) == nRows && any(isfinite(xj))
            Xall(:,j) = xj;
            keep(j) = true;
        end
    end

    X = Xall(:,keep);
end

function key = normalize_ticker(tickers)
% Normalize common ticker punctuation differences for matching only.
    key = upper(strtrim(string(tickers)));
    key = replace(key, '.', '-');
    key = replace(key, '/', '-');
    key = erase(key, ' ');
end

function [found, idx] = match_tickers(requested, available)
    req = normalize_ticker(requested);
    avail = normalize_ticker(available);
    [found, idx] = ismember(req, avail);
end

function [sectorLabels, sectorKnown] = load_sector_labels(sectorCsv, tickers)
    S = readtable(sectorCsv, 'VariableNamingRule', 'preserve');
    names = string(S.Properties.VariableNames);
    lowerNames = lower(names);

    symbolCol = find(contains(lowerNames,'symbol') | ...
                     contains(lowerNames,'ticker'), 1, 'first');
    if isempty(symbolCol)
        symbolCol = 1;
    end

    sectorCol = find(contains(lowerNames,'gics') & ...
                     contains(lowerNames,'sector'), 1, 'first');
    if isempty(sectorCol)
        sectorCol = find(contains(lowerNames,'sector'), 1, 'first');
    end
    if isempty(sectorCol)
        sectorCol = min(2,width(S));
    end

    symbols = string(S{:,symbolCol});
    sectors = string(S{:,sectorCol});
    symbolsKey = normalize_ticker(symbols);
    tickerKey = normalize_ticker(tickers);

    sectorLabels = strings(numel(tickers),1);
    sectorKnown = false(numel(tickers),1);

    for i = 1:numel(tickers)
        hit = find(symbolsKey == tickerKey(i), 1, 'first');
        if ~isempty(hit) && strlength(strtrim(sectors(hit))) > 0
            sectorLabels(i) = strtrim(sectors(hit));
            sectorKnown(i) = true;
        else
            sectorLabels(i) = "Unknown";
        end
    end
end

function [order, sectorNames, boundaries] = make_sector_order(sectorLabels)
    labels = string(sectorLabels(:));
    labels(strlength(labels)==0) = "Unknown";

    [sectorNames,~,groupId] = unique(labels,'stable');
    [groupSorted,order] = sort(groupId,'ascend');

    boundaries = [];
    for s = 1:numel(sectorNames)-1
        lastIdx = find(groupSorted == s, 1, 'last');
        if ~isempty(lastIdx)
            boundaries(end+1) = lastIdx + 0.5; %#ok<AGROW>
        end
    end
end


%% ========================================================================
%% ORIGINAL STABILITY VIEW-GRAPH / CLUSTERING HELPERS
%% ========================================================================
function viewResult = stability_view_graph_postprocess( ...
        fit,data,inputInfo,outDir,candidateId,cfg)

    B = double(fit.B);
    Theta = double(fit.Theta);
    Gamma = double(fit.Gamma);
    horizons = double(data.horizons(:));
    K = numel(horizons);

    sA = sum(Theta,1).';
    sC = sum(Gamma,1).';
    totalModeMass = 2*sA+sC;

    Gcross = B*diag(2*sA)*B.';
    Gcopy = B*diag(sC)*B.';
    Gview = Gcross+Gcopy;

    Gcross = stability_remove_diagonal( ...
        stability_symmetrize_nonnegative(Gcross));
    Gcopy = stability_remove_diagonal( ...
        stability_symmetrize_nonnegative(Gcopy));
    Gview = stability_remove_diagonal( ...
        stability_symmetrize_nonnegative(Gview));

    modalGraphs = zeros(K,K,size(B,2));
    for m = 1:size(B,2)
        modalGraphs(:,:,m) = stability_remove_diagonal( ...
            totalModeMass(m)*(B(:,m)*B(:,m).'));
    end

    savedRelativeError = NaN;
    if isfield(fit,'view_graph') && isequal(size(fit.view_graph),[K K])
        savedRelativeError = norm(Gview-double(fit.view_graph),'fro') / ...
            max(norm(Gview,'fro'),eps);
    end

    n = data.num_stocks;
    blockSum = zeros(K,K);
    for k = 1:K
        Ik = (k-1)*n+(1:n);
        for ell = 1:K
            if k==ell
                continue;
            end
            Il = (ell-1)*n+(1:n);
            blockSum(k,ell) = sum(fit.A_supra(Ik,Il),'all');
        end
    end
    blockSum = stability_remove_diagonal( ...
        stability_symmetrize_nonnegative(blockSum));
    blockRelativeError = norm(Gview-blockSum,'fro') / ...
        max(norm(blockSum,'fro'),eps);

    clusterInfo = stability_cluster_view_graph(Gview,horizons,fit.B,cfg);

    [edgeI,edgeJ] = find(triu(true(K),1));
    edgeIndex = sub2ind([K K],edgeI,edgeJ);
    totalWeight = Gview(edgeIndex);
    crossWeight = Gcross(edgeIndex);
    copyWeight = Gcopy(edgeIndex);
    logSeparation = abs(log(horizons(edgeI)./horizons(edgeJ)));
    rawSimilarity = nan(size(totalWeight));
    if isfield(inputInfo,'raw_view_similarity') && ...
            isequal(size(inputInfo.raw_view_similarity),[K K])
        rawSimilarity = inputInfo.raw_view_similarity(edgeIndex);
    end
    edgeTable = table(edgeI,edgeJ,horizons(edgeI),horizons(edgeJ), ...
        logSeparation,totalWeight,crossWeight,copyWeight,rawSimilarity, ...
        'VariableNames',{'view_index_1','view_index_2', ...
        'horizon_1_days','horizon_2_days','absolute_log_horizon_separation', ...
        'total_view_weight','crossnode_contribution', ...
        'copy_contribution','raw_same_stock_similarity'});
    edgeTable = sortrows(edgeTable,'total_view_weight','descend');

    stability_write_view_matrix(fullfile(outDir,'G_view_total.csv'), ...
        Gview,horizons);
    stability_write_view_matrix(fullfile(outDir, ...
        'G_view_crossnode_contribution.csv'),Gcross,horizons);
    stability_write_view_matrix(fullfile(outDir, ...
        'G_view_copy_contribution.csv'),Gcopy,horizons);
    stability_write_view_matrix(fullfile(outDir, ...
        'G_view_blocksum_verification.csv'),blockSum,horizons);
    writetable(edgeTable,fullfile(outDir,'view_graph_edge_list.csv'));

    massTable = table((1:numel(sA)).',sA,sC,totalModeMass, ...
        'VariableNames',{'mode','crossnode_state_mass_s_a', ...
        'copy_state_mass_s_c','total_view_mass_2s_a_plus_s_c'});
    writetable(massTable,fullfile(outDir,'modal_view_masses.csv'));

    metricsTable = stability_cluster_metrics_table(clusterInfo);
    writetable(metricsTable,fullfile(outDir, ...
        'spectral_cluster_metrics.csv'));

    assignments = table(horizons,'VariableNames',{'horizon_days'});
    for q = 1:numel(clusterInfo.sweep)
        c = clusterInfo.sweep(q).cluster_count;
        assignments.(sprintf('cluster_c%d',c)) = ...
            clusterInfo.sweep(q).labels;
    end
    writetable(assignments,fullfile(outDir, ...
        'view_cluster_assignments.csv'));

    stability_plot_view_graph_components( ...
        Gview,Gcross,Gcopy,B,sA,sC,horizons, ...
        savedRelativeError,blockRelativeError,candidateId, ...
        fullfile(outDir,'01_view_graph_components.png'),cfg);
    stability_plot_view_spectral_diagnostics( ...
        clusterInfo,candidateId, ...
        fullfile(outDir,'02_spectral_cluster_diagnostics.png'),cfg);
    stability_plot_all_view_clusterings( ...
        Gview,horizons,clusterInfo,candidateId, ...
        fullfile(outDir,'03_all_view_clusterings.png'),cfg);

    rankMatchedItem = stability_get_cluster_item( ...
        clusterInfo,clusterInfo.rank_matched_count);
    stability_plot_main_view_clustering( ...
        Gview,horizons,rankMatchedItem,candidateId,'rank-matched', ...
        fullfile(outDir,'04_rank_matched_clustering.png'),cfg);

    eigengapItem = stability_get_cluster_item( ...
        clusterInfo,clusterInfo.eigengap_suggested_count);
    stability_plot_main_view_clustering( ...
        Gview,horizons,eigengapItem,candidateId,'eigengap diagnostic', ...
        fullfile(outDir,'05_eigengap_clustering.png'),cfg);

    stability_plot_modal_view_graphs( ...
        modalGraphs,B,sA,sC,horizons,candidateId, ...
        fullfile(outDir,'06_modal_view_graph_decomposition.png'),cfg);

    viewResult = struct();
    viewResult.G_view = Gview;
    viewResult.G_cross = Gcross;
    viewResult.G_copy = Gcopy;
    viewResult.G_modal = modalGraphs;
    viewResult.state_cross_mass = sA;
    viewResult.state_copy_mass = sC;
    viewResult.state_total_mass = totalModeMass;
    viewResult.saved_view_graph_relative_error = savedRelativeError;
    viewResult.block_sum_relative_error = blockRelativeError;
    viewResult.edge_table = edgeTable;
    viewResult.clustering = clusterInfo;

    save(fullfile(outDir,'view_graph_reconstruction_and_clustering.mat'), ...
        'viewResult','-v7.3');
end


function clusterInfo = stability_cluster_view_graph(W,horizons,B,cfg)
    K = size(W,1);
    W = stability_remove_diagonal(stability_symmetrize_nonnegative(W));
    maxWeight = max(W(:));
    if ~isfinite(maxWeight) || maxWeight<=0
        error('The induced view graph is identically zero.');
    end
    Wnorm = W/maxWeight;

    degree = sum(Wnorm,2);
    invSqrtDegree = zeros(K,1);
    positive = degree>eps;
    invSqrtDegree(positive) = 1./sqrt(degree(positive));
    Lsym = eye(K) - ...
        (invSqrtDegree.*Wnorm).*invSqrtDegree.';
    Lsym = 0.5*(Lsym+Lsym.');

    [V,D] = eig(Lsym);
    eigenvalues = real(diag(D));
    [eigenvalues,order] = sort(eigenvalues,'ascend');
    V = real(V(:,order));

    counts = cfg.view_cluster_counts;
    counts = counts(counts>=2 & counts<=K-1);
    counts = unique(counts(:).');
    if isempty(counts)
        error('No valid view-cluster counts remain.');
    end

    sweep = repmat(stability_empty_cluster_result(),numel(counts),1);
    for q = 1:numel(counts)
        c = counts(q);
        embedding = V(:,1:c);
        embedding = stability_row_normalize(embedding);
        rng(cfg.seed+1709*c+31*size(B,2),'twister');
        [labels,centres,kmeansObjective] = stability_local_kmeans( ...
            embedding,c,cfg.view_kmeans_restarts, ...
            cfg.view_kmeans_max_iter);
        labels = stability_order_labels_by_horizon(labels,horizons);
        metrics = stability_view_cluster_metrics( ...
            Wnorm,embedding,labels,horizons);

        item = stability_empty_cluster_result();
        item.cluster_count = c;
        item.labels = labels;
        item.embedding = embedding;
        item.centres = centres;
        item.kmeans_objective = kmeansObjective;
        item.normalized_cut = metrics.normalized_cut;
        item.modularity = metrics.modularity;
        item.embedding_silhouette = metrics.embedding_silhouette;
        item.adjacent_horizon_same_cluster_fraction = ...
            metrics.adjacent_horizon_same_cluster_fraction;
        item.horizon_order_transition_count = ...
            metrics.horizon_order_transition_count;
        item.cluster_fragmentation = metrics.cluster_fragmentation;
        item.mean_within_cluster_log_span = ...
            metrics.mean_within_cluster_log_span;
        item.max_within_cluster_log_span = ...
            metrics.max_within_cluster_log_span;
        sweep(q) = item;
    end

    maxGapCount = min(cfg.max_eigengap_cluster_count,K-1);
    gapCounts = counts(counts<=maxGapCount);
    if isempty(gapCounts)
        gapCounts = counts;
    end
    eigengaps = nan(size(gapCounts));
    for q = 1:numel(gapCounts)
        c = gapCounts(q);
        eigengaps(q) = eigenvalues(c+1)-eigenvalues(c);
    end
    [~,gapIndex] = max(eigengaps);
    eigengapCount = gapCounts(gapIndex);

    rankMatched = min(max(size(B,2),2),K-1);
    if ~ismember(rankMatched,counts)
        [~,nearest] = min(abs(counts-rankMatched));
        rankMatched = counts(nearest);
    end

    clusterInfo = struct();
    clusterInfo.normalized_adjacency = Wnorm;
    clusterInfo.degree = degree;
    clusterInfo.normalized_laplacian = Lsym;
    clusterInfo.eigenvalues = eigenvalues;
    clusterInfo.eigenvectors = V;
    clusterInfo.eigengap_counts = gapCounts;
    clusterInfo.eigengaps = eigengaps;
    clusterInfo.eigengap_suggested_count = eigengapCount;
    clusterInfo.rank_matched_count = rankMatched;
    clusterInfo.sweep = sweep;
end


function item = stability_empty_cluster_result()
    item = struct('cluster_count',NaN,'labels',[], ...
        'embedding',[],'centres',[],'kmeans_objective',NaN, ...
        'normalized_cut',NaN,'modularity',NaN, ...
        'embedding_silhouette',NaN, ...
        'adjacent_horizon_same_cluster_fraction',NaN, ...
        'horizon_order_transition_count',NaN, ...
        'cluster_fragmentation',NaN, ...
        'mean_within_cluster_log_span',NaN, ...
        'max_within_cluster_log_span',NaN);
end


function X = stability_row_normalize(X)
    rowNorm = sqrt(sum(X.^2,2));
    X = X./max(rowNorm,eps);
end


function [bestLabels,bestCentres,bestObjective] = ...
        stability_local_kmeans(X,c,nStarts,maxIter)

    n = size(X,1);
    bestObjective = inf;
    bestLabels = [];
    bestCentres = [];

    for start = 1:nStarts
        centres = stability_initialize_centres(X,c,start);
        labels = ones(n,1);
        previousLabels = zeros(n,1);

        for iter = 1:maxIter
            D2 = stability_squared_distance(X,centres);
            [~,labels] = min(D2,[],2);
            if isequal(labels,previousLabels)
                break;
            end
            previousLabels = labels;

            for g = 1:c
                members = X(labels==g,:);
                if isempty(members)
                    nearest = min(D2,[],2);
                    [~,replacement] = max(nearest);
                    centres(g,:) = X(replacement,:);
                    labels(replacement) = g;
                else
                    centres(g,:) = mean(members,1);
                end
            end
        end

        D2 = stability_squared_distance(X,centres);
        linearIndex = sub2ind(size(D2),(1:n).',labels);
        objective = sum(D2(linearIndex));
        if objective < bestObjective
            bestObjective = objective;
            bestLabels = labels;
            bestCentres = centres;
        end
    end
end


function centres = stability_initialize_centres(X,c,start)
    n = size(X,1);
    centres = zeros(c,size(X,2));
    if start==1
        [~,first] = max(sum((X-mean(X,1)).^2,2));
    else
        first = randi(n);
    end
    chosen = false(n,1);
    chosen(first) = true;
    centres(1,:) = X(first,:);
    minD2 = sum((X-centres(1,:)).^2,2);

    for g = 2:c
        minD2(chosen) = -inf;
        [~,next] = max(minD2);
        if ~isfinite(minD2(next))
            next = find(~chosen,1);
        end
        chosen(next) = true;
        centres(g,:) = X(next,:);
        current = sum((X-centres(g,:)).^2,2);
        minD2 = min(minD2,current);
    end
end


function D2 = stability_squared_distance(X,C)
    D2 = sum(X.^2,2)+sum(C.^2,2).'-2*(X*C.');
    D2 = max(D2,0);
end


function labels = stability_order_labels_by_horizon(labels,horizons)
    uniqueLabels = unique(labels(:)).';
    centres = zeros(size(uniqueLabels));
    for q = 1:numel(uniqueLabels)
        centres(q) = mean(log(horizons(labels==uniqueLabels(q))));
    end
    [~,order] = sort(centres,'ascend');
    relabeled = zeros(size(labels));
    for q = 1:numel(order)
        relabeled(labels==uniqueLabels(order(q))) = q;
    end
    labels = relabeled;
end


function metrics = stability_view_cluster_metrics(W,embedding,labels,horizons)
    uniqueLabels = unique(labels(:)).';
    c = numel(uniqueLabels);
    degree = sum(W,2);
    totalVolume = sum(degree);
    ncut = 0;
    spans = zeros(c,1);

    [~,horizonOrder] = sort(horizons,'ascend');
    labelsSorted = labels(horizonOrder);
    segments = 0;

    for q = 1:c
        mask = labels==uniqueLabels(q);
        cutWeight = sum(W(mask,~mask),'all');
        volume = sum(degree(mask));
        ncut = ncut + cutWeight/max(volume,eps);
        localH = horizons(mask);
        spans(q) = log(max(localH)/min(localH));

        localMask = labelsSorted==uniqueLabels(q);
        padded = [false;localMask(:);false];
        segments = segments + nnz(diff(padded)==1);
    end

    if totalVolume<=eps
        modularity = NaN;
    else
        expected = degree*degree.'/totalVolume;
        same = labels==labels.';
        modularity = sum((W-expected).*same,'all')/totalVolume;
    end

    silhouetteValue = stability_mean_silhouette(embedding,labels);
    adjacentSame = mean(labelsSorted(1:end-1)==labelsSorted(2:end));
    transitionCount = nnz(labelsSorted(1:end-1)~=labelsSorted(2:end));

    metrics = struct();
    metrics.normalized_cut = ncut;
    metrics.modularity = modularity;
    metrics.embedding_silhouette = silhouetteValue;
    metrics.adjacent_horizon_same_cluster_fraction = adjacentSame;
    metrics.horizon_order_transition_count = transitionCount;
    metrics.cluster_fragmentation = segments-c;
    metrics.mean_within_cluster_log_span = mean(spans);
    metrics.max_within_cluster_log_span = max(spans);
end


function value = stability_mean_silhouette(X,labels)
    n = size(X,1);
    D = sqrt(stability_squared_distance(X,X));
    groups = unique(labels(:)).';
    if numel(groups)<2
        value = NaN;
        return;
    end

    s = zeros(n,1);
    for i = 1:n
        own = labels==labels(i);
        own(i) = false;
        if any(own)
            a = mean(D(i,own));
        else
            a = 0;
        end
        b = inf;
        for g = groups
            if g==labels(i)
                continue;
            end
            b = min(b,mean(D(i,labels==g)));
        end
        denom = max(a,b);
        if denom<=eps
            s(i) = 0;
        else
            s(i) = (b-a)/denom;
        end
    end
    value = mean(s);
end


function T = stability_cluster_metrics_table(clusterInfo)
    n = numel(clusterInfo.sweep);
    c = zeros(n,1);
    ncut = zeros(n,1);
    modularity = zeros(n,1);
    silhouette = zeros(n,1);
    adjacent = zeros(n,1);
    transitions = zeros(n,1);
    fragmentation = zeros(n,1);
    meanSpan = zeros(n,1);
    maxSpan = zeros(n,1);
    objective = zeros(n,1);

    for q = 1:n
        item = clusterInfo.sweep(q);
        c(q) = item.cluster_count;
        ncut(q) = item.normalized_cut;
        modularity(q) = item.modularity;
        silhouette(q) = item.embedding_silhouette;
        adjacent(q) = item.adjacent_horizon_same_cluster_fraction;
        transitions(q) = item.horizon_order_transition_count;
        fragmentation(q) = item.cluster_fragmentation;
        meanSpan(q) = item.mean_within_cluster_log_span;
        maxSpan(q) = item.max_within_cluster_log_span;
        objective(q) = item.kmeans_objective;
    end

    T = table(c,ncut,modularity,silhouette,adjacent,transitions, ...
        fragmentation,meanSpan,maxSpan,objective, ...
        'VariableNames',{'cluster_count','normalized_cut', ...
        'modularity','embedding_silhouette', ...
        'adjacent_horizon_same_cluster_fraction', ...
        'horizon_order_transition_count','cluster_fragmentation', ...
        'mean_within_cluster_log_span','max_within_cluster_log_span', ...
        'kmeans_objective'});
end


function stability_write_view_matrix(filePath,M,horizons)
    K = numel(horizons);
    variableNames = cell(1,K+1);
    variableNames{1} = 'horizon_days';
    for k = 1:K
        variableNames{k+1} = sprintf('h_%03dd',round(horizons(k)));
    end
    T = array2table([horizons(:),M], ...
        'VariableNames',variableNames);
    writetable(T,filePath);
end


function stability_plot_view_graph_components( ...
        Gview,Gcross,Gcopy,B,sA,sC,horizons,savedError,blockError, ...
        candidateId,filePath,cfg)

    f = finance_new_figure([1750,1050]);
    tl = tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');

    stability_view_heatmap(nexttile(tl,1),Gview,horizons, ...
        'Total induced view graph');
    stability_view_heatmap(nexttile(tl,2),Gcross,horizons, ...
        'Cross-stock contribution');
    stability_view_heatmap(nexttile(tl,3),Gcopy,horizons, ...
        'Same-stock copy contribution');

    ax = nexttile(tl,4);
    bar(ax,1:numel(sA),[2*sA,sC],'stacked');
    grid(ax,'on'); box(ax,'on');
    xlabel(ax,'mode'); ylabel(ax,'mass entering G_{view}');
    legend(ax,{'2s_a','s_c'},'Location','best');
    title(ax,'Modal mass decomposition');

    ax = nexttile(tl,5);
    plot(ax,horizons,B,'-o','LineWidth',1.5,'MarkerSize',5);
    set(ax,'XScale','log'); grid(ax,'on'); box(ax,'on');
    xlabel(ax,'return horizon (days)'); ylabel(ax,'B weight');
    title(ax,'Learned modal horizon profiles');
    legend(ax,arrayfun(@(m)sprintf('mode %d',m), ...
        1:size(B,2),'UniformOutput',false),'Location','best');

    ax = nexttile(tl,6);
    axis(ax,'off');
    text(ax,0.02,0.80,sprintf('saved G_{view} relative error: %.3e', ...
        savedError),'FontSize',12);
    text(ax,0.02,0.62,sprintf('supra block-sum relative error: %.3e', ...
        blockError),'FontSize',12);
    text(ax,0.02,0.44,sprintf('total view-graph mass: %.4g', ...
        sum(Gview(:))/2),'FontSize',12);
    text(ax,0.02,0.26,sprintf('cross contribution fraction: %.3f', ...
        sum(Gcross(:))/max(sum(Gview(:)),eps)),'FontSize',12);
    text(ax,0.02,0.08,sprintf('copy contribution fraction: %.3f', ...
        sum(Gcopy(:))/max(sum(Gview(:)),eps)),'FontSize',12);

    sgtitle(tl,sprintf('%s: induced graph over return horizons', ...
        candidateId),'Interpreter','none');
    finance_export_figure(f,filePath,cfg);
end


function stability_plot_view_spectral_diagnostics( ...
        clusterInfo,candidateId,filePath,cfg)

    T = stability_cluster_metrics_table(clusterInfo);
    f = finance_new_figure([1600,900]);
    tl = tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');

    ax = nexttile(tl,1);
    plot(ax,0:numel(clusterInfo.eigenvalues)-1, ...
        clusterInfo.eigenvalues,'-o','LineWidth',1.8);
    grid(ax,'on'); box(ax,'on');
    xlabel(ax,'eigenvalue index'); ylabel(ax,'\lambda');
    title(ax,'Normalized-Laplacian spectrum');

    ax = nexttile(tl,2);
    bar(ax,clusterInfo.eigengap_counts,clusterInfo.eigengaps);
    grid(ax,'on'); box(ax,'on');
    xlabel(ax,'candidate cluster count c');
    ylabel(ax,'\lambda_{c+1}-\lambda_c');
    title(ax,sprintf('Eigengap diagnostic: c=%d', ...
        clusterInfo.eigengap_suggested_count));

    ax = nexttile(tl,3);
    yyaxis(ax,'left');
    plot(ax,T.cluster_count,T.normalized_cut,'-o','LineWidth',1.6);
    ylabel(ax,'normalized cut (lower)');
    yyaxis(ax,'right');
    plot(ax,T.cluster_count,T.modularity,'-s','LineWidth',1.6);
    ylabel(ax,'modularity (higher)');
    grid(ax,'on'); box(ax,'on');
    xlabel(ax,'cluster count');
    title(ax,'Graph partition quality');

    ax = nexttile(tl,4);
    plot(ax,T.cluster_count,T.embedding_silhouette, ...
        '-o','LineWidth',1.6); hold(ax,'on');
    plot(ax,T.cluster_count, ...
        T.adjacent_horizon_same_cluster_fraction, ...
        '-s','LineWidth',1.6);
    plot(ax,T.cluster_count,1./(1+T.cluster_fragmentation), ...
        '-d','LineWidth',1.6);
    grid(ax,'on'); box(ax,'on');
    xlabel(ax,'cluster count'); ylabel(ax,'diagnostic value');
    legend(ax,{'embedding silhouette','adjacent horizons co-clustered', ...
        '1/(1+fragmentation)'},'Location','best');
    title(ax,'Horizon-cluster coherence');

    sgtitle(tl,sprintf('%s: view-graph clustering diagnostics', ...
        candidateId),'Interpreter','none');
    finance_export_figure(f,filePath,cfg);
end


function stability_plot_all_view_clusterings( ...
        W,horizons,clusterInfo,candidateId,filePath,cfg)

    nPlots = numel(clusterInfo.sweep);
    f = finance_new_figure([600*nPlots,600]);
    tl = tiledlayout(f,1,nPlots,'TileSpacing','compact','Padding','compact');
    for q = 1:nPlots
        item = clusterInfo.sweep(q);
        ax = nexttile(tl,q);
        stability_plot_cluster_reordered_heatmap( ...
            ax,W,horizons,item.labels);
        title(ax,sprintf('c=%d | Ncut %.3f | Q %.3f', ...
            item.cluster_count,item.normalized_cut,item.modularity));
    end
    sgtitle(tl,sprintf('%s: all requested horizon clusterings', ...
        candidateId),'Interpreter','none');
    finance_export_figure(f,filePath,cfg);
end


function stability_plot_main_view_clustering( ...
        W,horizons,item,candidateId,label,filePath,cfg)

    f = finance_new_figure([1500,650]);
    tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');

    ax = nexttile(tl,1);
    stability_plot_cluster_reordered_heatmap(ax,W,horizons,item.labels);
    title(ax,sprintf('%s clustering, c=%d',label,item.cluster_count));

    ax = nexttile(tl,2);
    [sortedH,order] = sort(horizons);
    labels = item.labels(order);
    imagesc(ax,1:numel(sortedH),1,labels.');
    colormap(ax,lines(item.cluster_count));
    set(ax,'YTick',[],'XTick',1:numel(sortedH), ...
        'XTickLabel',arrayfun(@(h)sprintf('%dd',h),sortedH, ...
        'UniformOutput',false));
    xlabel(ax,'return horizon');
    title(ax,sprintf(['ordered clusters | transitions %d | ', ...
        'silhouette %.3f'],item.horizon_order_transition_count, ...
        item.embedding_silhouette));

    sgtitle(tl,sprintf('%s: %s view-graph partition', ...
        candidateId,label),'Interpreter','none');
    finance_export_figure(f,filePath,cfg);
end


function stability_plot_modal_view_graphs( ...
        modalGraphs,B,sA,sC,horizons,candidateId,filePath,cfg)

    r = size(modalGraphs,3);
    f = finance_new_figure([620*r,620]);
    tl = tiledlayout(f,1,r,'TileSpacing','compact','Padding','compact');
    for m = 1:r
        ax = nexttile(tl,m);
        stability_view_heatmap(ax,modalGraphs(:,:,m),horizons, ...
            sprintf('mode %d | 2s_a+s_c=%.3g',m,2*sA(m)+sC(m)));
    end
    sgtitle(tl,sprintf('%s: modal contributions to G_{view}', ...
        candidateId),'Interpreter','none');
    finance_export_figure(f,filePath,cfg);
end


function stability_view_heatmap(ax,M,horizons,titleText)
    imagesc(ax,M);
    axis(ax,'image'); set(ax,'YDir','normal'); colorbar(ax);
    labels = arrayfun(@(h)sprintf('%dd',h),horizons, ...
        'UniformOutput',false);
    set(ax,'XTick',1:numel(horizons),'YTick',1:numel(horizons), ...
        'XTickLabel',labels,'YTickLabel',labels, ...
        'XTickLabelRotation',35);
    xlabel(ax,'horizon'); ylabel(ax,'horizon');
    title(ax,titleText);
end


function stability_plot_cluster_reordered_heatmap(ax,W,horizons,labels)
    [~,order] = sortrows([labels(:),log(horizons(:))],[1 2]);
    imagesc(ax,W(order,order));
    axis(ax,'image'); set(ax,'YDir','normal'); colorbar(ax);
    tickLabels = arrayfun(@(h,c)sprintf('%dd|c%d',h,c), ...
        horizons(order),labels(order),'UniformOutput',false);
    set(ax,'XTick',1:numel(order),'YTick',1:numel(order), ...
        'XTickLabel',tickLabels,'YTickLabel',tickLabels, ...
        'XTickLabelRotation',45);
    xlabel(ax,'cluster-ordered horizons');
    ylabel(ax,'cluster-ordered horizons');
end


function item = stability_get_cluster_item(clusterInfo,count)
    index = find([clusterInfo.sweep.cluster_count]==count,1);
    if isempty(index)
        error('Requested cluster count %d was not computed.',count);
    end
    item = clusterInfo.sweep(index);
end


function M = stability_symmetrize_nonnegative(M)
    M = max(0,0.5*(double(M)+double(M).'));
end


function M = stability_remove_diagonal(M)
    M(1:size(M,1)+1:end) = 0;
end


function stability_save_rank_representative_comparison( ...
        finalRankResults,consensusTable,data,outDir,modelRank,cfg)

    n = numel(finalRankResults);
    if n==0
        return;
    end

    f = finance_new_figure([620*n,1080]);
    tl = tiledlayout(f,2,n,'TileSpacing','compact','Padding','compact');
    for q = 1:n
        R = finalRankResults{q};
        ax = nexttile(tl,q);
        stability_view_heatmap(ax,R.view_graph.G_view,data.horizons, ...
            sprintf('representative %d: G_{view}',q));

        ax = nexttile(tl,n+q);
        item = stability_get_cluster_item( ...
            R.view_graph.clustering, ...
            R.view_graph.clustering.rank_matched_count);
        [sortedH,order] = sort(data.horizons);
        imagesc(ax,1:numel(sortedH),1,item.labels(order).');
        colormap(ax,lines(item.cluster_count));
        set(ax,'YTick',[],'XTick',1:numel(sortedH), ...
            'XTickLabel',arrayfun(@(h)sprintf('%dd',h),sortedH, ...
            'UniformOutput',false));
        title(ax,sprintf('rank-matched clusters, c=%d', ...
            item.cluster_count));
    end
    sgtitle(tl,sprintf('r=%d: three consensus-representative fits', ...
        modelRank));
    finance_export_figure(f,fullfile(outDir, ...
        'representative_view_graph_comparison.png'),cfg);

    nConsensus = min(cfg.consensus_relation_count,height(consensusTable));
    C = consensusTable(1:nConsensus,:);
    labels = strings(nConsensus,1);
    for q = 1:nConsensus
        labels(q) = sprintf('%s--%s | %dd--%dd', ...
            C.stock_a_ticker(q),C.stock_b_ticker(q), ...
            C.horizon_low_days(q),C.horizon_high_days(q));
    end

    matchMatrix = zeros(nConsensus,n);
    for q = 1:n
        T = finalRankResults{q}.consensus_match;
        matchMatrix(:,q) = double(T.present_in_final_refit);
    end
    matchTable = array2table(matchMatrix, ...
        'VariableNames',arrayfun(@(q)sprintf('representative_%d',q), ...
        1:n,'UniformOutput',false));
    matchTable.relation = labels;
    matchTable = movevars(matchTable,'relation','Before',1);
    writetable(matchTable,fullfile(outDir, ...
        'representative_consensus_presence.csv'));
end


function stability_write_readme(filePath,dataFile,outputRoot,cfg, ...
        inputInfo,screeningTable,consensusResults)

    fid = fopen(filePath,'w');
    if fid<0
        return;
    end
    cleanup = onCleanup(@()fclose(fid)); %#ok<NASGU>

    fprintf(fid,'MOSAIC FINANCE STABILITY-CONSENSUS RESULTS\n');
    fprintf(fid,'==========================================\n\n');
    fprintf(fid,'Input MAT-file: %s\n',dataFile);
    fprintf(fid,'Output folder : %s\n\n',outputRoot);
    fprintf(fid,'Data: n=%d stocks, K=%d horizons, p=%d coordinates.\n', ...
        inputInfo.n,inputInfo.K,inputInfo.p);
    fprintf(fid,'Horizons: %s trading days.\n\n', ...
        mat2str(inputInfo.horizons));
    fprintf(fid,'Ranks screened: %s.\n',mat2str(cfg.rank_candidates));
    fprintf(fid,'Profiles per rank: %d.\n',cfg.n_profiles_per_rank);
    fprintf(fid,'Total screening candidates: %d.\n',height(screeningTable));
    fprintf(fid,'Top relations retained per candidate: %d.\n', ...
        cfg.top_relations_per_candidate);
    fprintf(fid,'Consensus relations reported per rank: %d.\n', ...
        cfg.consensus_relation_count);
    fprintf(fid,'Representative fits retained per rank: %d.\n\n', ...
        cfg.n_representatives_per_rank);

    fprintf(fid,['SELECTION MEANING\nEach usable run contributes its ', ...
        'global top %d cross-stock/cross-horizon relations, with no ', ...
        'horizon-pair quota. Relation weights are normalized by the ', ...
        'strongest relation in the same run. A relation is ranked by its ', ...
        'average normalized strength across all usable runs, assigning ', ...
        'zero when absent. The top %d relations form the consensus set. ', ...
        'Each candidate is scored by its mean normalized weight over that ', ...
        'set. Three representatives are chosen from the %.1f%% score ', ...
        'plateau, preferring distinct profile families. The selected fits ', ...
        'are not oracle-optimal or predictive-best because no graph truth ', ...
        'is available. Sector labels are used only after fitting.\n\n'], ...
        cfg.top_relations_per_candidate,cfg.consensus_relation_count, ...
        100*cfg.representative_plateau_relative_tolerance);

    for q = 1:numel(consensusResults)
        C = consensusResults{q};
        fprintf(fid,'Rank r=%d: %d usable candidates entered consensus.\n', ...
            C.rank,numel(C.usable_candidate_indices));
        fprintf(fid,'Selected candidate indices: %s.\n', ...
            mat2str(C.selected_table.candidate_index.'));
    end

    fprintf(fid,'\nOUTPUT STRUCTURE\n');
    fprintf(fid,'00_input_diagnostics       prepared-signal checks\n');
    fprintf(fid,'01_screening_fits          compact fit for all screened candidates\n');
    fprintf(fid,'02_screening_tables        top cross relations for every candidate\n');
    fprintf(fid,'03_consensus               per-rank relation and candidate stability\n');
    fprintf(fid,'04_selected_representatives detailed top-1 fits for r=3 and r=4\n');
    fprintf(fid,'05_summary                 configuration and diagnostic catalogues\n');
end


function stability_write_error_report(filePath,ME,candidateId)
    fid = fopen(filePath,'w');
    if fid<0
        return;
    end
    cleanup = onCleanup(@()fclose(fid)); %#ok<NASGU>
    fprintf(fid,'Candidate: %s\n',candidateId);
    fprintf(fid,'Identifier: %s\n',ME.identifier);
    fprintf(fid,'Message: %s\n\n',ME.message);
    for q = 1:numel(ME.stack)
        fprintf(fid,'%s, line %d\n',ME.stack(q).name,ME.stack(q).line);
    end
end



%% ========================================================================
%% ORIGINAL FINANCE ANALYSIS + EMBEDDED MOSAIC PALM SOLVER
%% ========================================================================
function [data,Y,observed,inputInfo] = finance_prepare_observed(data,cfg)
    required = {'X','horizons'};
    for q = 1:numel(required)
        if ~isfield(data,required{q})
            error('data.%s is required.',required{q});
        end
    end
    if ~iscell(data.X) || isempty(data.X)
        error('data.X must be a nonempty cell array.');
    end

    K = numel(data.X);
    horizons = double(data.horizons(:)');
    if numel(horizons) ~= K
        error('data.horizons must contain one horizon per view.');
    end
    if any(horizons <= 0) || any(diff(horizons) <= 0)
        error('data.horizons must be strictly increasing positive values.');
    end

    n = size(data.X{1},1);
    p = size(data.X{1},2);
    for k = 1:K
        if ~isequal(size(data.X{k}),[n,p])
            error('Every data.X{k} must have the same n-by-p dimensions.');
        end
        if any(~isfinite(data.X{k}(:)))
            error('data.X{%d} contains non-finite values.',k);
        end
    end

    if isfield(data,'tickers')
        tickers = string(data.tickers(:));
    elseif isfield(data,'node_names')
        tickers = string(data.node_names(:));
    else
        tickers = "stock_" + string((1:n)');
    end
    if numel(tickers) ~= n
        error('The number of ticker/node names must equal n.');
    end

    if isfield(data,'sector_labels')
        sectorLabels = string(data.sector_labels(:));
    else
        sectorLabels = repmat("Unknown",n,1);
    end
    if numel(sectorLabels) ~= n
        error('data.sector_labels must contain n entries.');
    end
    sectorLabels(strlength(strtrim(sectorLabels))==0) = "Unknown";

    if isfield(data,'sector_known')
        sectorKnown = logical(data.sector_known(:));
    else
        sectorKnown = sectorLabels ~= "Unknown";
    end
    if numel(sectorKnown) ~= n
        sectorKnown = sectorLabels ~= "Unknown";
    end

    if isfield(data,'view_names') && numel(data.view_names)==K
        viewNames = cellstr(string(data.view_names(:)));
    else
        viewNames = arrayfun(@(h)sprintf('%d-day',h),horizons, ...
            'UniformOutput',false);
    end

    [sectorOrder,sectorNames,sectorBoundaries,sectorIndex] = ...
        finance_sector_order(sectorLabels);

    data.num_views = K;
    data.num_stocks = n;
    data.signal_dimension = p;
    data.horizons = horizons;
    data.tickers = tickers;
    data.node_names = tickers;
    data.sector_labels = sectorLabels;
    data.sector_known = sectorKnown;
    data.sector_order = sectorOrder;
    data.sector_names = sectorNames;
    data.sector_boundaries = sectorBoundaries;
    data.sector_index = sectorIndex;
    data.view_names = viewNames;

    Y = vertcat(data.X{:});
    if cfg.center_each_signal
        Y = Y - mean(Y,2);
    end
    if cfg.l2_normalize_each_signal
        rowNorms = sqrt(sum(Y.^2,2));
        Y = Y ./ max(rowNorms,eps);
    end

    sampleCount = size(Y,2);
    d = sum(Y.^2,2)/sampleCount;
    Sraw = max(d+d'-2*(Y*Y'/sampleCount),0);
    Sraw = 0.5*(Sraw+Sraw');
    Sraw(1:size(Sraw,1)+1:end) = 0;
    upperValues = Sraw(triu(true(size(Sraw)),1));
    positiveValues = upperValues(upperValues>0 & isfinite(upperValues));
    if isempty(positiveValues)
        error('All empirical pairwise signal distances are zero.');
    end
    empiricalScale = median(positiveValues);
    Snorm = Sraw/empiricalScale;

    observed = struct();
    observed.n = n;
    observed.K = K;
    observed.Y = Y;
    observed.S = Snorm;
    observed.cost_scale = 1;

    rawCopyDistance = zeros(K,K);
    rawCrossDistance = zeros(K,K);
    rawBlockDistance = zeros(K,K);
    offMask = ~eye(n);
    for k = 1:K
        Ik = (k-1)*n+(1:n);
        for ell = 1:K
            Il = (ell-1)*n+(1:n);
            block = Sraw(Ik,Il);
            rawCopyDistance(k,ell) = mean(diag(block));
            rawCrossDistance(k,ell) = mean(block(offMask));
            rawBlockDistance(k,ell) = mean(block(:));
        end
    end

    positiveCopy = rawCopyDistance(triu(true(K),1));
    positiveCopy = positiveCopy(positiveCopy>0 & isfinite(positiveCopy));
    if isempty(positiveCopy)
        copyScale = 1;
    else
        copyScale = median(positiveCopy);
    end
    rawViewSimilarity = exp(-rawCopyDistance/max(copyScale,eps));
    rawViewSimilarity(1:K+1:end) = 0;

    inputInfo = struct();
    inputInfo.n = n;
    inputInfo.K = K;
    inputInfo.p = p;
    inputInfo.horizons = horizons;
    inputInfo.view_names = viewNames;
    inputInfo.tickers = tickers;
    inputInfo.sector_labels = sectorLabels;
    inputInfo.sector_known = sectorKnown;
    inputInfo.sector_order = sectorOrder;
    inputInfo.sector_names = sectorNames;
    inputInfo.sector_boundaries = sectorBoundaries;
    inputInfo.sector_index = sectorIndex;
    inputInfo.empirical_distance_scale = empiricalScale;
    inputInfo.raw_distance_matrix = Sraw;
    inputInfo.normalized_distance_matrix = Snorm;
    inputInfo.raw_copy_distance_by_view = rawCopyDistance;
    inputInfo.raw_crossnode_distance_by_view = rawCrossDistance;
    inputInfo.raw_block_distance_by_view = rawBlockDistance;
    inputInfo.raw_view_similarity = rawViewSimilarity;
    inputInfo.center_each_signal = cfg.center_each_signal;
    inputInfo.l2_normalize_each_signal = cfg.l2_normalize_each_signal;
end


function finance_save_input_diagnostics(data,Y,inputInfo,outDir,cfg)
    n = inputInfo.n;
    K = inputInfo.K;
    horizons = inputInfo.horizons;
    labels = finance_view_labels(data);

    % Sector-sorted input-signal heatmaps.
    nCols = ceil(sqrt(K));
    nRows = ceil(K/nCols);
    f = finance_new_figure([540*nCols,430*nRows]);
    tl = tiledlayout(f,nRows,nCols,'TileSpacing','compact','Padding','compact');
    lim = quantile(abs(Y(:)),0.995);
    lim = max(lim,eps);
    for k = 1:K
        ax = nexttile(tl,k);
        imagesc(ax,data.X{k}(data.sector_order,:),[-lim lim]);
        colorbar(ax); xlabel(ax,'aligned ending-date coordinate');
        ylabel(ax,'stocks sorted by GICS sector');
        title(ax,labels{k});
        finance_draw_sector_boundaries(ax,data.sector_boundaries,'y');
    end
    sgtitle(tl,'Prepared standardized return signals (sector-sorted rows)');
    finance_export_figure(f,fullfile(outDir,'01_input_signal_heatmaps.png'),cfg);

    % Raw cross-horizon distance/similarity structure.
    f = finance_new_figure([1550,520]);
    tl = tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
    ax = nexttile(tl,1);
    imagesc(ax,inputInfo.raw_copy_distance_by_view);
    axis(ax,'image'); colorbar(ax); set(ax,'YDir','normal');
    finance_set_view_ticks(ax,labels);
    title(ax,'Mean same-stock squared distance');

    ax = nexttile(tl,2);
    imagesc(ax,inputInfo.raw_crossnode_distance_by_view);
    axis(ax,'image'); colorbar(ax); set(ax,'YDir','normal');
    finance_set_view_ticks(ax,labels);
    title(ax,'Mean different-stock squared distance');

    ax = nexttile(tl,3);
    imagesc(ax,inputInfo.raw_view_similarity,[0 1]);
    axis(ax,'image'); colorbar(ax); set(ax,'YDir','normal');
    finance_set_view_ticks(ax,labels);
    title(ax,'Raw same-stock horizon similarity');
    sgtitle(tl,'Input cross-horizon diagnostics');
    finance_export_figure(f,fullfile(outDir,'02_input_horizon_diagnostics.png'),cfg);

    % Preliminary finance diagnostics retained from the dataset builder.
    if isfield(data,'view_diagnostics') && istable(data.view_diagnostics)
        T = data.view_diagnostics;
        f = finance_new_figure([1450,850]);
        tl = tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
        ax = nexttile(tl,1);
        plot(ax,horizons,T.median_pairwise_distance,'-o','LineWidth',1.7);
        set(ax,'XScale','log'); grid(ax,'on'); box(ax,'on');
        xlabel(ax,'return horizon (days)'); ylabel(ax,'median pairwise distance');
        title(ax,'Within-view signal dispersion');
        ax = nexttile(tl,2);
        plot(ax,horizons,T.mean_lag1_autocorr,'-o','LineWidth',1.7);
        set(ax,'XScale','log'); grid(ax,'on'); box(ax,'on');
        xlabel(ax,'return horizon (days)'); ylabel(ax,'mean lag-1 correlation');
        title(ax,'Coordinate overlap / temporal dependence');
        ax = nexttile(tl,3);
        plot(ax,horizons,T.mean_within_sector_corr,'-o','LineWidth',1.7); hold(ax,'on');
        plot(ax,horizons,T.mean_between_sector_corr,'-s','LineWidth',1.7);
        set(ax,'XScale','log'); grid(ax,'on'); box(ax,'on');
        legend(ax,{'within sector','between sector'},'Location','best');
        xlabel(ax,'return horizon (days)'); ylabel(ax,'mean signal correlation');
        title(ax,'Sector signal structure');
        ax = nexttile(tl,4);
        plot(ax,horizons,T.within_minus_between_corr,'-o', ...
            'LineWidth',1.7,'MarkerSize',6);
        set(ax,'XScale','log'); grid(ax,'on'); box(ax,'on');
        yline(ax,0,'--','LineWidth',0.8);
        xlabel(ax,'return horizon (days)'); ylabel(ax,'within - between correlation');
        title(ax,'Sector separation in input signals');
        sgtitle(tl,'Prepared finance-signal diagnostics');
        finance_export_figure(f,fullfile(outDir,'03_input_finance_diagnostics.png'),cfg);
    end

    inputSummary = struct();
    inputSummary.Y_min = min(Y(:));
    inputSummary.Y_max = max(Y(:));
    inputSummary.Y_mean = mean(Y(:));
    inputSummary.Y_rms = sqrt(mean(Y(:).^2));
    inputSummary.row_norm_min = min(sqrt(sum(Y.^2,2)));
    inputSummary.row_norm_max = max(sqrt(sum(Y.^2,2)));
    inputSummary.empirical_distance_scale = inputInfo.empirical_distance_scale;
    inputSummary.n = n;
    inputSummary.K = K;
    inputSummary.p = inputInfo.p;
    save(fullfile(outDir,'input_diagnostics.mat'),'inputSummary','inputInfo');
end


function analysis = finance_analyze_fit(fit,data,inputInfo,cfg)
    n = inputInfo.n;
    K = inputInfo.K;
    horizons = inputInfo.horizons;
    r = size(fit.B,2);

    within = zeros(n,n,K);
    copyMean = zeros(K,K);
    crossMean = zeros(K,K);
    blockMean = zeros(K,K);
    crossBlocks = cell(K,K);
    copyVectors = cell(K,K);

    for k = 1:K
        for m = 1:r
            within(:,:,k) = within(:,:,k) + ...
                fit.B(k,m)^2*fit.A_state(:,:,m);
        end
        for ell = 1:K
            if k==ell
                C = within(:,:,k);
            else
                C = zeros(n);
                for m = 1:r
                    C = C + fit.B(k,m)*fit.B(ell,m)* ...
                        (fit.A_state(:,:,m)+diag(fit.Gamma(:,m)));
                end
            end
            copyMean(k,ell) = mean(diag(C));
            offMask = ~eye(n);
            crossMean(k,ell) = mean(C(offMask));
            blockMean(k,ell) = mean(C(:));
            if k~=ell
                Ccross = C;
                Ccross(1:n+1:end) = 0;
                crossBlocks{k,ell} = Ccross;
                copyVectors{k,ell} = diag(C);
            end
        end
    end
    copyMean(1:K+1:end) = 0;

    [withinValues,crossValues,copyValues] = ...
        finance_unique_component_values(fit.A_supra,n,K);
    totalMass = sum(withinValues)+sum(crossValues)+sum(copyValues);
    if totalMass<=0, totalMass=1; end

    threshold = cfg.display_relative_threshold;
    withinDensity = finance_relative_density(withinValues,threshold);
    crossDensity = finance_relative_density(crossValues,threshold);
    copyDensity = finance_relative_density(copyValues,threshold);

    [viewCorr,viewCount] = finance_log_horizon_correlation(fit.view_graph,horizons);
    [copyCorr,copyCount] = finance_log_horizon_correlation(copyMean,horizons);
    [crossCorr,crossCount] = finance_log_horizon_correlation(crossMean,horizons);

    Bcosine = finance_max_column_cosine(fit.B);
    Bentropy = zeros(1,r);
    Broughness = zeros(1,r);
    Bcentres = zeros(1,r);
    logH = log(horizons(:));
    for m = 1:r
        b = max(fit.B(:,m),eps);
        Bentropy(m) = -sum(b.*log(b))/log(K);
        Broughness(m) = sum(diff(fit.B(:,m)).^2);
        Bcentres(m) = exp(sum(fit.B(:,m).*logH));
    end

    strongestCrossPair = finance_strongest_upper_pair(crossMean);
    strongestCopyPair = finance_strongest_upper_pair(copyMean);
    strongestViewPair = finance_strongest_upper_pair(fit.view_graph);

    withinStats = finance_within_sector_stats(within,data);
    crossStats = finance_cross_sector_stats(crossBlocks,data);
    modeStats = finance_mode_sector_stats(fit.A_state,data);
    stockStrength = finance_stock_strength_table(fit.A_supra,data);
    sectorStrength = finance_sector_strength_table(stockStrength,data);

    analysis = struct();
    analysis.within_view_adjacencies = within;
    analysis.copy_mean_by_view_pair = copyMean;
    analysis.crossnode_mean_by_view_pair = crossMean;
    analysis.block_mean_by_view_pair = blockMean;
    analysis.cross_blocks = crossBlocks;
    analysis.copy_vectors = copyVectors;
    analysis.within_mass = sum(withinValues);
    analysis.crossnode_mass = sum(crossValues);
    analysis.copy_mass = sum(copyValues);
    analysis.within_mass_fraction = sum(withinValues)/totalMass;
    analysis.crossnode_mass_fraction = sum(crossValues)/totalMass;
    analysis.copy_mass_fraction = sum(copyValues)/totalMass;
    analysis.within_display_density = withinDensity;
    analysis.crossnode_display_density = crossDensity;
    analysis.copy_display_density = copyDensity;
    analysis.view_graph_log_horizon_correlation = viewCorr;
    analysis.copy_log_horizon_correlation = copyCorr;
    analysis.crossnode_log_horizon_correlation = crossCorr;
    analysis.view_correlation_edge_count = viewCount;
    analysis.copy_correlation_edge_count = copyCount;
    analysis.cross_correlation_edge_count = crossCount;
    analysis.B_maximum_column_cosine = Bcosine;
    analysis.B_mean_normalized_entropy = mean(Bentropy);
    analysis.B_column_entropy = Bentropy;
    analysis.B_column_roughness = Broughness;
    analysis.B_geometric_horizon_centres = Bcentres;
    analysis.strongest_cross_view_pair = strongestCrossPair;
    analysis.strongest_copy_view_pair = strongestCopyPair;
    analysis.strongest_view_graph_pair = strongestViewPair;
    analysis.raw_view_similarity = inputInfo.raw_view_similarity;
    analysis.within_sector_stats = withinStats;
    analysis.cross_sector_stats = crossStats;
    analysis.mode_sector_stats = modeStats;
    analysis.stock_strength_table = stockStrength;
    analysis.sector_strength_table = sectorStrength;
    analysis.mean_within_sector_ratio = mean(withinStats.within_between_ratio,'omitnan');
    analysis.mean_within_sector_modularity = mean(withinStats.sector_modularity,'omitnan');
    analysis.mean_cross_same_sector_ratio = mean(crossStats.same_between_ratio,'omitnan');
end


function compact = finance_compact_fit(fit)
    compact = fit;
    if isfield(compact,'observed'), compact = rmfield(compact,'observed'); end
    if isfield(compact,'all_restarts')
        restartSummary = repmat(struct('restart',NaN,'status','', ...
            'objective',NaN,'stationarity_combined',NaN,'iterations',NaN), ...
            numel(compact.all_restarts),1);
        for q = 1:numel(compact.all_restarts)
            rr = compact.all_restarts{q};
            restartSummary(q).restart = q;
            restartSummary(q).status = rr.status;
            restartSummary(q).objective = rr.objective;
            restartSummary(q).stationarity_combined = ...
                rr.stationarity_combined;
            restartSummary(q).iterations = rr.iterations;
        end
        compact.restart_summary = restartSummary;
        compact = rmfield(compact,'all_restarts');
    end
end


function T = finance_rows_to_table(rows)
    if isempty(rows)
        T = table();
    elseif isscalar(rows)
        T = struct2table(rows,'AsArray',true);
    else
        T = struct2table(rows);
    end
end


function row = finance_fill_success_row(row,fit,analysis,runtimeSeconds)
    row.status = fit.status;
    row.error_message = '';
    row.runtime_seconds = runtimeSeconds;
    row.objective = fit.objective;
    row.stationarity_residual = fit.stationarity_combined;
    row.iterations = fit.iterations;
    row.best_restart = fit.best_restart;
    row.minimum_intra_degree = fit.minimum_intra_degree;
    row.minimum_inter_degree = fit.minimum_inter_degree;
    row.minimum_copy_degree = fit.minimum_copy_degree;
    row.minimum_supra_degree = fit.minimum_supra_degree;
    row.minimum_view_degree = fit.minimum_view_degree;
    row.B_maximum_column_cosine = analysis.B_maximum_column_cosine;
    row.B_mean_normalized_entropy = analysis.B_mean_normalized_entropy;
    row.view_graph_log_horizon_correlation = ...
        analysis.view_graph_log_horizon_correlation;
    row.copy_log_horizon_correlation = analysis.copy_log_horizon_correlation;
    row.crossnode_log_horizon_correlation = ...
        analysis.crossnode_log_horizon_correlation;
    row.within_mass_fraction = analysis.within_mass_fraction;
    row.crossnode_mass_fraction = analysis.crossnode_mass_fraction;
    row.copy_mass_fraction = analysis.copy_mass_fraction;
    row.within_display_density = analysis.within_display_density;
    row.crossnode_display_density = analysis.crossnode_display_density;
    row.copy_display_density = analysis.copy_display_density;
    row.mean_within_sector_ratio = analysis.mean_within_sector_ratio;
    row.mean_within_sector_modularity = analysis.mean_within_sector_modularity;
    row.mean_cross_same_sector_ratio = analysis.mean_cross_same_sector_ratio;
end


function row = finance_empty_candidate_row()
    row = struct( ...
        'candidate_index',NaN, ...
        'candidate_id','', ...
        'rank',NaN, ...
        'profile_id',NaN, ...
        'profile_name','', ...
        'status','not_run', ...
        'error_message','', ...
        'runtime_seconds',NaN, ...
        'objective',NaN, ...
        'stationarity_residual',NaN, ...
        'iterations',NaN, ...
        'best_restart',NaN, ...
        'alpha_a',NaN,'alpha_c',NaN, ...
        'beta_a',NaN,'beta_c',NaN, ...
        'tau_intra',NaN,'tau_inter',NaN, ...
        'eta_c',NaN,'lambda_B',NaN, ...
        'epsilon_intra',NaN,'epsilon_inter',NaN, ...
        'n_restarts',NaN,'max_iter',NaN, ...
        'minimum_intra_degree',NaN, ...
        'minimum_inter_degree',NaN, ...
        'minimum_copy_degree',NaN, ...
        'minimum_supra_degree',NaN, ...
        'minimum_view_degree',NaN, ...
        'B_maximum_column_cosine',NaN, ...
        'B_mean_normalized_entropy',NaN, ...
        'view_graph_log_horizon_correlation',NaN, ...
        'copy_log_horizon_correlation',NaN, ...
        'crossnode_log_horizon_correlation',NaN, ...
        'within_mass_fraction',NaN, ...
        'crossnode_mass_fraction',NaN, ...
        'copy_mass_fraction',NaN, ...
        'within_display_density',NaN, ...
        'crossnode_display_density',NaN, ...
        'copy_display_density',NaN, ...
        'mean_within_sector_ratio',NaN, ...
        'mean_within_sector_modularity',NaN, ...
        'mean_cross_same_sector_ratio',NaN);
end


function finance_save_candidate_tables(data,fit,analysis,outDir,candidateId,cfg)
    writetable(analysis.within_sector_stats, ...
        fullfile(outDir,'within_horizon_sector_statistics.csv'));
    writetable(analysis.cross_sector_stats, ...
        fullfile(outDir,'cross_horizon_sector_statistics.csv'));
    writetable(analysis.mode_sector_stats, ...
        fullfile(outDir,'mode_sector_statistics.csv'));
    writetable(analysis.stock_strength_table, ...
        fullfile(outDir,'stock_strengths.csv'));
    writetable(analysis.sector_strength_table, ...
        fullfile(outDir,'sector_strengths.csv'));

    crossTable = finance_top_cross_edge_table(fit.A_supra,data,cfg.top_cross_edges);
    copyTable = finance_top_copy_edge_table(fit.A_supra,data,cfg.top_copy_edges);
    withinTable = finance_top_within_edge_table( ...
        analysis.within_view_adjacencies,data,cfg.top_within_edges_per_view);
    writetable(crossTable,fullfile(outDir,'top_cross_stock_cross_horizon_edges.csv'));
    writetable(copyTable,fullfile(outDir,'top_same_stock_copy_edges.csv'));
    writetable(withinTable,fullfile(outDir,'top_within_horizon_edges.csv'));

    modeTable = table((1:size(fit.B,2))', ...
        analysis.B_column_entropy(:),analysis.B_column_roughness(:), ...
        analysis.B_geometric_horizon_centres(:), ...
        'VariableNames',{'mode','normalized_entropy','profile_roughness', ...
        'geometric_horizon_centre_days'});
    writetable(modeTable,fullfile(outDir,'mode_profile_diagnostics.csv'));

    context = struct('candidate_id',candidateId,'analysis',analysis);
    save(fullfile(outDir,'analysis_tables.mat'),'context','-v7.3');
end


function finance_save_candidate_figures(data,Y,fit,analysis,outDir,candidateId,cfg)
    finance_plot_convergence(fit,fullfile(outDir,'01_convergence.png'), ...
        candidateId,cfg);
    finance_plot_profiles_states(data,fit,analysis, ...
        fullfile(outDir,'02_modal_profiles_sector_graphs_copy.png'), ...
        candidateId,cfg);
    finance_plot_within_heatmaps(data,analysis, ...
        fullfile(outDir,'03_within_horizon_adjacencies.png'),candidateId,cfg);
    finance_plot_supra_overview(data,fit,analysis, ...
        fullfile(outDir,'04_supra_overview.png'),candidateId,cfg);
    finance_plot_horizon_diagnostics(data,fit,analysis, ...
        fullfile(outDir,'05_horizon_dependence.png'),candidateId,cfg);
    finance_plot_top_edge_bars(data,fit,analysis,'cross', ...
        fullfile(outDir,'06_top_cross_stock_cross_horizon_edges.png'), ...
        candidateId,cfg);
    finance_plot_top_edge_bars(data,fit,analysis,'copy', ...
        fullfile(outDir,'07_top_same_stock_copy_edges.png'), ...
        candidateId,cfg);
    finance_plot_sector_diagnostics(data,fit,analysis, ...
        fullfile(outDir,'08_sector_diagnostics.png'),candidateId,cfg);
    finance_plot_strongest_cross_sector_matrix(data,analysis, ...
        fullfile(outDir,'09_strongest_cross_horizon_sector_matrix.png'), ...
        candidateId,cfg);
    finance_plot_stock_strengths(data,analysis, ...
        fullfile(outDir,'10_stock_and_sector_strengths.png'),candidateId,cfg);

    signalChecksum = struct('sum',sum(Y(:)),'sum_squares',sum(Y(:).^2), ...
        'minimum',min(Y(:)),'maximum',max(Y(:)));
    save(fullfile(outDir,'figure_context.mat'),'signalChecksum','analysis');
end


function finance_plot_convergence(fit,filePath,candidateId,cfg)
    f = finance_new_figure([1300,520]);
    tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
    ax = nexttile(tl,1);
    objective = fit.objective_history(:);
    plot(ax,0:numel(objective)-1,objective,'LineWidth',1.7);
    grid(ax,'on'); box(ax,'on');
    xlabel(ax,'PALM iteration'); ylabel(ax,'objective');
    title(ax,'Objective history');

    ax = nexttile(tl,2);
    R = fit.residual_history;
    if size(R,1)>1
        semilogy(ax,0:size(R,1)-1,max(R(:,1),realmin),'-', ...
            'LineWidth',1.5); hold(ax,'on');
        semilogy(ax,0:size(R,1)-1,max(R(:,2),realmin),'--', ...
            'LineWidth',1.5);
        semilogy(ax,0:size(R,1)-1,max(R(:,3),realmin),'-.', ...
            'LineWidth',1.5);
        legend(ax,{'Theta','Gamma','B'},'Location','best');
    end
    grid(ax,'on'); box(ax,'on');
    xlabel(ax,'PALM iteration'); ylabel(ax,'projected residual');
    title(ax,sprintf('Stationarity: %.3e',fit.stationarity_combined));
    sgtitle(tl,sprintf('%s: convergence diagnostics',candidateId), ...
        'Interpreter','none');
    finance_export_figure(f,filePath,cfg);
end


function finance_plot_profiles_states(data,fit,analysis,filePath,candidateId,cfg)
    r = size(fit.B,2);
    labels = finance_view_labels(data);
    f = finance_new_figure([470*r,1050]);
    tl = tiledlayout(f,3,r,'TileSpacing','compact','Padding','compact');
    for m = 1:r
        ax = nexttile(tl,m);
        bar(ax,1:data.num_views,fit.B(:,m));
        grid(ax,'on'); box(ax,'on');
        set(ax,'XTick',1:data.num_views,'XTickLabel',labels, ...
            'XTickLabelRotation',35);
        ylabel(ax,'B weight');
        title(ax,sprintf('Mode %d: centre %.1f days',m, ...
            analysis.B_geometric_horizon_centres(m)));

        ax = nexttile(tl,r+m);
        A = fit.A_state(data.sector_order,data.sector_order,m);
        imagesc(ax,A); axis(ax,'image'); colorbar(ax);
        finance_draw_sector_boundaries(ax,data.sector_boundaries,'both');
        xlabel(ax,'stocks sorted by sector'); ylabel(ax,'stocks sorted by sector');
        title(ax,sprintf('Mode %d stock-state graph',m));

        ax = nexttile(tl,2*r+m);
        [vals,idx] = maxk(fit.Gamma(:,m),min(15,data.num_stocks));
        barh(ax,vals(end:-1:1));
        set(ax,'YTick',1:numel(idx),'YTickLabel', ...
            cellstr(data.tickers(idx(end:-1:1))));
        grid(ax,'on'); box(ax,'on');
        xlabel(ax,'copy coefficient');
        title(ax,sprintf('Mode %d top stock-copy coefficients',m));
    end
    sgtitle(tl,sprintf('%s: modal horizon profiles and finance structure', ...
        candidateId),'Interpreter','none');
    finance_export_figure(f,filePath,cfg);
end


function finance_plot_within_heatmaps(data,analysis,filePath,candidateId,cfg)
    K = data.num_views;
    labels = finance_view_labels(data);
    nCols = ceil(sqrt(K));
    nRows = ceil(K/nCols);
    f = finance_new_figure([540*nCols,480*nRows]);
    tl = tiledlayout(f,nRows,nCols,'TileSpacing','compact','Padding','compact');
    commonMax = max(analysis.within_view_adjacencies(:));
    commonMax = max(commonMax,eps);
    for k = 1:K
        ax = nexttile(tl,k);
        A = analysis.within_view_adjacencies( ...
            data.sector_order,data.sector_order,k);
        imagesc(ax,A,[0 commonMax]); axis(ax,'image');
        finance_draw_sector_boundaries(ax,data.sector_boundaries,'both');
        title(ax,labels{k});
        xlabel(ax,'stocks sorted by GICS sector');
        ylabel(ax,'stocks sorted by GICS sector');
        if k==K, colorbar(ax); end
    end
    sgtitle(tl,sprintf('%s: learned within-horizon stock graphs',candidateId), ...
        'Interpreter','none');
    finance_export_figure(f,filePath,cfg);
end


function finance_plot_supra_overview(data,fit,analysis,filePath,candidateId,cfg)
    n=data.num_stocks; K=data.num_views;
    labels = finance_view_labels(data);
    f=finance_new_figure([1700,570]);
    tl=tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');

    ax=nexttile(tl,1);
    imagesc(ax,fit.A_supra); axis(ax,'image'); colorbar(ax);
    hold(ax,'on');
    for k=1:K-1
        xline(ax,k*n+0.5,'w:','LineWidth',0.6);
        yline(ax,k*n+0.5,'w:','LineWidth',0.6);
    end
    xlabel(ax,'stock--horizon replica'); ylabel(ax,'stock--horizon replica');
    title(ax,'Full learned supra-adjacency');

    ax=nexttile(tl,2);
    imagesc(ax,analysis.block_mean_by_view_pair);
    axis(ax,'image'); set(ax,'YDir','normal'); colorbar(ax);
    finance_set_view_ticks(ax,labels);
    title(ax,'Mean weight of each horizon block');

    ax=nexttile(tl,3);
    masses=[analysis.within_mass_fraction, ...
        analysis.crossnode_mass_fraction,analysis.copy_mass_fraction];
    bar(ax,masses); ylim(ax,[0 1]); grid(ax,'on'); box(ax,'on');
    set(ax,'XTick',1:3,'XTickLabel',{'within','cross-stock','copy'});
    ylabel(ax,'fraction of unique supra-edge mass');
    title(ax,'Relation-type mass composition');

    sgtitle(tl,sprintf('%s: supra-graph overview',candidateId), ...
        'Interpreter','none');
    finance_export_figure(f,filePath,cfg);
end


function finance_plot_horizon_diagnostics(data,fit,analysis,filePath,candidateId,cfg)
    labels = finance_view_labels(data);
    f = finance_new_figure([1600,1050]);
    tl = tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');

    ax = nexttile(tl,1);
    imagesc(ax,fit.view_graph); axis(ax,'image'); set(ax,'YDir','normal');
    colorbar(ax); finance_set_view_ticks(ax,labels);
    title(ax,sprintf('Induced horizon graph, corr(-log separation)=%.3f', ...
        analysis.view_graph_log_horizon_correlation));

    ax = nexttile(tl,2);
    imagesc(ax,analysis.copy_mean_by_view_pair);
    axis(ax,'image'); set(ax,'YDir','normal'); colorbar(ax);
    finance_set_view_ticks(ax,labels);
    title(ax,sprintf('Mean same-stock copy weight, corr=%.3f', ...
        analysis.copy_log_horizon_correlation));

    ax = nexttile(tl,3);
    imagesc(ax,analysis.crossnode_mean_by_view_pair);
    axis(ax,'image'); set(ax,'YDir','normal'); colorbar(ax);
    finance_set_view_ticks(ax,labels);
    title(ax,sprintf('Mean cross-stock/cross-horizon weight, corr=%.3f', ...
        analysis.crossnode_log_horizon_correlation));

    ax = nexttile(tl,4);
    Mraw = finance_normalize_matrix_offdiag(analysis.raw_view_similarity);
    Mview = finance_normalize_matrix_offdiag(fit.view_graph);
    Mcopy = finance_normalize_matrix_offdiag(analysis.copy_mean_by_view_pair);
    Mcross = finance_normalize_matrix_offdiag(analysis.crossnode_mean_by_view_pair);
    pairLabels = finance_pair_labels(data.horizons);
    upper = triu(true(data.num_views),1);
    vRaw = finance_upper_pair_values(Mraw);
    vView = finance_upper_pair_values(Mview);
    vCopy = finance_upper_pair_values(Mcopy);
    vCross = finance_upper_pair_values(Mcross);
    plot(ax,1:numel(vRaw),vRaw,'-o','LineWidth',1.4); hold(ax,'on');
    plot(ax,1:numel(vView),vView,'-s','LineWidth',1.4);
    plot(ax,1:numel(vCopy),vCopy,'-^','LineWidth',1.4);
    plot(ax,1:numel(vCross),vCross,'-d','LineWidth',1.4);
    grid(ax,'on'); box(ax,'on');
    set(ax,'XTick',1:numel(pairLabels),'XTickLabel',pairLabels, ...
        'XTickLabelRotation',45);
    ylabel(ax,'off-diagonal value / maximum');
    legend(ax,{'raw same-stock similarity','view graph', ...
        'copy edges','cross-stock edges'},'Location','best');
    title(ax,'Learned versus raw horizon-pair structure');

    sgtitle(tl,sprintf('%s: return-horizon diagnostics',candidateId), ...
        'Interpreter','none');
    finance_export_figure(f,filePath,cfg);
end


function finance_plot_top_edge_bars(data,fit,analysis,kind,filePath,candidateId,cfg)
    if strcmp(kind,'cross')
        T = finance_top_cross_edge_table(fit.A_supra,data, ...
            min(15,cfg.top_cross_edges));
        heading = 'Top cross-stock / cross-horizon connections';
    elseif strcmp(kind,'copy')
        T = finance_top_copy_edge_table(fit.A_supra,data, ...
            min(15,cfg.top_copy_edges));
        heading = 'Top same-stock copy connections';
    else
        error('Unknown edge kind: %s',kind);
    end
    if isempty(T), return; end

    f = finance_new_figure([1500,max(650,42*height(T)+260)]);
    ax = axes('Parent',f,'Position',[0.36 0.11 0.59 0.78]);
    values = T.weight(end:-1:1);
    barh(ax,values); grid(ax,'on'); box(ax,'on');
    labels = strings(height(T),1);
    for q=1:height(T)
        rr = height(T)-q+1;
        labels(q) = finance_edge_label(T,rr,kind);
    end
    set(ax,'YTick',1:height(T),'YTickLabel',cellstr(labels),'FontSize',8);
    xlabel(ax,'learned edge weight');
    title(ax,heading);
    annotation(f,'textbox',[0.08 0.92 0.84 0.05], ...
        'String',sprintf('%s: %s',candidateId,heading), ...
        'Interpreter','none','HorizontalAlignment','center', ...
        'FontWeight','bold','EdgeColor','none','FontSize',12);
    finance_export_figure(f,filePath,cfg);
end


function finance_plot_sector_diagnostics(data,fit,analysis,filePath,candidateId,cfg)
    W = analysis.within_sector_stats;
    C = analysis.cross_sector_stats;
    M = analysis.mode_sector_stats;
    f = finance_new_figure([1600,950]);
    tl = tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');

    ax = nexttile(tl,1);
    yyaxis(ax,'left');
    plot(ax,W.horizon_days,W.within_between_ratio,'-o','LineWidth',1.7);
    ylabel(ax,'within / between mean weight');
    yyaxis(ax,'right');
    plot(ax,W.horizon_days,W.sector_modularity,'-s','LineWidth',1.7);
    ylabel(ax,'GICS modularity');
    set(ax,'XScale','log'); grid(ax,'on'); box(ax,'on');
    xlabel(ax,'return horizon (days)');
    title(ax,'Within-horizon sector structure');

    ax = nexttile(tl,2);
    if ~isempty(C)
        scatter(ax,C.log_horizon_separation,C.same_between_ratio, ...
            65,C.mean_cross_weight,'filled'); colorbar(ax);
        grid(ax,'on'); box(ax,'on');
        xlabel(ax,'absolute log horizon separation');
        ylabel(ax,'same-sector / between-sector cross weight');
        title(ax,'Cross-horizon sector consistency');
    end

    ax = nexttile(tl,3);
    bar(ax,M.mode,M.within_between_ratio); grid(ax,'on'); box(ax,'on');
    xlabel(ax,'mode'); ylabel(ax,'within / between mean weight');
    title(ax,'Mode-specific sector concentration');

    ax = nexttile(tl,4);
    bar(ax,1:height(analysis.sector_strength_table), ...
        analysis.sector_strength_table.mean_total_strength);
    set(ax,'XTick',1:height(analysis.sector_strength_table), ...
        'XTickLabel',cellstr(analysis.sector_strength_table.sector), ...
        'XTickLabelRotation',45);
    grid(ax,'on'); box(ax,'on');
    ylabel(ax,'mean supra strength per stock');
    title(ax,'Sector-level learned connectivity');

    sgtitle(tl,sprintf('%s: sector-based post-hoc diagnostics',candidateId), ...
        'Interpreter','none');
    finance_export_figure(f,filePath,cfg);
end


function finance_plot_strongest_cross_sector_matrix(data,analysis,filePath,candidateId,cfg)
    pair = analysis.strongest_cross_view_pair;
    if any(~isfinite(pair)), return; end
    k=pair(1); ell=pair(2);
    C = analysis.cross_blocks{k,ell};
    if isempty(C), return; end
    M = finance_sector_pair_matrix(C,data.sector_index,numel(data.sector_names),false);

    f = finance_new_figure([1450,650]);
    tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
    ax = nexttile(tl,1);
    imagesc(ax,M); axis(ax,'image'); set(ax,'YDir','normal'); colorbar(ax);
    set(ax,'XTick',1:numel(data.sector_names),'YTick',1:numel(data.sector_names), ...
        'XTickLabel',cellstr(data.sector_names),'YTickLabel',cellstr(data.sector_names), ...
        'XTickLabelRotation',45);
    xlabel(ax,sprintf('%d-day destination sector',data.horizons(ell)));
    ylabel(ax,sprintf('%d-day source sector',data.horizons(k)));
    title(ax,'Mean cross-stock edge weight by sector pair');

    ax = nexttile(tl,2);
    topT = finance_top_cross_edges_for_pair_table(C,data,k,ell,15);
    vals = topT.weight(end:-1:1);
    barh(ax,vals); grid(ax,'on'); box(ax,'on');
    labs = strings(height(topT),1);
    for q=1:height(topT)
        rr=height(topT)-q+1;
        labs(q)=sprintf('%s [%s] -- %s [%s]', ...
            char(topT.source_ticker(rr)),char(topT.source_sector(rr)), ...
            char(topT.target_ticker(rr)),char(topT.target_sector(rr)));
    end
    set(ax,'YTick',1:height(topT),'YTickLabel',cellstr(labs),'FontSize',8);
    xlabel(ax,'edge weight');
    title(ax,'Strongest stock pairs for this horizon pair');

    sgtitle(tl,sprintf('%s: strongest mean cross-stock pair %d-day to %d-day', ...
        candidateId,data.horizons(k),data.horizons(ell)), ...
        'Interpreter','none');
    finance_export_figure(f,filePath,cfg);
end


function finance_plot_stock_strengths(data,analysis,filePath,candidateId,cfg)
    T = analysis.stock_strength_table;
    topN = min(cfg.top_stock_strengths,height(T));
    T = sortrows(T,'total_strength','descend');
    top = T(1:topN,:);
    f = finance_new_figure([1550,720]);
    tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');

    ax=nexttile(tl,1);
    b=barh(ax,[top.within_strength,top.cross_stock_strength,top.copy_strength], ...
        'stacked'); %#ok<NASGU>
    set(ax,'YDir','reverse','YTick',1:topN, ...
        'YTickLabel',cellstr(top.ticker));
    grid(ax,'on'); box(ax,'on');
    xlabel(ax,'aggregate learned strength over all horizons');
    legend(ax,{'within','cross-stock','copy'},'Location','best');
    title(ax,'Top stocks by supra-graph strength');

    ax=nexttile(tl,2);
    S=analysis.sector_strength_table;
    bar(ax,[S.mean_within_strength,S.mean_cross_stock_strength,S.mean_copy_strength], ...
        'stacked');
    set(ax,'XTick',1:height(S),'XTickLabel',cellstr(S.sector), ...
        'XTickLabelRotation',45);
    grid(ax,'on'); box(ax,'on');
    ylabel(ax,'mean aggregate strength per stock');
    legend(ax,{'within','cross-stock','copy'},'Location','best');
    title(ax,'Sector-average strength composition');

    sgtitle(tl,sprintf('%s: stock and sector connectivity',candidateId), ...
        'Interpreter','none');
    finance_export_figure(f,filePath,cfg);
end


function T = finance_within_sector_stats(within,data)
    K = data.num_views;
    meanWithin=zeros(K,1); meanBetween=zeros(K,1); ratio=zeros(K,1);
    modularity=zeros(K,1); totalMass=zeros(K,1);
    density=zeros(K,1);
    known = data.sector_known;
    same = (data.sector_labels==data.sector_labels.') & (known & known.');
    between = (data.sector_labels~=data.sector_labels.') & (known & known.');
    upper=triu(true(data.num_stocks),1);
    same=same&upper; between=between&upper;
    for k=1:K
        A=within(:,:,k);
        meanWithin(k)=mean(A(same),'omitnan');
        meanBetween(k)=mean(A(between),'omitnan');
        ratio(k)=meanWithin(k)/max(meanBetween(k),eps);
        modularity(k)=finance_unsigned_modularity(A,data.sector_index);
        totalMass(k)=sum(A(:))/2;
        vals=A(upper);
        density(k)=nnz(vals>0)/numel(vals);
    end
    T=table((1:K)',data.horizons(:),string(finance_view_labels(data))', ...
        meanWithin,meanBetween,ratio,modularity,totalMass,density, ...
        'VariableNames',{'view_index','horizon_days','view_name', ...
        'mean_within_sector_weight','mean_between_sector_weight', ...
        'within_between_ratio','sector_modularity','total_edge_mass','positive_density'});
end


function T = finance_cross_sector_stats(crossBlocks,data)
    K=data.num_views; nPairs=K*(K-1)/2;
    rows=repmat(struct('source_view',NaN,'target_view',NaN, ...
        'source_horizon',NaN,'target_horizon',NaN,'log_horizon_separation',NaN, ...
        'mean_same_sector_weight',NaN,'mean_between_sector_weight',NaN, ...
        'same_between_ratio',NaN,'mean_cross_weight',NaN,'total_cross_mass',NaN), ...
        nPairs,1);
    known=data.sector_known;
    same=(data.sector_labels==data.sector_labels.')&(known&known.');
    between=(data.sector_labels~=data.sector_labels.')&(known&known.');
    same(1:data.num_stocks+1:end)=false;
    between(1:data.num_stocks+1:end)=false;
    q=0;
    for k=1:K
        for ell=k+1:K
            q=q+1; C=crossBlocks{k,ell};
            rows(q).source_view=k; rows(q).target_view=ell;
            rows(q).source_horizon=data.horizons(k);
            rows(q).target_horizon=data.horizons(ell);
            rows(q).log_horizon_separation=abs(log(data.horizons(k)/data.horizons(ell)));
            rows(q).mean_same_sector_weight=mean(C(same),'omitnan');
            rows(q).mean_between_sector_weight=mean(C(between),'omitnan');
            rows(q).same_between_ratio=rows(q).mean_same_sector_weight / ...
                max(rows(q).mean_between_sector_weight,eps);
            off=~eye(data.num_stocks);
            rows(q).mean_cross_weight=mean(C(off),'omitnan');
            rows(q).total_cross_mass=sum(C(:));
        end
    end
    T=struct2table(rows);
end


function T = finance_mode_sector_stats(Astate,data)
    r=size(Astate,3);
    meanWithin=zeros(r,1); meanBetween=zeros(r,1); ratio=zeros(r,1);
    modularity=zeros(r,1); mass=zeros(r,1);
    known=data.sector_known;
    same=(data.sector_labels==data.sector_labels.')&(known&known.');
    between=(data.sector_labels~=data.sector_labels.')&(known&known.');
    upper=triu(true(data.num_stocks),1);
    same=same&upper; between=between&upper;
    for m=1:r
        A=Astate(:,:,m);
        meanWithin(m)=mean(A(same),'omitnan');
        meanBetween(m)=mean(A(between),'omitnan');
        ratio(m)=meanWithin(m)/max(meanBetween(m),eps);
        modularity(m)=finance_unsigned_modularity(A,data.sector_index);
        mass(m)=sum(A(:))/2;
    end
    T=table((1:r)',meanWithin,meanBetween,ratio,modularity,mass, ...
        'VariableNames',{'mode','mean_within_sector_weight', ...
        'mean_between_sector_weight','within_between_ratio', ...
        'sector_modularity','state_edge_mass'});
end


function T = finance_stock_strength_table(A,data)
    n=data.num_stocks; K=data.num_views;
    within=zeros(n,1); cross=zeros(n,1); copy=zeros(n,1);
    for k=1:K
        Ik=(k-1)*n+(1:n);
        within=within+sum(A(Ik,Ik),2);
        for ell=k+1:K
            Il=(ell-1)*n+(1:n);
            C=A(Ik,Il);
            copy=copy+diag(C)+diag(C');
            Ccross=C; Ccross(1:n+1:end)=0;
            cross=cross+sum(Ccross,2)+sum(Ccross,1)';
        end
    end
    total=within+cross+copy;
    T=table((1:n)',data.tickers(:),data.sector_labels(:), ...
        within,cross,copy,total, ...
        'VariableNames',{'node_index','ticker','sector','within_strength', ...
        'cross_stock_strength','copy_strength','total_strength'});
end


function T = finance_sector_strength_table(stockT,data)
    sectors=data.sector_names(:); S=numel(sectors);
    count=zeros(S,1); wi=zeros(S,1); cr=zeros(S,1); cp=zeros(S,1); tot=zeros(S,1);
    for s=1:S
        mask=stockT.sector==sectors(s);
        count(s)=nnz(mask);
        wi(s)=mean(stockT.within_strength(mask),'omitnan');
        cr(s)=mean(stockT.cross_stock_strength(mask),'omitnan');
        cp(s)=mean(stockT.copy_strength(mask),'omitnan');
        tot(s)=mean(stockT.total_strength(mask),'omitnan');
    end
    T=table(sectors,count,wi,cr,cp,tot, ...
        'VariableNames',{'sector','stock_count','mean_within_strength', ...
        'mean_cross_stock_strength','mean_copy_strength','mean_total_strength'});
end


function T = finance_top_cross_edge_table(A,data,maxEdges)
    n=data.num_stocks; K=data.num_views;
    records=repmat(struct('source_stock_index',NaN,'target_stock_index',NaN, ...
        'source_view',NaN,'target_view',NaN,'weight',NaN),0,1);
    for k=1:K
        Ik=(k-1)*n+(1:n);
        for ell=k+1:K
            Il=(ell-1)*n+(1:n); C=A(Ik,Il);
            for i=1:n
                for j=1:n
                    if i~=j && C(i,j)>0
                        records(end+1,1)=struct('source_stock_index',i, ...
                            'target_stock_index',j,'source_view',k, ...
                            'target_view',ell,'weight',C(i,j)); %#ok<AGROW>
                    end
                end
            end
        end
    end
    records=finance_diverse_cross_records(records,maxEdges);
    T=finance_cross_records_to_table(records,data);
end


function records = finance_diverse_cross_records(records,maxEdges)
    if isempty(records), return; end
    [~,ord]=sort([records.weight],'descend'); records=records(ord);
    selected=false(size(records)); pairCount=containers.Map('KeyType','char','ValueType','double');
    stockPairCount=containers.Map('KeyType','char','ValueType','double');
    nSel=0;
    for q=1:numel(records)
        r=records(q);
        vk=sprintf('%d_%d',r.source_view,r.target_view);
        sk=sprintf('%d_%d',min(r.source_stock_index,r.target_stock_index), ...
            max(r.source_stock_index,r.target_stock_index));
        if ~isKey(pairCount,vk), pairCount(vk)=0; end
        if ~isKey(stockPairCount,sk), stockPairCount(sk)=0; end
        if pairCount(vk)>=3 || stockPairCount(sk)>=2, continue; end
        selected(q)=true; nSel=nSel+1;
        pairCount(vk)=pairCount(vk)+1;
        stockPairCount(sk)=stockPairCount(sk)+1;
        if nSel>=maxEdges, break; end
    end
    records=records(selected);
end


function T = finance_cross_records_to_table(records,data)
    if isempty(records)
        T=table(); return;
    end
    N=numel(records);
    rank=(1:N)'; si=zeros(N,1); tj=zeros(N,1); sv=zeros(N,1); tv=zeros(N,1); w=zeros(N,1);
    for q=1:N
        si(q)=records(q).source_stock_index; tj(q)=records(q).target_stock_index;
        sv(q)=records(q).source_view; tv(q)=records(q).target_view; w(q)=records(q).weight;
    end
    T=table(rank,si,data.tickers(si),data.sector_labels(si), ...
        sv,reshape(data.horizons(sv),[],1),tj,data.tickers(tj),data.sector_labels(tj), ...
        tv,reshape(data.horizons(tv),[],1),w, ...
        'VariableNames',{'rank','source_stock_index','source_ticker','source_sector', ...
        'source_view','source_horizon_days','target_stock_index','target_ticker', ...
        'target_sector','target_view','target_horizon_days','weight'});
end


function T = finance_top_copy_edge_table(A,data,maxEdges)
    n=data.num_stocks; K=data.num_views;
    rows=repmat(struct('stock_index',NaN,'source_view',NaN, ...
        'target_view',NaN,'weight',NaN),0,1);
    for k=1:K
        Ik=(k-1)*n+(1:n);
        for ell=k+1:K
            Il=(ell-1)*n+(1:n); C=A(Ik,Il);
            for i=1:n
                if C(i,i)>0
                    rows(end+1,1)=struct('stock_index',i,'source_view',k, ...
                        'target_view',ell,'weight',C(i,i)); %#ok<AGROW>
                end
            end
        end
    end
    if isempty(rows), T=table(); return; end
    [~,ord]=sort([rows.weight],'descend'); rows=rows(ord(1:min(maxEdges,numel(ord))));
    N=numel(rows); rank=(1:N)'; idx=zeros(N,1); sv=zeros(N,1); tv=zeros(N,1); w=zeros(N,1);
    for q=1:N
        idx(q)=rows(q).stock_index; sv(q)=rows(q).source_view;
        tv(q)=rows(q).target_view; w(q)=rows(q).weight;
    end
    T=table(rank,idx,data.tickers(idx),data.sector_labels(idx), ...
        sv,reshape(data.horizons(sv),[],1),tv,reshape(data.horizons(tv),[],1),w, ...
        'VariableNames',{'rank','stock_index','ticker','sector','source_view', ...
        'source_horizon_days','target_view','target_horizon_days','weight'});
end


function T = finance_top_within_edge_table(within,data,maxPerView)
    K=data.num_views; n=data.num_stocks;
    allT=cell(K,1);
    upper=triu(true(n),1);
    [ii,jj]=find(upper);
    for k=1:K
        A=within(:,:,k); vals=A(upper);
        [sorted,ord]=sort(vals,'descend');
        take=ord(1:min(maxPerView,numel(ord)));
        take=take(sorted(1:min(maxPerView,numel(ord)))>0);
        N=numel(take);
        allT{k}=table(repmat(k,N,1),repmat(data.horizons(k),N,1), ...
            ii(take),data.tickers(ii(take)),data.sector_labels(ii(take)), ...
            jj(take),data.tickers(jj(take)),data.sector_labels(jj(take)), ...
            vals(take), ...
            'VariableNames',{'view','horizon_days','stock_i','ticker_i', ...
            'sector_i','stock_j','ticker_j','sector_j','weight'});
    end
    T=vertcat(allT{:});
end


function T = finance_top_cross_edges_for_pair_table(C,data,k,ell,maxEdges)
    n=data.num_stocks; C(1:n+1:end)=0;
    [vals,lin]=sort(C(:),'descend');
    keep=find(vals>0, min(maxEdges,nnz(vals>0)),'first');
    vals=vals(keep); lin=lin(keep);
    [i,j]=ind2sub([n,n],lin);
    T=table((1:numel(vals))',i,data.tickers(i),data.sector_labels(i), ...
        repmat(k,numel(vals),1),repmat(data.horizons(k),numel(vals),1), ...
        j,data.tickers(j),data.sector_labels(j), ...
        repmat(ell,numel(vals),1),repmat(data.horizons(ell),numel(vals),1),vals, ...
        'VariableNames',{'rank','source_stock_index','source_ticker','source_sector', ...
        'source_view','source_horizon_days','target_stock_index','target_ticker', ...
        'target_sector','target_view','target_horizon_days','weight'});
end


function M = finance_sector_pair_matrix(C,sectorIndex,S,excludeDiagonal)
    M=nan(S,S);
    for a=1:S
        for b=1:S
            sub=C(sectorIndex==a,sectorIndex==b);
            if excludeDiagonal && a==b
                sub=sub(~eye(size(sub)));
            end
            if isempty(sub), M(a,b)=NaN; else, M(a,b)=mean(sub(:),'omitnan'); end
        end
    end
end


function label = finance_edge_label(T,row,kind)
    if strcmp(kind,'cross')
        label=sprintf('%s [%s], %dd -- %s [%s], %dd', ...
            char(T.source_ticker(row)),char(T.source_sector(row)),T.source_horizon_days(row), ...
            char(T.target_ticker(row)),char(T.target_sector(row)),T.target_horizon_days(row));
    else
        label=sprintf('%s [%s]: %dd -- %dd',char(T.ticker(row)),char(T.sector(row)), ...
            T.source_horizon_days(row),T.target_horizon_days(row));
    end
end


function [withinValues,crossValues,copyValues] = ...
        finance_unique_component_values(A,n,K)
    withinValues=[]; crossValues=[]; copyValues=[];
    upper=triu(true(n),1);
    off=~eye(n);
    for k=1:K
        Ik=(k-1)*n+(1:n);
        block=A(Ik,Ik);
        withinValues=[withinValues;block(upper)]; %#ok<AGROW>
        for ell=k+1:K
            Il=(ell-1)*n+(1:n);
            C=A(Ik,Il);
            copyValues=[copyValues;diag(C)]; %#ok<AGROW>
            crossValues=[crossValues;C(off)]; %#ok<AGROW>
        end
    end
end


function density = finance_relative_density(values,relativeThreshold)
    if isempty(values) || max(values)<=0
        density=0; return;
    end
    density=nnz(values>relativeThreshold*max(values))/numel(values);
end


function [value,count] = finance_log_horizon_correlation(M,horizons)
    K=numel(horizons); weights=[]; closeness=[];
    for k=1:K
        for ell=k+1:K
            if isfinite(M(k,ell))
                weights(end+1,1)=M(k,ell); %#ok<AGROW>
                closeness(end+1,1)=-abs(log(horizons(k)/horizons(ell))); %#ok<AGROW>
            end
        end
    end
    count=numel(weights);
    if count<3 || std(weights)<=eps || std(closeness)<=eps
        value=NaN;
    else
        C=corrcoef(weights,closeness); value=C(1,2);
    end
end


function value = finance_max_column_cosine(B)
    r=size(B,2); value=0;
    for m=1:r
        for q=m+1:r
            c=(B(:,m)'*B(:,q))/max(norm(B(:,m))*norm(B(:,q)),eps);
            value=max(value,c);
        end
    end
end


function pair = finance_strongest_upper_pair(M)
    K=size(M,1); mask=triu(true(K),1); vals=M(mask);
    if isempty(vals) || all(~isfinite(vals))
        pair=[NaN NaN]; return;
    end
    vals(~isfinite(vals))=-inf;
    [~,q]=max(vals); [ii,jj]=find(mask); pair=[ii(q),jj(q)];
end


function Q = finance_unsigned_modularity(W,comm)
    comm=double(comm(:)); k=sum(W,2); m2=sum(k);
    if m2<eps, Q=0; return; end
    Q=0;
    for c=unique(comm(:))'
        in=(comm==c); Wcc=sum(sum(W(in,in))); kc=sum(k(in));
        Q=Q+Wcc/m2-(kc/m2)^2;
    end
end


function M = finance_normalize_matrix_offdiag(M)
    d=diag(M); M(1:size(M,1)+1:end)=0;
    mx=max(M(:)); if mx>0, M=M/mx; end
    M(1:size(M,1)+1:end)=d*0;
end


function v = finance_upper_pair_values(M)
    K=size(M,1); v=zeros(K*(K-1)/2,1); q=0;
    for k=1:K
        for ell=k+1:K
            q=q+1; v(q)=M(k,ell);
        end
    end
end


function labels = finance_pair_labels(horizons)
    labels={};
    for k=1:numel(horizons)
        for ell=k+1:numel(horizons)
            labels{end+1}=sprintf('%d-%dd',horizons(k),horizons(ell)); %#ok<AGROW>
        end
    end
end


function labels = finance_view_labels(data)
    labels=cell(1,data.num_views);
    for k=1:data.num_views
        labels{k}=sprintf('%d-day',data.horizons(k));
    end
end


function finance_set_view_ticks(ax,labels)
    K=numel(labels);
    set(ax,'XTick',1:K,'YTick',1:K,'XTickLabel',labels, ...
        'YTickLabel',labels,'XTickLabelRotation',35);
end


function finance_draw_sector_boundaries(ax,boundaries,whichAxis)
    hold(ax,'on');
    for b=boundaries(:)'
        if strcmp(whichAxis,'both') || strcmp(whichAxis,'x')
            xline(ax,b,'k-','LineWidth',0.35);
        end
        if strcmp(whichAxis,'both') || strcmp(whichAxis,'y')
            yline(ax,b,'k-','LineWidth',0.35);
        end
    end
end


function [order,names,boundaries,index] = finance_sector_order(labels)
    labels=string(labels(:)); labels(labels=="")="Unknown";
    [names,~,index]=unique(labels,'stable');
    [~,order]=sort(index,'ascend');
    sortedIndex=index(order); boundaries=[];
    for s=1:numel(names)-1
        last=find(sortedIndex==s,1,'last');
        if ~isempty(last), boundaries(end+1)=last+0.5; end %#ok<AGROW>
    end
end


function finance_plot_candidate_summary(T,outDir,cfg)
    ok=~startsWith(string(T.status),'error:') & string(T.status)~="not_run";
    if ~any(ok), return; end
    S=T(ok,:); x=1:height(S);
    f=finance_new_figure([1750,1150]);
    tl=tiledlayout(f,3,3,'TileSpacing','compact','Padding','compact');

    finance_summary_panel(nexttile(tl,1),x,S.stationarity_residual, ...
        'stationarity residual',true);
    finance_summary_panel(nexttile(tl,2),x,S.B_maximum_column_cosine, ...
        'maximum B-column cosine',false);
    finance_summary_panel(nexttile(tl,3),x,S.B_mean_normalized_entropy, ...
        'mean normalized B entropy',false);

    ax=nexttile(tl,4);
    bar(ax,[S.within_mass_fraction,S.crossnode_mass_fraction,S.copy_mass_fraction], ...
        'stacked'); grid(ax,'on'); box(ax,'on');
    ylabel(ax,'mass fraction'); title(ax,'Supra-edge mass composition');
    legend(ax,{'within','cross-stock','copy'},'Location','best');

    finance_summary_panel(nexttile(tl,5),x,S.mean_within_sector_ratio, ...
        'mean within/between sector ratio',false);
    finance_summary_panel(nexttile(tl,6),x,S.mean_within_sector_modularity, ...
        'mean GICS modularity',false);
    finance_summary_panel(nexttile(tl,7),x,S.mean_cross_same_sector_ratio, ...
        'mean cross-view same/between sector ratio',false);
    finance_summary_panel(nexttile(tl,8),x,S.view_graph_log_horizon_correlation, ...
        'view graph corr(-log horizon separation)',false);
    finance_summary_panel(nexttile(tl,9),x,S.runtime_seconds, ...
        'runtime (seconds)',false);

    sgtitle(tl,'Candidate diagnostic catalogue (not a model-selection ranking)');
    finance_export_figure(f,fullfile(outDir,'candidate_diagnostic_catalogue.png'),cfg);

    labelTable=table((1:height(S))',string(S.candidate_id),S.rank,S.profile_id, ...
        'VariableNames',{'plot_index','candidate_id','rank','profile_id'});
    writetable(labelTable,fullfile(outDir,'candidate_plot_index.csv'));
end


function finance_summary_panel(ax,x,y,titleText,useLog)
    if useLog
        semilogy(ax,x,max(y,realmin),'-o','LineWidth',1.1,'MarkerSize',4);
    else
        plot(ax,x,y,'-o','LineWidth',1.1,'MarkerSize',4);
    end
    grid(ax,'on'); box(ax,'on'); xlabel(ax,'successful candidate index');
    title(ax,titleText);
end


function finance_write_profile_table(profiles,ids,filePath)
    P=profiles(ids);
    T=struct2table(P);
    T.profile_id=ids(:);
    T=movevars(T,'profile_id','Before',1);
    writetable(T,filePath);
end


function finance_write_readme(filePath,dataFile,outputRoot,cfg,inputInfo,T)
    fid=fopen(filePath,'w');
    if fid<0, return; end
    cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
    fprintf(fid,'MOSAIC FINANCE MULTI-HORIZON RESULTS\n');
    fprintf(fid,'=====================================\n\n');
    fprintf(fid,'Input MAT-file: %s\n',dataFile);
    fprintf(fid,'Output folder : %s\n\n',outputRoot);
    fprintf(fid,'Data dimensions: n=%d stocks, K=%d horizons, p=%d aligned coordinates.\n', ...
        inputInfo.n,inputInfo.K,inputInfo.p);
    fprintf(fid,'Horizons: %s trading days.\n\n',mat2str(inputInfo.horizons));
    fprintf(fid,'Screen: ranks %s, profile IDs %s, %d restarts, %d maximum iterations.\n', ...
        mat2str(cfg.rank_candidates),mat2str(cfg.profile_ids), ...
        cfg.n_restarts,cfg.max_iter);
    fprintf(fid,'Total candidate combinations: %d.\n\n',height(T));
    fprintf(fid,['No candidate is automatically selected. Objective values are not ', ...
        'comparable across different regularization profiles.\n']);
    fprintf(fid,['Use convergence/stationarity, nondegeneracy, stability, held-out ', ...
        'signal validation, and sector/stock interpretation after freezing ', ...
        'the grid.\n\n']);
    fprintf(fid,'Per-candidate outputs:\n');
    fprintf(fid,'  01_candidate_fits       compact fit and complete analysis MAT files\n');
    fprintf(fid,'  02_candidate_figures    finance-specific figures for each candidate\n');
    fprintf(fid,'  03_candidate_tables     top edges, sector statistics, stock strengths\n');
    fprintf(fid,'  04_summary              cross-candidate catalogue and profile table\n');
end


function finance_write_error_report(filePath,ME,candidateId)
    fid=fopen(filePath,'w'); if fid<0, return; end
    cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
    fprintf(fid,'Candidate: %s\n',candidateId);
    fprintf(fid,'Identifier: %s\n',ME.identifier);
    fprintf(fid,'Message: %s\n\n',ME.message);
    fprintf(fid,'%s\n',getReport(ME,'extended','hyperlinks','off'));
end


function f = finance_new_figure(pixelSize)
    f=figure('Visible','off','Color','w','Position',[50 50 pixelSize]);
end


function finance_export_figure(f,filePath,cfg)
    folder=fileparts(filePath); finance_ensure_folder(folder);
    try
        exportgraphics(f,filePath,'Resolution',cfg.figure_resolution);
    catch
        print(f,filePath,'-dpng',sprintf('-r%d',cfg.figure_resolution));
    end
    close(f);
end


function finance_ensure_folder(folderPath)
    if ~exist(folderPath,'dir'), mkdir(folderPath); end
end


function out = finance_merge_struct(base,override)
    out=base;
    if isempty(override), return; end
    f=fieldnames(override);
    for q=1:numel(f), out.(f{q})=override.(f{q}); end
end

function result = solve_mosaic_crossview_embedded(observed, opts)
%SOLVE_MOSAIC_CROSSVIEW Three-block PALM solver for no-mass MOSAIC.
%
% Variables
%   Theta >= 0 : E-by-r cross-node modal edge coefficients
%   Gamma >= 0 : n-by-r same-node copy coefficients
%   B           : K-by-r, each column in the probability simplex
%
% The objective is the simplified no-prescribed-mass MOSAIC objective:
% cross-node and copy smoothness, l1/l2 coefficient penalties, two
% complementary smoothed log-degree penalties, and modal-profile
% diversity.  The two penalized degrees are
%   D_intra = (H*Theta)*(B.^2)'
%   D_inter = (H*Theta+Gamma)*(B.*(1-B))'.
% D_supra=D_intra+D_inter and d_view=sum(D_inter,1)' are derived
% diagnostics and are not penalized again.
%
% Important: epsilon-smoothed logarithms are penalties, not strict barriers.
% With beta_a,beta_c>0, positive epsilons, positive line-search margins,
% sufficient-decrease backtracking, and bounded simplex B, PALM theory
% supports convergence of the iterate sequence to a critical point (not a
% global optimum). Initial profiles are random interior-simplex points and
% contain no planted modal-center information.

if nargin < 2, opts = struct(); end
opts = defaults(opts);
[problem, observed] = prepare_problem(observed, opts);
if strcmpi(opts.operation, 'gradient_check')
    result = gradient_check(problem, opts);
    return
end

old_rng = rng;
cleanup = onCleanup(@() rng(old_rng));
rng(opts.seed, 'twister');

records = cell(opts.n_restarts,1);
bestJ = inf;
best = 0;
for s = 1:opts.n_restarts
    [Theta,Gamma,B] = initialize(problem,opts,s);
    records{s} = one_restart(Theta,Gamma,B,problem,opts,s);
    if isfinite(records{s}.objective) && records{s}.objective < bestJ
        bestJ = records{s}.objective;
        best = s;
    end
    if opts.verbose
        fprintf('MOSAIC restart %d/%d: J=% .6e, residual=%.3e, %s\n', ...
            s,opts.n_restarts,records{s}.objective, ...
            records{s}.stationarity_combined,records{s}.status);
    end
end
if best == 0, error('Every MOSAIC restart failed.'); end
result = records{best};
result.best_restart = best;
result.all_restarts = records;
result.options = opts;
result.cost_scale = problem.cost_scale;
result.observed = observed;
end

function rec = one_restart(Theta,Gamma,B,P,o,restart)
Lth=o.L_theta; Lga=o.L_gamma; LB=o.L_B;
B_initial=B;
hist=nan(o.max_iter+1,1);
res_hist=nan(o.max_iter+1,3);
[J,Gth,~,~,~]=objective_gradient(Theta,Gamma,B,P,o);
hist(1)=J;
status='maximum_iterations';

for it=1:o.max_iter
    J0=J;

    % Theta update.
    [Theta_new,Lth,ok]=backtrack_nonnegative(Theta,Gth,Lth, ...
        @(X)objective_only(X,Gamma,B,P,o),J,o,o.sigma_theta);
    if ~ok, status='theta_line_search_failure'; break; end
    Theta=Theta_new;

    % Gamma update at the new Theta.
    [Jmid,~,Gga,~,~]=objective_gradient(Theta,Gamma,B,P,o);
    [Gamma_new,Lga,ok]=backtrack_nonnegative(Gamma,Gga,Lga, ...
        @(X)objective_only(Theta,X,B,P,o),Jmid,o,o.sigma_gamma);
    if ~ok, status='gamma_line_search_failure'; break; end
    Gamma=Gamma_new;

    % B update at the new Theta and Gamma.
    [Jmid,~,~,GB,~]=objective_gradient(Theta,Gamma,B,P,o);
    [B_new,LB,ok]=backtrack_simplex(B,GB,LB, ...
        @(X)objective_only(Theta,Gamma,X,P,o),Jmid,o,o.sigma_B);
    if ~ok, status='B_line_search_failure'; break; end
    B=B_new;

    [J,Gth,Gga,GB,~]=objective_gradient(Theta,Gamma,B,P,o);
    hist(it+1)=J;
    rth=projected_residual_nonnegative(Theta,Gth,Lth);
    rga=projected_residual_nonnegative(Gamma,Gga,Lga);
    rB=projected_residual_simplex(B,GB,LB);
    res_hist(it+1,:)=[rth,rga,rB];
    relJ=abs(J-J0)/max(1,abs(J0));
    if max([rth,rga,rB]) <= o.stationarity_tol && relJ <= o.objective_tol
        status='converged';
        break
    end
    if ~isfinite(J), status='nonfinite_objective'; break; end
end

last=find(isfinite(hist),1,'last');
if isempty(last), last=1; end
hist=hist(1:last);
res_hist=res_hist(1:last,:);
[J,Gth,Gga,GB,aux]=objective_gradient(Theta,Gamma,B,P,o);
rth=projected_residual_nonnegative(Theta,Gth,Lth);
rga=projected_residual_nonnegative(Gamma,Gga,Lga);
rB=projected_residual_simplex(B,GB,LB);
graph=reconstruct(Theta,Gamma,B,P);

rec=graph;
rec.Theta=Theta; rec.Gamma=Gamma; rec.B=B;
rec.initial_B=B_initial;
rec.objective=J; rec.objective_components=aux.components;
rec.status=status; rec.restart=restart; rec.iterations=last-1;
rec.stationarity_theta=rth;
rec.stationarity_gamma=rga;
rec.stationarity_B=rB;
rec.stationarity_combined=sqrt(rth^2+rga^2+rB^2);
rec.objective_history=hist;
rec.residual_history=res_hist;
rec.minimum_intra_degree=min(aux.D_intra(:));
rec.minimum_inter_degree=min(aux.D_inter(:));
rec.minimum_copy_degree=min(aux.D_copy(:));
rec.minimum_supra_degree=min(aux.D_supra(:));
rec.minimum_view_degree=min(aux.d_view(:));
rec.feasibility=struct( ...
    'theta_nonnegativity',max(0,-min(Theta(:))), ...
    'gamma_nonnegativity',max(0,-min(Gamma(:))), ...
    'B_nonnegativity',max(0,-min(B(:))), ...
    'B_simplex',max(abs(sum(B,1)-1)));
end

function [Xnew,L,ok]=backtrack_nonnegative(X,G,L,fun,J,o,sigma)
L=max([o.L_min,L/o.bt_relax,sigma*(1+1e-12)]);
ok=false;
for q=1:o.max_backtracks
    Xnew=max(X-G/L,0);
    d=Xnew-X;
    Jnew=fun(Xnew);
    upper=J+sum(G(:).*d(:))+0.5*(L-sigma)*sum(d(:).^2);
    if isfinite(Jnew) && Jnew <= upper+o.bt_tol
        ok=true; return
    end
    L=L*o.bt_factor;
end
Xnew=X;
end

function [Bnew,L,ok]=backtrack_simplex(B,G,L,fun,J,o,sigma)
L=max([o.L_min,L/o.bt_relax,sigma*(1+1e-12)]);
ok=false;
for q=1:o.max_backtracks
    Bnew=project_simplex_columns(B-G/L);
    d=Bnew-B;
    Jnew=fun(Bnew);
    upper=J+sum(G(:).*d(:))+0.5*(L-sigma)*sum(d(:).^2);
    if isfinite(Jnew) && Jnew <= upper+o.bt_tol
        ok=true; return
    end
    L=L*o.bt_factor;
end
Bnew=B;
end

function [J,Gth,Gga,GB,A]=objective_gradient(Theta,Gamma,B,P,o)
Q=P.H*Theta;
U=B.^2;
W=B.*(1-B);
sa=sum(Theta,1)';
sc=sum(Gamma,1)';
h=2*sa+sc;
D_intra=Q*U';
D_cross=Q*W';
D_copy=Gamma*W';
D_inter=D_cross+D_copy;
D_supra=D_intra+D_inter;
d_view=W*h;

Da=zeros(P.E,P.r); Dc=zeros(P.n,P.r);
Ta=zeros(P.K,P.K,P.r); Tc=zeros(P.K,P.K,P.r);
for m=1:P.r
    b=B(:,m);
    for e=1:P.E
        Se=P.S_edge(:,:,e);
        Da(e,m)=b'*Se*b;
        Ta(:,:,m)=Ta(:,:,m)+Theta(e,m)*Se;
    end
    for i=1:P.n
        Ri=P.R_copy(:,:,i);
        Dc(i,m)=0.5*o.eta_c*(b'*Ri*b);
        Tc(:,:,m)=Tc(:,:,m)+Gamma(i,m)*Ri;
    end
end

div=0;
for m=1:P.r
    for q=m+1:P.r, div=div+(B(:,m)'*B(:,q))^2; end
end
c_data=sum(Theta(:).*Da(:))+sum(Gamma(:).*Dc(:));
c_l1=o.alpha_a*sum(Theta(:))+o.alpha_c*sum(Gamma(:));
c_l2=0.5*o.beta_a*sum(Theta(:).^2)+0.5*o.beta_c*sum(Gamma(:).^2);
c_intra=-o.tau_intra*sum(log(D_intra(:)+o.epsilon_intra));
c_inter=-o.tau_inter*sum(log(D_inter(:)+o.epsilon_inter));
c_B=o.lambda_B*div;
J=c_data+c_l1+c_l2+c_intra+c_inter+c_B;

Uintra=1./(D_intra+o.epsilon_intra);
Uinter=1./(D_inter+o.epsilon_inter);
Gth=Da+o.alpha_a+o.beta_a*Theta ...
    -o.tau_intra*P.H'*(Uintra*U) ...
    -o.tau_inter*P.H'*(Uinter*W);
Gga=Dc+o.alpha_c+o.beta_c*Gamma ...
    -o.tau_inter*Uinter*W;

GB=zeros(P.K,P.r);
for m=1:P.r
    b=B(:,m);
    q_intra=sum((Q(:,m).*ones(1,P.K)).*Uintra,1)';
    q_inter=sum(((Q(:,m)+Gamma(:,m)).*ones(1,P.K)).*Uinter,1)';
    GB(:,m)=2*Ta(:,:,m)*b+o.eta_c*Tc(:,:,m)*b ...
        -2*o.tau_intra*b.*q_intra ...
        -o.tau_inter*(1-2*b).*q_inter;
    for q=1:P.r
        if q~=m
            GB(:,m)=GB(:,m)+2*o.lambda_B*(b'*B(:,q))*B(:,q);
        end
    end
end

A=struct();
A.D_intra=D_intra; A.D_cross=D_cross;
A.D_copy=D_copy; A.D_inter=D_inter;
A.D_supra=D_supra; A.d_view=d_view;
A.components=struct('data',c_data,'l1',c_l1,'l2',c_l2, ...
    'log_intra',c_intra,'log_inter',c_inter,'diversity',c_B);
end

function J=objective_only(Theta,Gamma,B,P,o)
[J,~,~,~,~]=objective_gradient(Theta,Gamma,B,P,o);
end

function [P,observed]=prepare_problem(observed,o)
required={'n','K','r'};
for q=1:numel(required)
    if ~isfield(observed,required{q}), error('observed.%s is required.',required{q}); end
end
n=observed.n; K=observed.K; r=observed.r;
[ii,jj]=find(triu(true(n),1)); E=numel(ii);
H=sparse([ii;jj],[(1:E)';(1:E)'],1,n,E);

if isfield(observed,'S_edge') && isfield(observed,'R_copy')
    Sedge=observed.S_edge; Rcopy=observed.R_copy;
elseif isfield(observed,'S') && isequal(size(observed.S),[n*K,n*K])
    Sedge=zeros(K,K,E); Rcopy=zeros(K,K,n);
    Sfull=0.5*(observed.S+observed.S');
    for e=1:E
        for k=1:K
            for ell=1:K
                Sedge(k,ell,e)=Sfull((k-1)*n+ii(e),(ell-1)*n+jj(e));
            end
        end
        Sedge(:,:,e)=0.5*(Sedge(:,:,e)+Sedge(:,:,e)');
    end
    for i=1:n
        for k=1:K
            for ell=1:K
                if k~=ell
                    Rcopy(k,ell,i)=Sfull((k-1)*n+i,(ell-1)*n+i);
                end
            end
        end
        Ri=0.5*(Rcopy(:,:,i)+Rcopy(:,:,i)');
        Ri(1:K+1:end)=0;
        Rcopy(:,:,i)=Ri;
    end
else
    if ~isfield(observed,'Y'), error('Provide observed.S or observed.Y.'); end
    Y=observed.Y;
    if size(Y,1)~=n*K
        error('observed.Y must have n*K rows.');
    end
    if isempty(Y)||size(Y,2)<1||any(~isfinite(Y(:)))
        error('observed.Y must contain at least one finite signal sample.');
    end
    d=sum(Y.^2,2)/size(Y,2);
    Sfull=max(d+d'-2*(Y*Y'/size(Y,2)),0);
    observed.S=Sfull;
    [P,observed]=prepare_problem(observed,o);
    return
end
if ~isequal(size(Sedge),[K,K,E]), error('S_edge must be K-by-K-by-E.'); end
if ~isequal(size(Rcopy),[K,K,n]), error('R_copy must be K-by-K-by-n.'); end
if any(~isfinite(Sedge(:)))||any(~isfinite(Rcopy(:)))
    error('S_edge and R_copy must be finite.');
end
for e=1:E
    Sedge(:,:,e)=max(0,0.5*(Sedge(:,:,e)+Sedge(:,:,e)'));
end
for i=1:n
    Ri=max(0,0.5*(Rcopy(:,:,i)+Rcopy(:,:,i)'));
    Ri(1:K+1:end)=0;
    Rcopy(:,:,i)=Ri;
end

if isfield(observed,'cost_scale') && ~isempty(observed.cost_scale)
    scale=observed.cost_scale;
elseif ~isempty(o.cost_scale)
    scale=o.cost_scale;
else
    vals=[Sedge(:);Rcopy(Rcopy>0)];
    vals=vals(isfinite(vals)&vals>0);
    scale=median(vals);
end
if ~isfinite(scale)||scale<=0, error('The cost scale must be positive.'); end
Sedge=Sedge/scale; Rcopy=Rcopy/scale;
P=struct('n',n,'K',K,'r',r,'E',E,'edge_i',ii,'edge_j',jj, ...
    'H',H,'S_edge',Sedge,'R_copy',Rcopy,'cost_scale',scale);
end

function [Theta,Gamma,B]=initialize(P,o,s)
if s==1 && ~isempty(o.initial_Theta), Theta=max(o.initial_Theta,0);
else, Theta=o.theta_scale*(0.25+rand(P.E,P.r)); end
if s==1 && ~isempty(o.initial_Gamma), Gamma=max(o.initial_Gamma,0);
else, Gamma=o.gamma_scale*(0.25+rand(P.n,P.r)); end
if s==1 && ~isempty(o.initial_B), B=project_simplex_columns(o.initial_B);
else
    B=random_interior_simplex(P.K,P.r,o.initial_uniform_mix);
end
end

function B=random_interior_simplex(K,r,uniform_mix)
% Independent Dirichlet(1) columns, mildly mixed with the barycenter.
% This is truth-independent and does not constrain later B iterates.
Z=-log(max(rand(K,r),realmin('double')));
B=Z./sum(Z,1);
B=(1-uniform_mix)*B+uniform_mix/K;
B=B./sum(B,1);
end

function G=reconstruct(Theta,Gamma,B,P)
Astate=zeros(P.n,P.n,P.r);
for m=1:P.r
    A=zeros(P.n);
    A(sub2ind([P.n,P.n],P.edge_i,P.edge_j))=Theta(:,m);
    A=A+A';
    Astate(:,:,m)=A;
end
Asupra=zeros(P.n*P.K);
for m=1:P.r
    bb=B(:,m)*B(:,m)';
    M=bb-diag(diag(bb));
    Asupra=Asupra+kron(bb,Astate(:,:,m))+kron(M,diag(Gamma(:,m)));
end
Asupra=0.5*(Asupra+Asupra');
Asupra(1:size(Asupra,1)+1:end)=0;
sa=sum(Theta,1)'; sc=sum(Gamma,1)';
view=B*diag(2*sa+sc)*B';
view=view-diag(diag(view));
G=struct('A_state',Astate,'A_supra',Asupra,'supra',Asupra, ...
    'view_graph',view,'state_cross_mass',sa,'state_copy_mass',sc);
end

function audit=gradient_check(P,o)
rng(o.seed+991,'twister');
Theta=.2+rand(P.E,P.r); Gamma=.2+rand(P.n,P.r);
B=.2+rand(P.K,P.r); B=B./sum(B,1);
[~,Gt,Gg,Gb]=objective_gradient(Theta,Gamma,B,P,o);
audit=struct();
[audit.theta_fd,audit.theta_analytic,audit.theta_relative_error]= ...
    directional_check(Theta,Gt,@(X)objective_only(X,Gamma,B,P,o),false);
[audit.gamma_fd,audit.gamma_analytic,audit.gamma_relative_error]= ...
    directional_check(Gamma,Gg,@(X)objective_only(Theta,X,B,P,o),false);
[audit.B_fd,audit.B_analytic,audit.B_relative_error]= ...
    directional_check(B,Gb,@(X)objective_only(Theta,Gamma,X,P,o),true);
audit.passed=max([audit.theta_relative_error,audit.gamma_relative_error, ...
    audit.B_relative_error])<o.gradient_check_tol;
end

function [fd,an,err]=directional_check(X,G,fun,tangent)
D=randn(size(X));
if tangent, D=D-mean(D,1); end
D=D/norm(D,'fro'); h=1e-6;
fd=(fun(X+h*D)-fun(X-h*D))/(2*h);
an=sum(G(:).*D(:));
err=abs(fd-an)/max([1,abs(fd),abs(an)]);
end

function r=projected_residual_nonnegative(X,G,L)
r=L*norm(X-max(X-G/L,0),'fro')/max(1,norm(X,'fro'));
end
function r=projected_residual_simplex(B,G,L)
r=L*norm(B-project_simplex_columns(B-G/L),'fro')/max(1,norm(B,'fro'));
end
function B=project_simplex_columns(Y)
B=zeros(size(Y));
for m=1:size(Y,2), B(:,m)=project_simplex(Y(:,m)); end
end
function x=project_simplex(y)
y=y(:); u=sort(y,'descend'); cssv=cumsum(u)-1;
rho=find(u-cssv./(1:numel(y))'>0,1,'last');
if isempty(rho), x=ones(size(y))/numel(y); return; end
x=max(y-cssv(rho)/rho,0);
x=x/sum(x);
end

function o=defaults(o)
d=struct('operation','solve','seed',1,'n_restarts',8,'max_iter',600, ...
 'eta_c',1,'alpha_a',0.01,'alpha_c',0.01,'beta_a',0.10,'beta_c',0.10, ...
 'tau_intra',0.10,'tau_inter',0.10,'lambda_B',0.01, ...
 'epsilon_intra',1e-5,'epsilon_inter',1e-5, ...
 'L_theta',1,'L_gamma',1,'L_B',1,'L_min',1e-10,'bt_factor',2, ...
 'sigma_theta',1e-6,'sigma_gamma',1e-6,'sigma_B',1e-6, ...
 'bt_relax',1.25,'max_backtracks',60,'bt_tol',0, ...
 'stationarity_tol',1e-6,'objective_tol',1e-9, ...
 'theta_scale',0.20,'gamma_scale',0.20,'initial_uniform_mix',0.10, ...
 'initial_Theta',[],'initial_Gamma',[],'initial_B',[], ...
 'cost_scale',[],'gradient_check_tol',5e-5,'verbose',false);
o=merge(d,o);
if o.beta_a<=0||o.beta_c<=0, error('beta_a and beta_c must be positive.'); end
epsnames={'epsilon_intra','epsilon_inter'};
for q=1:numel(epsnames)
    if o.(epsnames{q})<=0, error('%s must be positive.',epsnames{q}); end
end
nonneg={'eta_c','alpha_a','alpha_c','tau_intra','tau_inter','lambda_B'};
for q=1:numel(nonneg)
    if o.(nonneg{q})<0, error('%s must be nonnegative.',nonneg{q}); end
end
positive={'L_min','sigma_theta','sigma_gamma','sigma_B'};
for q=1:numel(positive)
    if o.(positive{q})<=0, error('%s must be positive.',positive{q}); end
end
if o.bt_tol~=0
    error('bt_tol must be zero for the stated sufficient-decrease guarantee.');
end
if o.initial_uniform_mix<0||o.initial_uniform_mix>=1
    error('initial_uniform_mix must lie in [0,1).');
end
end
function a=merge(a,b)
f=fieldnames(b);
for q=1:numel(f), a.(f{q})=b.(f{q}); end
end

%% ========================================================================
%% ORIGINAL CROSS-vs-WITHIN / MODE POSTPROCESSING HELPERS (V3)
%% ========================================================================
function data = pcm_prepare_data(data)
    if ~isfield(data,'X') || ~iscell(data.X) || isempty(data.X)
        error('data.X must be a nonempty cell array.');
    end
    if ~isfield(data,'horizons')
        error('data.horizons is required.');
    end

    data.num_views = numel(data.X);
    data.num_stocks = size(data.X{1},1);
    data.horizons = double(data.horizons(:).');
    if numel(data.horizons) ~= data.num_views
        error('data.horizons must contain one entry per view.');
    end

    if isfield(data,'tickers')
        data.tickers = string(data.tickers(:));
    elseif isfield(data,'node_names')
        data.tickers = string(data.node_names(:));
    else
        data.tickers = "stock_" + string((1:data.num_stocks)');
    end
    if numel(data.tickers) ~= data.num_stocks
        error('Ticker count does not equal the number of stocks.');
    end

    if isfield(data,'sector_labels')
        data.sector_labels = string(data.sector_labels(:));
    else
        data.sector_labels = repmat("Unknown",data.num_stocks,1);
    end
    if numel(data.sector_labels) ~= data.num_stocks
        data.sector_labels = repmat("Unknown",data.num_stocks,1);
    end
    data.sector_labels(strlength(strtrim(data.sector_labels))==0) = "Unknown";

    [~,order] = sortrows(table(data.sector_labels,data.tickers),[1 2]);
    data.sector_order = order;

    sortedSectors = data.sector_labels(order);
    [sectorNames,~,sectorIndexSorted] = unique(sortedSectors,'stable');
    counts = accumarray(sectorIndexSorted,1);
    data.sector_names = sectorNames;
    data.sector_boundaries = cumsum(counts(1:end-1))+0.5;

    data.view_labels = strings(data.num_views,1);
    for k = 1:data.num_views
        data.view_labels(k) = sprintf('%d-day',data.horizons(k));
    end
end


% ========================================================================
% Representative loading and reconstruction
% ========================================================================

function fit = pcm_extract_fit(saved,fitFile)
    if isfield(saved,'compactFit') && isstruct(saved.compactFit)
        fit = saved.compactFit;
    elseif isfield(saved,'fit') && isstruct(saved.fit)
        fit = saved.fit;
    else
        error('No compactFit or fit struct found in:\n  %s',fitFile);
    end
end


function pcm_validate_fit(fit,data,fitFile)
    required = {'A_state','A_supra','B','Gamma'};
    for q = 1:numel(required)
        if ~isfield(fit,required{q})
            error('fit.%s is missing in:\n  %s',required{q},fitFile);
        end
    end

    n = data.num_stocks;
    K = data.num_views;
    if size(fit.A_state,1)~=n || size(fit.A_state,2)~=n
        error('A_state stock dimensions do not match the dataset.');
    end
    if size(fit.B,1)~=K
        error('B has %d rows but the dataset has %d views.',size(fit.B,1),K);
    end
    if ~isequal(size(fit.A_supra),[n*K,n*K])
        error('A_supra must be (nK)-by-(nK).');
    end
    if any(~isfinite(fit.A_state(:))) || any(~isfinite(fit.A_supra(:))) || ...
            any(~isfinite(fit.B(:))) || any(~isfinite(fit.Gamma(:)))
        error('Saved fit contains nonfinite values.');
    end
end


function within = pcm_reconstruct_within(fit,data)
    n = data.num_stocks;
    K = data.num_views;
    r = size(fit.B,2);
    within = zeros(n,n,K);
    for k = 1:K
        for m = 1:r
            within(:,:,k) = within(:,:,k) + ...
                fit.B(k,m)^2 * fit.A_state(:,:,m);
        end
        within(:,:,k) = pcm_symmetrize_nonnegative(within(:,:,k));
        W = within(:,:,k);
        W(1:n+1:end) = 0;
        within(:,:,k) = W;
    end
end


function T = pcm_get_top_cross_relations(saved,fit,data,maxRelations)
    if isfield(saved,'finalRelations') && istable(saved.finalRelations) && ...
            ~isempty(saved.finalRelations)
        T = saved.finalRelations;
        if ismember('rank',T.Properties.VariableNames)
            T = sortrows(T,'rank','ascend');
        elseif ismember('max_weight',T.Properties.VariableNames)
            T = sortrows(T,'max_weight','descend');
        end
        T = T(1:min(maxRelations,height(T)),:);
        return;
    end
    T = pcm_compute_top_cross_relations(fit.A_supra,data,maxRelations);
end


function T = pcm_compute_top_cross_relations(A,data,maxRelations)
    n = data.num_stocks;
    K = data.num_views;
    relationMap = containers.Map('KeyType','char','ValueType','any');

    for k = 1:K
        Ik = (k-1)*n+(1:n);
        for ell = k+1:K
            Il = (ell-1)*n+(1:n);
            C = double(A(Ik,Il));
            C(1:n+1:end) = 0;
            [ii,jj,ww] = find(C);
            for q = 1:numel(ww)
                if ii(q)==jj(q) || ww(q)<=0 || ~isfinite(ww(q))
                    continue;
                end
                stockA = min(ii(q),jj(q));
                stockB = max(ii(q),jj(q));
                hLow = min(data.horizons(k),data.horizons(ell));
                hHigh = max(data.horizons(k),data.horizons(ell));
                key = sprintf('s%04d_s%04d_h%04d_h%04d', ...
                    stockA,stockB,hLow,hHigh);

                if ~isKey(relationMap,key)
                    rec = struct();
                    rec.relation_key = key;
                    rec.stock_a_index = stockA;
                    rec.stock_b_index = stockB;
                    rec.horizon_low_days = hLow;
                    rec.horizon_high_days = hHigh;
                    rec.max_weight = ww(q);
                    rec.best_source_stock_index = ii(q);
                    rec.best_target_stock_index = jj(q);
                    rec.best_source_view = k;
                    rec.best_target_view = ell;
                else
                    rec = relationMap(key);
                    if ww(q)>rec.max_weight
                        rec.max_weight = ww(q);
                        rec.best_source_stock_index = ii(q);
                        rec.best_target_stock_index = jj(q);
                        rec.best_source_view = k;
                        rec.best_target_view = ell;
                    end
                end
                relationMap(key) = rec;
            end
        end
    end

    keysList = relationMap.keys;
    if isempty(keysList)
        T = table();
        return;
    end
    recs = repmat(relationMap(keysList{1}),numel(keysList),1);
    for q = 1:numel(keysList)
        recs(q) = relationMap(keysList{q});
    end
    [~,order] = sort([recs.max_weight],'descend');
    recs = recs(order(1:min(maxRelations,numel(order))));

    N = numel(recs);
    rank = (1:N)';
    relation_key = strings(N,1);
    relation_label = strings(N,1);
    stock_a_index = zeros(N,1);
    stock_b_index = zeros(N,1);
    horizon_low_days = zeros(N,1);
    horizon_high_days = zeros(N,1);
    max_weight = zeros(N,1);
    best_source_stock_index = zeros(N,1);
    best_target_stock_index = zeros(N,1);
    best_source_view = zeros(N,1);
    best_target_view = zeros(N,1);
    for q = 1:N
        rec = recs(q);
        relation_key(q) = string(rec.relation_key);
        stock_a_index(q) = rec.stock_a_index;
        stock_b_index(q) = rec.stock_b_index;
        horizon_low_days(q) = rec.horizon_low_days;
        horizon_high_days(q) = rec.horizon_high_days;
        max_weight(q) = rec.max_weight;
        best_source_stock_index(q) = rec.best_source_stock_index;
        best_target_stock_index(q) = rec.best_target_stock_index;
        best_source_view(q) = rec.best_source_view;
        best_target_view(q) = rec.best_target_view;
        relation_label(q) = sprintf('%s--%s | %dd--%dd', ...
            data.tickers(stock_a_index(q)),data.tickers(stock_b_index(q)), ...
            horizon_low_days(q),horizon_high_days(q));
    end
    normalized_weight_to_top = max_weight/max(max_weight(1),eps);
    T = table(rank,relation_key,relation_label,stock_a_index, ...
        data.tickers(stock_a_index),data.sector_labels(stock_a_index), ...
        stock_b_index,data.tickers(stock_b_index), ...
        data.sector_labels(stock_b_index),horizon_low_days, ...
        horizon_high_days,max_weight,normalized_weight_to_top, ...
        best_source_stock_index,best_target_stock_index,best_source_view, ...
        best_target_view, ...
        'VariableNames',{'rank','relation_key','relation_label', ...
        'stock_a_index','stock_a_ticker','stock_a_sector', ...
        'stock_b_index','stock_b_ticker','stock_b_sector', ...
        'horizon_low_days','horizon_high_days','max_weight', ...
        'normalized_weight_to_top','best_source_stock_index', ...
        'best_target_stock_index','best_source_view','best_target_view'});
end


% ========================================================================
% Cross-vs-within and mode decomposition
% ========================================================================

function [summaryT,modeLong,couplingLong,pairMatrices] = ...
        pcm_analyze_cross_relations(T,fit,within,data,representativeId,rankValue)
    if isempty(T)
        summaryT = table();
        modeLong = table();
        couplingLong = table();
        pairMatrices = {};
        return;
    end

    nRel = height(T);
    K = data.num_views;
    r = size(fit.B,2);
    horizons = data.horizons;

    % Ranking reference populations. Every within-view graph and every
    % cross-horizon block contains the same n-choose-2 unordered stock
    % pairs. The global cross population additionally pools all unordered
    % horizon pairs.
    rankReference = pcm_build_rank_reference(fit.A_supra,within,data);
    nStockPairs = rankReference.n_stock_pairs;
    nGlobalCrossRelations = rankReference.n_global_cross_relations;

    % Core columns.
    representative_id = repmat(string(representativeId),nRel,1);
    model_rank = repmat(rankValue,nRel,1);
    relation_rank = (1:nRel)';
    if ismember('rank',T.Properties.VariableNames)
        relation_rank = T.rank;
    end
    relation_key = strings(nRel,1);
    relation_label = strings(nRel,1);
    stock_a_index = zeros(nRel,1);
    stock_b_index = zeros(nRel,1);
    stock_a_ticker = strings(nRel,1);
    stock_b_ticker = strings(nRel,1);
    stock_a_sector = strings(nRel,1);
    stock_b_sector = strings(nRel,1);
    same_sector = false(nRel,1);
    horizon_low_days = zeros(nRel,1);
    horizon_high_days = zeros(nRel,1);
    cross_weight = zeros(nRel,1);
    reciprocal_cross_weight = zeros(nRel,1);
    reconstructed_cross_weight_from_modes = zeros(nRel,1);
    cross_reconstruction_absolute_error = zeros(nRel,1);
    cross_reconstruction_relative_error = zeros(nRel,1);
    cross_weight_normalized_to_representative_top = zeros(nRel,1);
    within_low_horizon_weight = zeros(nRel,1);
    within_high_horizon_weight = zeros(nRel,1);
    within_endpoint_mean = zeros(nRel,1);
    within_endpoint_geometric_mean = zeros(nRel,1);
    within_endpoint_maximum = zeros(nRel,1);
    strongest_within_horizon_days = zeros(nRel,1);
    strongest_within_weight = zeros(nRel,1);
    cross_to_endpoint_mean_ratio = zeros(nRel,1);
    cross_to_endpoint_geometric_mean_ratio = zeros(nRel,1);
    cross_to_endpoint_max_ratio = zeros(nRel,1);
    cross_to_strongest_within_ratio = zeros(nRel,1);
    cross_minus_endpoint_maximum = zeros(nRel,1);
    cross_is_larger_than_both_endpoint_within = false(nRel,1);
    selected_cross_fraction_of_pair_maximum = zeros(nRel,1);
    dominant_mode = zeros(nRel,1);
    dominant_mode_cross_fraction = zeros(nRel,1);

    % Rank-prominence diagnostics.
    cross_global_rank = zeros(nRel,1);
    cross_global_tie_count = zeros(nRel,1);
    cross_global_total_relations = repmat(nGlobalCrossRelations,nRel,1);
    cross_global_top_fraction_percent = zeros(nRel,1);
    cross_global_prominence_percentile = zeros(nRel,1);

    cross_block_rank = zeros(nRel,1);
    cross_block_tie_count = zeros(nRel,1);
    cross_block_total_stock_pairs = repmat(nStockPairs,nRel,1);
    cross_block_top_fraction_percent = zeros(nRel,1);
    cross_block_prominence_percentile = zeros(nRel,1);

    within_low_rank = zeros(nRel,1);
    within_low_tie_count = zeros(nRel,1);
    within_low_total_stock_pairs = repmat(nStockPairs,nRel,1);
    within_low_top_fraction_percent = zeros(nRel,1);
    within_low_prominence_percentile = zeros(nRel,1);

    within_high_rank = zeros(nRel,1);
    within_high_tie_count = zeros(nRel,1);
    within_high_total_stock_pairs = repmat(nStockPairs,nRel,1);
    within_high_top_fraction_percent = zeros(nRel,1);
    within_high_prominence_percentile = zeros(nRel,1);

    best_endpoint_within_rank = zeros(nRel,1);
    best_endpoint_within_prominence_percentile = zeros(nRel,1);
    cross_block_rank_advantage_over_best_within = zeros(nRel,1);
    cross_block_prominence_advantage_over_best_within = zeros(nRel,1);
    cross_block_rank_ratio_to_best_within = zeros(nRel,1);
    cross_more_prominent_than_both_endpoint_within = false(nRel,1);
    cross_more_prominent_than_at_least_one_endpoint_within = false(nRel,1);

    withinByHorizon = zeros(nRel,K);
    AModeByRelation = zeros(nRel,r);
    crossModeContribution = zeros(nRel,r);
    pairMatrices = cell(nRel,1);

    modeRows = cell(nRel*r,1);
    couplingRows = cell(nRel*K*K,1);
    modeCounter = 0;
    couplingCounter = 0;

    for q = 1:nRel
        [i,j,hLow,hHigh,key,label] = pcm_relation_identity(T,q,data);
        k = pcm_horizon_index(horizons,hLow);
        ell = pcm_horizon_index(horizons,hHigh);

        relation_key(q) = key;
        relation_label(q) = label;
        stock_a_index(q) = i;
        stock_b_index(q) = j;
        stock_a_ticker(q) = data.tickers(i);
        stock_b_ticker(q) = data.tickers(j);
        stock_a_sector(q) = data.sector_labels(i);
        stock_b_sector(q) = data.sector_labels(j);
        same_sector(q) = stock_a_sector(q)==stock_b_sector(q);
        horizon_low_days(q) = hLow;
        horizon_high_days(q) = hHigh;

        % Complete pair-specific K-by-K matrix. Its diagonal equals the
        % within-view edge weights of this same stock pair.
        aMode = reshape(fit.A_state(i,j,:),1,r);
        Ppair = fit.B * diag(aMode) * fit.B';
        Ppair = pcm_symmetrize_nonnegative(Ppair);
        pairMatrices{q} = Ppair;

        orientation1 = Ppair(k,ell);
        orientation2 = Ppair(ell,k);
        reconstructed_cross_weight_from_modes(q) = max(orientation1,orientation2);

        supra1 = fit.A_supra((k-1)*data.num_stocks+i, ...
            (ell-1)*data.num_stocks+j);
        supra2 = fit.A_supra((k-1)*data.num_stocks+j, ...
            (ell-1)*data.num_stocks+i);
        cross_weight(q) = max(supra1,supra2);
        reciprocal_cross_weight(q) = min(supra1,supra2);
        cross_reconstruction_absolute_error(q) = abs(cross_weight(q) - ...
            reconstructed_cross_weight_from_modes(q));
        cross_reconstruction_relative_error(q) = ...
            cross_reconstruction_absolute_error(q) / max(cross_weight(q),eps);

        % Rank this same stock pair in (i) all cross relations, (ii) its
        % selected horizon-pair block, and (iii) each endpoint within graph.
        [cross_global_rank(q),cross_global_tie_count(q)] = ...
            pcm_descending_competition_rank( ...
            cross_weight(q),rankReference.global_cross_weights);
        [cross_block_rank(q),cross_block_tie_count(q)] = ...
            pcm_descending_competition_rank( ...
            cross_weight(q),rankReference.cross_block_weights{k,ell});

        withinLowValue = within(i,j,k);
        withinHighValue = within(i,j,ell);
        [within_low_rank(q),within_low_tie_count(q)] = ...
            pcm_descending_competition_rank( ...
            withinLowValue,rankReference.within_weights{k});
        [within_high_rank(q),within_high_tie_count(q)] = ...
            pcm_descending_competition_rank( ...
            withinHighValue,rankReference.within_weights{ell});

        cross_global_top_fraction_percent(q) = ...
            pcm_top_fraction_percent(cross_global_rank(q), ...
            nGlobalCrossRelations);
        cross_global_prominence_percentile(q) = ...
            pcm_prominence_percentile(cross_global_rank(q), ...
            nGlobalCrossRelations);
        cross_block_top_fraction_percent(q) = ...
            pcm_top_fraction_percent(cross_block_rank(q),nStockPairs);
        cross_block_prominence_percentile(q) = ...
            pcm_prominence_percentile(cross_block_rank(q),nStockPairs);
        within_low_top_fraction_percent(q) = ...
            pcm_top_fraction_percent(within_low_rank(q),nStockPairs);
        within_low_prominence_percentile(q) = ...
            pcm_prominence_percentile(within_low_rank(q),nStockPairs);
        within_high_top_fraction_percent(q) = ...
            pcm_top_fraction_percent(within_high_rank(q),nStockPairs);
        within_high_prominence_percentile(q) = ...
            pcm_prominence_percentile(within_high_rank(q),nStockPairs);

        best_endpoint_within_rank(q) = min( ...
            within_low_rank(q),within_high_rank(q));
        best_endpoint_within_prominence_percentile(q) = max( ...
            within_low_prominence_percentile(q), ...
            within_high_prominence_percentile(q));
        cross_block_rank_advantage_over_best_within(q) = ...
            best_endpoint_within_rank(q)-cross_block_rank(q);
        cross_block_prominence_advantage_over_best_within(q) = ...
            cross_block_prominence_percentile(q)- ...
            best_endpoint_within_prominence_percentile(q);
        cross_block_rank_ratio_to_best_within(q) = ...
            cross_block_rank(q)/max(best_endpoint_within_rank(q),1);
        cross_more_prominent_than_both_endpoint_within(q) = ...
            cross_block_rank(q)<within_low_rank(q) && ...
            cross_block_rank(q)<within_high_rank(q);
        cross_more_prominent_than_at_least_one_endpoint_within(q) = ...
            cross_block_rank(q)<within_low_rank(q) || ...
            cross_block_rank(q)<within_high_rank(q);

        withinVector = diag(Ppair).';
        withinByHorizon(q,:) = withinVector;
        within_low_horizon_weight(q) = withinVector(k);
        within_high_horizon_weight(q) = withinVector(ell);
        within_endpoint_mean(q) = mean([withinVector(k),withinVector(ell)]);
        within_endpoint_geometric_mean(q) = sqrt(max(withinVector(k),0) * ...
            max(withinVector(ell),0));
        within_endpoint_maximum(q) = max(withinVector(k),withinVector(ell));
        [strongest_within_weight(q),bestWithin] = max(withinVector);
        strongest_within_horizon_days(q) = horizons(bestWithin);

        cross_to_endpoint_mean_ratio(q) = cross_weight(q) / ...
            max(within_endpoint_mean(q),eps);
        cross_to_endpoint_geometric_mean_ratio(q) = cross_weight(q) / ...
            max(within_endpoint_geometric_mean(q),eps);
        cross_to_endpoint_max_ratio(q) = cross_weight(q) / ...
            max(within_endpoint_maximum(q),eps);
        cross_to_strongest_within_ratio(q) = cross_weight(q) / ...
            max(strongest_within_weight(q),eps);
        cross_minus_endpoint_maximum(q) = cross_weight(q) - ...
            within_endpoint_maximum(q);
        cross_is_larger_than_both_endpoint_within(q) = ...
            cross_weight(q) > within_endpoint_maximum(q);
        selected_cross_fraction_of_pair_maximum(q) = cross_weight(q) / ...
            max(max(Ppair(:)),eps);

        AModeByRelation(q,:) = aMode;
        contributions = aMode .* (fit.B(k,:).*fit.B(ell,:));
        crossModeContribution(q,:) = contributions;
        [dominantValue,dominant_mode(q)] = max(contributions);
        dominant_mode_cross_fraction(q) = dominantValue / ...
            max(sum(contributions),eps);

        for m = 1:r
            modeCounter = modeCounter+1;
            rec = struct();
            rec.representative_id = string(representativeId);
            rec.model_rank = rankValue;
            rec.relation_rank = relation_rank(q);
            rec.relation_key = key;
            rec.relation_label = label;
            rec.stock_a_ticker = data.tickers(i);
            rec.stock_b_ticker = data.tickers(j);
            rec.horizon_low_days = hLow;
            rec.horizon_high_days = hHigh;
            rec.mode = m;
            rec.A_mode_stock_pair_weight = aMode(m);
            rec.B_low_horizon_weight = fit.B(k,m);
            rec.B_high_horizon_weight = fit.B(ell,m);
            rec.cross_horizon_mode_contribution = contributions(m);
            rec.cross_horizon_mode_fraction = contributions(m) / ...
                max(sum(contributions),eps);
            rec.within_low_mode_contribution = aMode(m)*fit.B(k,m)^2;
            rec.within_high_mode_contribution = aMode(m)*fit.B(ell,m)^2;
            modeRows{modeCounter} = rec;
        end

        for kk = 1:K
            for ll = 1:K
                couplingCounter = couplingCounter+1;
                rec = struct();
                rec.representative_id = string(representativeId);
                rec.model_rank = rankValue;
                rec.relation_rank = relation_rank(q);
                rec.relation_key = key;
                rec.relation_label = label;
                rec.stock_a_ticker = data.tickers(i);
                rec.stock_b_ticker = data.tickers(j);
                rec.row_horizon_days = horizons(kk);
                rec.column_horizon_days = horizons(ll);
                rec.weight = Ppair(kk,ll);
                rec.normalized_by_pair_maximum = Ppair(kk,ll) / ...
                    max(max(Ppair(:)),eps);
                rec.is_within_view_diagonal = kk==ll;
                rec.is_selected_cross_horizon_pair = ...
                    (kk==k && ll==ell) || (kk==ell && ll==k);
                couplingRows{couplingCounter} = rec;
            end
        end
    end

    topCross = max(cross_weight(1),eps);
    if ismember('normalized_weight_to_top',T.Properties.VariableNames)
        cross_weight_normalized_to_representative_top = ...
            T.normalized_weight_to_top;
    else
        cross_weight_normalized_to_representative_top = cross_weight/topCross;
    end

    summaryT = table(representative_id,model_rank,relation_rank, ...
        relation_key,relation_label,stock_a_index,stock_a_ticker, ...
        stock_a_sector,stock_b_index,stock_b_ticker,stock_b_sector, ...
        same_sector,horizon_low_days,horizon_high_days,cross_weight, ...
        reciprocal_cross_weight,reconstructed_cross_weight_from_modes, ...
        cross_reconstruction_absolute_error,cross_reconstruction_relative_error, ...
        cross_weight_normalized_to_representative_top, ...
        within_low_horizon_weight,within_high_horizon_weight, ...
        within_endpoint_mean,within_endpoint_geometric_mean, ...
        within_endpoint_maximum,strongest_within_horizon_days, ...
        strongest_within_weight,cross_to_endpoint_mean_ratio, ...
        cross_to_endpoint_geometric_mean_ratio,cross_to_endpoint_max_ratio, ...
        cross_to_strongest_within_ratio,cross_minus_endpoint_maximum, ...
        cross_is_larger_than_both_endpoint_within, ...
        selected_cross_fraction_of_pair_maximum, ...
        cross_global_rank,cross_global_tie_count, ...
        cross_global_total_relations,cross_global_top_fraction_percent, ...
        cross_global_prominence_percentile, ...
        cross_block_rank,cross_block_tie_count, ...
        cross_block_total_stock_pairs,cross_block_top_fraction_percent, ...
        cross_block_prominence_percentile, ...
        within_low_rank,within_low_tie_count, ...
        within_low_total_stock_pairs,within_low_top_fraction_percent, ...
        within_low_prominence_percentile, ...
        within_high_rank,within_high_tie_count, ...
        within_high_total_stock_pairs,within_high_top_fraction_percent, ...
        within_high_prominence_percentile, ...
        best_endpoint_within_rank, ...
        best_endpoint_within_prominence_percentile, ...
        cross_block_rank_advantage_over_best_within, ...
        cross_block_prominence_advantage_over_best_within, ...
        cross_block_rank_ratio_to_best_within, ...
        cross_more_prominent_than_both_endpoint_within, ...
        cross_more_prominent_than_at_least_one_endpoint_within, ...
        dominant_mode,dominant_mode_cross_fraction);

    % Append one within-view column for every horizon and one A^(m) column
    % for every latent mode.
    for k = 1:K
        name = matlab.lang.makeValidName( ...
            sprintf('within_%dd_weight',horizons(k)));
        summaryT.(name) = withinByHorizon(:,k);
    end
    for m = 1:r
        nameA = sprintf('A_mode_%d_stock_pair_weight',m);
        nameC = sprintf('mode_%d_cross_contribution',m);
        summaryT.(nameA) = AModeByRelation(:,m);
        summaryT.(nameC) = crossModeContribution(:,m);
    end

    modeRows = modeRows(1:modeCounter);
    modeLong = struct2table(vertcat(modeRows{:}));
    couplingRows = couplingRows(1:couplingCounter);
    couplingLong = struct2table(vertcat(couplingRows{:}));
end


function R = pcm_build_rank_reference(A_supra,within,data)
% Construct the exact reference populations used for ranking.
% - within_weights{k}: all n-choose-2 unordered stock pairs at horizon k.
% - cross_block_weights{k,l}: all n-choose-2 unordered stock pairs in the
%   k-l cross-horizon block, after collapsing reciprocal orientations by
%   their maximum.
% - global_cross_weights: all cross-block populations concatenated.
    n = data.num_stocks;
    K = data.num_views;
    [ii,jj] = find(triu(true(n),1));
    nPairs = numel(ii);

    withinWeights = cell(K,1);
    for k = 1:K
        A = double(within(:,:,k));
        values = A(sub2ind([n n],ii,jj));
        values(~isfinite(values)) = -Inf;
        withinWeights{k} = values(:);
    end

    crossBlockWeights = cell(K,K);
    nViewPairs = K*(K-1)/2;
    globalWeights = zeros(nPairs*nViewPairs,1);
    cursor = 0;
    for k = 1:K
        Ik = (k-1)*n+(1:n);
        for ell = k+1:K
            Il = (ell-1)*n+(1:n);
            C = double(A_supra(Ik,Il));
            forward = C(sub2ind([n n],ii,jj));
            reverse = C(sub2ind([n n],jj,ii));
            values = max(forward,reverse);
            values(~isfinite(values)) = -Inf;
            crossBlockWeights{k,ell} = values(:);
            crossBlockWeights{ell,k} = values(:);
            globalWeights(cursor+(1:nPairs)) = values(:);
            cursor = cursor+nPairs;
        end
    end
    globalWeights = globalWeights(1:cursor);

    R = struct();
    R.within_weights = withinWeights;
    R.cross_block_weights = crossBlockWeights;
    R.global_cross_weights = globalWeights;
    R.n_stock_pairs = nPairs;
    R.n_global_cross_relations = numel(globalWeights);
end


function [rankValue,tieCount] = pcm_descending_competition_rank(value,values)
% Rank 1 is strongest. Equal values receive the same competition rank:
% rank = 1 + number of reference values strictly greater than value.
    values = double(values(:));
    value = double(value);
    if ~isfinite(value)
        rankValue = NaN;
        tieCount = 0;
        return;
    end
    finiteMask = isfinite(values);
    values = values(finiteMask);
    if isempty(values)
        rankValue = NaN;
        tieCount = 0;
        return;
    end
    tolerance = 1e-12*max(1,abs(value));
    rankValue = 1 + nnz(values > value+tolerance);
    tieCount = nnz(abs(values-value) <= tolerance);
end


function value = pcm_top_fraction_percent(rankValue,totalCount)
% Lower is better: 0.10 means approximately the top 0.10 percent.
    if ~isfinite(rankValue) || totalCount<=0
        value = NaN;
    else
        value = 100*rankValue/totalCount;
    end
end


function value = pcm_prominence_percentile(rankValue,totalCount)
% Higher is better: 100 is the strongest possible rank.
    if ~isfinite(rankValue) || totalCount<=0
        value = NaN;
    else
        value = 100*(totalCount-rankValue+1)/totalCount;
    end
end


function T = pcm_extract_ranking_table(summaryT)
% A compact, easy-to-find table containing only the rank comparison.
    if isempty(summaryT)
        T = table();
        return;
    end
    names = { ...
        'representative_id','model_rank','relation_rank', ...
        'relation_key','relation_label', ...
        'stock_a_index','stock_a_ticker','stock_a_sector', ...
        'stock_b_index','stock_b_ticker','stock_b_sector', ...
        'horizon_low_days','horizon_high_days', ...
        'cross_weight','within_low_horizon_weight', ...
        'within_high_horizon_weight', ...
        'cross_global_rank','cross_global_tie_count', ...
        'cross_global_total_relations', ...
        'cross_global_top_fraction_percent', ...
        'cross_global_prominence_percentile', ...
        'cross_block_rank','cross_block_tie_count', ...
        'cross_block_total_stock_pairs', ...
        'cross_block_top_fraction_percent', ...
        'cross_block_prominence_percentile', ...
        'within_low_rank','within_low_tie_count', ...
        'within_low_total_stock_pairs', ...
        'within_low_top_fraction_percent', ...
        'within_low_prominence_percentile', ...
        'within_high_rank','within_high_tie_count', ...
        'within_high_total_stock_pairs', ...
        'within_high_top_fraction_percent', ...
        'within_high_prominence_percentile', ...
        'best_endpoint_within_rank', ...
        'best_endpoint_within_prominence_percentile', ...
        'cross_block_rank_advantage_over_best_within', ...
        'cross_block_prominence_advantage_over_best_within', ...
        'cross_block_rank_ratio_to_best_within', ...
        'cross_more_prominent_than_both_endpoint_within', ...
        'cross_more_prominent_than_at_least_one_endpoint_within'};
    missingNames = names(~ismember(names,summaryT.Properties.VariableNames));
    if ~isempty(missingNames)
        error('Ranking columns were not constructed: %s', ...
            strjoin(missingNames,', '));
    end
    T = summaryT(:,names);
end


function [i,j,hLow,hHigh,key,label] = pcm_relation_identity(T,q,data)
    vars = T.Properties.VariableNames;
    if ismember('stock_a_index',vars)
        i = T.stock_a_index(q);
        j = T.stock_b_index(q);
    elseif ismember('best_source_stock_index',vars)
        i = min(T.best_source_stock_index(q),T.best_target_stock_index(q));
        j = max(T.best_source_stock_index(q),T.best_target_stock_index(q));
    else
        error('Relation table lacks stock indices.');
    end
    stock1 = i;
    stock2 = j;
    i = min(stock1,stock2);
    j = max(stock1,stock2);

    if ismember('horizon_low_days',vars)
        hLow = T.horizon_low_days(q);
        hHigh = T.horizon_high_days(q);
    elseif ismember('best_source_horizon_days',vars)
        hLow = min(T.best_source_horizon_days(q),T.best_target_horizon_days(q));
        hHigh = max(T.best_source_horizon_days(q),T.best_target_horizon_days(q));
    else
        error('Relation table lacks horizon information.');
    end

    if ismember('relation_key',vars)
        key = string(T.relation_key(q));
    else
        key = string(sprintf('s%04d_s%04d_h%04d_h%04d',i,j,hLow,hHigh));
    end
    label = string(sprintf('%s--%s | %dd--%dd', ...
        data.tickers(i),data.tickers(j),hLow,hHigh));
end


function idx = pcm_horizon_index(horizons,h)
    idx = find(horizons==h,1);
    if isempty(idx)
        error('Horizon %g was not found in data.horizons.',h);
    end
end


% ========================================================================
% Complete matrix and edge-table outputs
% ========================================================================

function pcm_save_complete_matrices(fit,within,data,matrixDir)
    horizons = data.horizons;
    r = size(fit.A_state,3);

    writematrix(fit.B,fullfile(matrixDir,'B_view_mode_profiles.csv'));
    writematrix(fit.Gamma,fullfile(matrixDir,'Gamma_stock_copy_coefficients.csv'));
    writematrix(fit.view_graph,fullfile(matrixDir,'learned_view_graph.csv'));

    stockMeta = table((1:data.num_stocks)',data.tickers,data.sector_labels, ...
        'VariableNames',{'stock_index','ticker','sector'});
    writetable(stockMeta,fullfile(matrixDir,'stock_index_metadata.csv'));
    viewMeta = table((1:data.num_views)',data.horizons(:),data.view_labels, ...
        'VariableNames',{'view_index','horizon_days','view_label'});
    writetable(viewMeta,fullfile(matrixDir,'view_index_metadata.csv'));

    for k = 1:data.num_views
        fileName = sprintf('within_A_h%03dd.csv',horizons(k));
        writematrix(within(:,:,k),fullfile(matrixDir,fileName));
        sortedName = sprintf('within_A_h%03dd_sector_sorted.csv',horizons(k));
        writematrix(within(data.sector_order,data.sector_order,k), ...
            fullfile(matrixDir,sortedName));
    end

    for m = 1:r
        fileName = sprintf('A_mode_%02d.csv',m);
        writematrix(fit.A_state(:,:,m),fullfile(matrixDir,fileName));
        sortedName = sprintf('A_mode_%02d_sector_sorted.csv',m);
        writematrix(fit.A_state(data.sector_order,data.sector_order,m), ...
            fullfile(matrixDir,sortedName));
    end
end


function T = pcm_top_within_edges(within,data,maxPerHorizon, ...
        representativeId,rankValue)
    rows = cell(data.num_views*maxPerHorizon,1);
    counter = 0;
    for k = 1:data.num_views
        A = within(:,:,k);
        [ii,jj] = find(triu(true(data.num_stocks),1));
        weights = A(sub2ind(size(A),ii,jj));
        [weights,order] = sort(weights,'descend');
        ii = ii(order); jj = jj(order);
        keep = find(weights>0 & isfinite(weights),min(maxPerHorizon,numel(weights)),'first');
        for qq = 1:numel(keep)
            z = keep(qq);
            counter = counter+1;
            rec = struct();
            rec.representative_id = string(representativeId);
            rec.model_rank = rankValue;
            rec.horizon_days = data.horizons(k);
            rec.within_rank = qq;
            rec.stock_a_index = ii(z);
            rec.stock_a_ticker = data.tickers(ii(z));
            rec.stock_a_sector = data.sector_labels(ii(z));
            rec.stock_b_index = jj(z);
            rec.stock_b_ticker = data.tickers(jj(z));
            rec.stock_b_sector = data.sector_labels(jj(z));
            rec.same_sector = rec.stock_a_sector==rec.stock_b_sector;
            rec.within_weight = weights(z);
            rows{counter} = rec;
        end
    end
    rows = rows(1:counter);
    if isempty(rows), T=table(); else, T=struct2table(vertcat(rows{:})); end
end


function T = pcm_top_mode_edges(Astate,data,maxPerMode, ...
        representativeId,rankValue)
    r = size(Astate,3);
    rows = cell(r*maxPerMode,1);
    counter = 0;
    for m = 1:r
        A = Astate(:,:,m);
        [ii,jj] = find(triu(true(data.num_stocks),1));
        weights = A(sub2ind(size(A),ii,jj));
        [weights,order] = sort(weights,'descend');
        ii = ii(order); jj = jj(order);
        keep = find(weights>0 & isfinite(weights),min(maxPerMode,numel(weights)),'first');
        for qq = 1:numel(keep)
            z = keep(qq);
            counter = counter+1;
            rec = struct();
            rec.representative_id = string(representativeId);
            rec.model_rank = rankValue;
            rec.mode = m;
            rec.mode_edge_rank = qq;
            rec.stock_a_index = ii(z);
            rec.stock_a_ticker = data.tickers(ii(z));
            rec.stock_a_sector = data.sector_labels(ii(z));
            rec.stock_b_index = jj(z);
            rec.stock_b_ticker = data.tickers(jj(z));
            rec.stock_b_sector = data.sector_labels(jj(z));
            rec.same_sector = rec.stock_a_sector==rec.stock_b_sector;
            rec.A_mode_weight = weights(z);
            rows{counter} = rec;
        end
    end
    rows = rows(1:counter);
    if isempty(rows), T=table(); else, T=struct2table(vertcat(rows{:})); end
end


% ========================================================================
% Representative-level figures
% ========================================================================

function pcm_plot_cross_vs_within(T,figureDir,representativeId,cfg)
    if isempty(T), return; end
    n = height(T);
    f = pcm_new_figure([1700,max(760,48*n+280)]);
    ax = axes('Parent',f,'Position',[0.35 0.10 0.62 0.82]);
    values = [T.within_low_horizon_weight, ...
        T.within_high_horizon_weight,T.cross_weight];
    barh(ax,values(end:-1:1,:),'grouped');
    labels = cellstr(T.relation_label(end:-1:1));
    set(ax,'YTick',1:n,'YTickLabel',labels,'FontSize',9);
    grid(ax,'on'); box(ax,'on');
    xlabel(ax,'learned edge weight');
    legend(ax,{'within at lower horizon','within at higher horizon', ...
        'selected cross-horizon'},'Location','southoutside', ...
        'Orientation','horizontal');
    title(ax,sprintf('%s: cross-horizon edge versus same-pair within-view edges', ...
        representativeId),'Interpreter','none');
    pcm_export_figure(f,fullfile(figureDir, ...
        '01_top_cross_vs_within_edge_weights.png'),cfg);

    f = pcm_new_figure([1550,650]);
    tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
    ax = nexttile(tl,1);
    scatter(ax,T.within_endpoint_maximum,T.cross_weight,70, ...
        T.relation_rank,'filled');
    hold(ax,'on');
    lim = max([T.within_endpoint_maximum;T.cross_weight]);
    plot(ax,[0 lim],[0 lim],'k--','LineWidth',1.2);
    grid(ax,'on'); box(ax,'on'); axis(ax,'square');
    xlabel(ax,'max within weight at the two endpoint horizons');
    ylabel(ax,'cross-horizon weight');
    title(ax,'Above diagonal = cross exceeds both endpoint within edges');
    colorbar(ax);

    ax = nexttile(tl,2);
    barh(ax,T.cross_to_endpoint_max_ratio(end:-1:1));
    set(ax,'YTick',1:n,'YTickLabel',labels,'FontSize',8);
    xline(ax,1,'k--','LineWidth',1.2);
    grid(ax,'on'); box(ax,'on');
    xlabel(ax,'cross / max(endpoint within)');
    title(ax,'Cross-horizon specificity ratio');
    sgtitle(tl,sprintf('%s: cross-horizon specificity diagnostics', ...
        representativeId),'Interpreter','none');
    pcm_export_figure(f,fullfile(figureDir, ...
        '02_cross_horizon_specificity_ratios.png'),cfg);
end


function pcm_plot_cross_vs_within_ranks(T,figureDir,representativeId,cfg)
% Compare rank prominence inside the selected cross block against the same
% stock pair's ranks in the two endpoint within-view graphs.
    if isempty(T), return; end
    n = height(T);
    labels = cellstr(T.relation_label(end:-1:1));

    f = pcm_new_figure([1750,max(760,50*n+300)]);
    ax = axes('Parent',f,'Position',[0.35 0.12 0.62 0.80]);
    values = [T.within_low_prominence_percentile, ...
        T.within_high_prominence_percentile, ...
        T.cross_block_prominence_percentile];
    barh(ax,values(end:-1:1,:),'grouped');
    set(ax,'YTick',1:n,'YTickLabel',labels,'FontSize',8);
    xlim(ax,[0 100]); grid(ax,'on'); box(ax,'on');
    xlabel(ax,'prominence percentile within the corresponding block (higher is better)');
    legend(ax,{'within at lower horizon','within at higher horizon', ...
        'selected cross-horizon block'},'Location','southoutside', ...
        'Orientation','horizontal');
    title(ax,sprintf('%s: same stock-pair rank comparison', ...
        representativeId),'Interpreter','none');
    pcm_export_figure(f,fullfile(figureDir, ...
        '02b_cross_vs_within_rank_percentiles.png'),cfg);

    f = pcm_new_figure([1550,650]);
    tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');

    ax = nexttile(tl,1);
    scatter(ax,T.best_endpoint_within_prominence_percentile, ...
        T.cross_block_prominence_percentile,80,T.relation_rank,'filled');
    hold(ax,'on');
    plot(ax,[0 100],[0 100],'k--','LineWidth',1.2);
    xlim(ax,[0 100]); ylim(ax,[0 100]); axis(ax,'square');
    grid(ax,'on'); box(ax,'on');
    xlabel(ax,'best endpoint within-view prominence percentile');
    ylabel(ax,'cross-block prominence percentile');
    title(ax,'Above diagonal = relatively more prominent across horizons');
    colorbar(ax);

    ax = nexttile(tl,2);
    advantage = T.cross_block_prominence_advantage_over_best_within;
    barh(ax,advantage(end:-1:1));
    set(ax,'YTick',1:n,'YTickLabel',labels,'FontSize',8);
    xline(ax,0,'k--','LineWidth',1.2);
    grid(ax,'on'); box(ax,'on');
    xlabel(ax,'cross-block percentile minus best endpoint-within percentile');
    title(ax,'Relative rank advantage of the cross-horizon relation');

    sgtitle(tl,sprintf('%s: cross-vs-within rank diagnostics', ...
        representativeId),'Interpreter','none');
    pcm_export_figure(f,fullfile(figureDir, ...
        '02c_cross_vs_within_rank_advantage.png'),cfg);
end


function pcm_plot_pair_coupling_grid(pairMatrices,T,data, ...
        figureDir,representativeId,cfg)
    if isempty(T), return; end
    n = height(T);
    nCols = min(5,ceil(sqrt(n)));
    nRows = ceil(n/nCols);
    f = pcm_new_figure([400*nCols,370*nRows+100]);
    tl = tiledlayout(f,nRows,nCols,'TileSpacing','compact','Padding','compact');
    for q = 1:n
        ax = nexttile(tl,q);
        P = pairMatrices{q};
        imagesc(ax,P/max(max(P(:)),eps),[0 1]);
        axis(ax,'image'); set(ax,'YDir','normal');
        set(ax,'XTick',1:data.num_views,'YTick',1:data.num_views, ...
            'XTickLabel',cellstr(data.view_labels), ...
            'YTickLabel',cellstr(data.view_labels), ...
            'XTickLabelRotation',45,'FontSize',7);
        title(ax,sprintf('#%d %s--%s',T.relation_rank(q), ...
            T.stock_a_ticker(q),T.stock_b_ticker(q)), ...
            'Interpreter','none','FontSize',9);
        hLow = T.horizon_low_days(q);
        hHigh = T.horizon_high_days(q);
        k = pcm_horizon_index(data.horizons,hLow);
        ell = pcm_horizon_index(data.horizons,hHigh);
        hold(ax,'on');
        rectangle('Parent',ax,'Position',[ell-0.5,k-0.5,1,1], ...
            'EdgeColor','r','LineWidth',1.5);
        rectangle('Parent',ax,'Position',[k-0.5,ell-0.5,1,1], ...
            'EdgeColor','r','LineWidth',1.5);
        if q==n, colorbar(ax); end
    end
    sgtitle(tl,sprintf(['%s: complete stock-pair horizon-coupling matrices ', ...
        '(each normalized by its own maximum)'],representativeId), ...
        'Interpreter','none');
    pcm_export_figure(f,fullfile(figureDir, ...
        '03_top_pair_complete_horizon_coupling_matrices.png'),cfg);
end


function pcm_plot_within_graphs(within,data,figureDir,representativeId,cfg)
    K = data.num_views;
    nCols = ceil(sqrt(K));
    nRows = ceil(K/nCols);
    f = pcm_new_figure([540*nCols,480*nRows]);
    tl = tiledlayout(f,nRows,nCols,'TileSpacing','compact','Padding','compact');
    commonMax = max(within(:));
    commonMax = max(commonMax,eps);
    for k = 1:K
        ax = nexttile(tl,k);
        A = within(data.sector_order,data.sector_order,k);
        imagesc(ax,A,[0 commonMax]); axis(ax,'image');
        pcm_draw_sector_boundaries(ax,data.sector_boundaries);
        title(ax,data.view_labels(k));
        xlabel(ax,'stocks sorted by sector');
        ylabel(ax,'stocks sorted by sector');
        if k==K, colorbar(ax); end
    end
    sgtitle(tl,sprintf('%s: complete learned within-view graphs', ...
        representativeId),'Interpreter','none');
    pcm_export_figure(f,fullfile(figureDir, ...
        '04_all_within_view_adjacencies.png'),cfg);
end


function pcm_plot_mode_graphs(fit,data,figureDir,representativeId,cfg)
    r = size(fit.A_state,3);
    f = pcm_new_figure([520*r,900]);
    tl = tiledlayout(f,2,r,'TileSpacing','compact','Padding','compact');
    commonMax = max(fit.A_state(:));
    commonMax = max(commonMax,eps);
    for m = 1:r
        ax = nexttile(tl,m);
        bar(ax,1:data.num_views,fit.B(:,m));
        set(ax,'XTick',1:data.num_views, ...
            'XTickLabel',cellstr(data.view_labels), ...
            'XTickLabelRotation',40);
        grid(ax,'on'); box(ax,'on');
        ylabel(ax,'B weight');
        title(ax,sprintf('Mode %d view profile',m));

        ax = nexttile(tl,r+m);
        A = fit.A_state(data.sector_order,data.sector_order,m);
        imagesc(ax,A,[0 commonMax]); axis(ax,'image');
        pcm_draw_sector_boundaries(ax,data.sector_boundaries);
        xlabel(ax,'stocks sorted by sector');
        ylabel(ax,'stocks sorted by sector');
        title(ax,sprintf('A^{(%d)} stock graph',m));
        if m==r, colorbar(ax); end
    end
    sgtitle(tl,sprintf('%s: latent mode profiles and A^{(m)} graphs', ...
        representativeId),'Interpreter','none');
    pcm_export_figure(f,fullfile(figureDir, ...
        '05_all_mode_A_graphs_and_B_profiles.png'),cfg);
end


function pcm_plot_mode_decomposition(modeLong,T,figureDir, ...
        representativeId,rankValue,cfg)
    if isempty(modeLong) || isempty(T), return; end
    n = height(T);
    contribution = zeros(n,rankValue);
    Aweight = zeros(n,rankValue);
    for q = 1:n
        rows = modeLong.relation_rank==T.relation_rank(q);
        M = modeLong(rows,:);
        for z = 1:height(M)
            m = M.mode(z);
            contribution(q,m) = M.cross_horizon_mode_contribution(z);
            Aweight(q,m) = M.A_mode_stock_pair_weight(z);
        end
    end

    f = pcm_new_figure([1600,850]);
    tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
    ax = nexttile(tl,1);
    imagesc(ax,contribution); colorbar(ax);
    set(ax,'XTick',1:rankValue,'YTick',1:n, ...
        'YTickLabel',cellstr(T.relation_label),'FontSize',8);
    xlabel(ax,'mode'); ylabel(ax,'top cross relation');
    title(ax,'Mode contribution to selected cross-horizon weight');

    ax = nexttile(tl,2);
    imagesc(ax,Aweight); colorbar(ax);
    set(ax,'XTick',1:rankValue,'YTick',1:n, ...
        'YTickLabel',cellstr(T.relation_label),'FontSize',8);
    xlabel(ax,'mode'); ylabel(ax,'top cross relation');
    title(ax,'A^{(m)} stock-pair coefficient');
    sgtitle(tl,sprintf('%s: latent-mode decomposition of top relations', ...
        representativeId),'Interpreter','none');
    pcm_export_figure(f,fullfile(figureDir, ...
        '06_top_relations_mode_decomposition.png'),cfg);
end


function pcm_save_individual_pair_figures(pairMatrices,T,modeLong,data, ...
        pairFigureDir,representativeId,cfg)
    n = min(height(T),cfg.max_pair_figures);
    for q = 1:n
        rows = modeLong.relation_rank==T.relation_rank(q);
        M = modeLong(rows,:);
        P = pairMatrices{q};
        withinVector = diag(P);

        f = pcm_new_figure([1650,520]);
        tl = tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
        ax = nexttile(tl,1);
        imagesc(ax,P/max(max(P(:)),eps),[0 1]);
        axis(ax,'image'); set(ax,'YDir','normal'); colorbar(ax);
        set(ax,'XTick',1:data.num_views,'YTick',1:data.num_views, ...
            'XTickLabel',cellstr(data.view_labels), ...
            'YTickLabel',cellstr(data.view_labels), ...
            'XTickLabelRotation',40);
        xlabel(ax,'horizon of stock B'); ylabel(ax,'horizon of stock A');
        title(ax,'Complete pair coupling (normalized)');

        ax = nexttile(tl,2);
        plot(ax,data.horizons,withinVector,'-o','LineWidth',1.7);
        hold(ax,'on');
        yline(ax,T.cross_weight(q),'--','LineWidth',1.3);
        set(ax,'XScale','log','XTick',data.horizons, ...
            'XTickLabel',string(data.horizons));
        grid(ax,'on'); box(ax,'on');
        xlabel(ax,'horizon (trading days)'); ylabel(ax,'edge weight');
        legend(ax,{'same-pair within-view weight', ...
            'selected cross-horizon weight'},'Location','best');
        title(ax,'Within-view profile versus cross weight');

        ax = nexttile(tl,3);
        vals = [M.A_mode_stock_pair_weight, ...
            M.cross_horizon_mode_contribution];
        bar(ax,M.mode,vals,'grouped');
        grid(ax,'on'); box(ax,'on');
        xlabel(ax,'mode'); ylabel(ax,'weight');
        legend(ax,{'A^{(m)} pair coefficient', ...
            'cross contribution'},'Location','best');
        title(ax,'Mode decomposition');

        sgtitle(tl,sprintf('%s | #%d %s',representativeId, ...
            T.relation_rank(q),T.relation_label(q)), ...
            'Interpreter','none');
        safeLabel = pcm_safe_filename(char(T.relation_label(q)));
        pcm_export_figure(f,fullfile(pairFigureDir, ...
            sprintf('rank_%02d_%s.png',T.relation_rank(q),safeLabel)),cfg);
    end
end


% ========================================================================
% Across-representative aggregation
% ========================================================================

function T = pcm_relation_consensus(allRelations,nRepresentatives)
    if isempty(allRelations)
        T = table();
        return;
    end
    keys = unique(allRelations.relation_key,'stable');
    rows = cell(numel(keys),1);
    for q = 1:numel(keys)
        R = allRelations(allRelations.relation_key==keys(q),:);
        rec = struct();
        rec.relation_key = keys(q);
        rec.relation_label = R.relation_label(1);
        rec.stock_a_ticker = R.stock_a_ticker(1);
        rec.stock_a_sector = R.stock_a_sector(1);
        rec.stock_b_ticker = R.stock_b_ticker(1);
        rec.stock_b_sector = R.stock_b_sector(1);
        rec.horizon_low_days = R.horizon_low_days(1);
        rec.horizon_high_days = R.horizon_high_days(1);
        rec.representative_count = numel(unique(R.representative_id));
        rec.representative_frequency = rec.representative_count / ...
            max(nRepresentatives,1);
        rec.mean_normalized_cross_weight = mean( ...
            R.cross_weight_normalized_to_representative_top,'omitnan');
        rec.median_normalized_cross_weight = median( ...
            R.cross_weight_normalized_to_representative_top,'omitnan');
        rec.mean_cross_to_endpoint_max_ratio = mean( ...
            R.cross_to_endpoint_max_ratio,'omitnan');
        rec.median_cross_to_endpoint_max_ratio = median( ...
            R.cross_to_endpoint_max_ratio,'omitnan');
        rec.mean_cross_to_strongest_within_ratio = mean( ...
            R.cross_to_strongest_within_ratio,'omitnan');
        rec.fraction_cross_exceeds_both_endpoint_within = mean( ...
            double(R.cross_is_larger_than_both_endpoint_within),'omitnan');
        rec.mean_selected_cross_fraction_of_pair_maximum = mean( ...
            R.selected_cross_fraction_of_pair_maximum,'omitnan');
        rec.mean_relation_rank = mean(R.relation_rank,'omitnan');
        rec.ranks_present = strjoin(string(unique(R.model_rank))',"|");
        rec.representatives = strjoin(unique(R.representative_id)',"|");
        rec.consensus_cross_specificity_score = ...
            rec.representative_frequency * ...
            rec.mean_normalized_cross_weight * ...
            rec.mean_selected_cross_fraction_of_pair_maximum;
        rows{q} = rec;
    end
    T = struct2table(vertcat(rows{:}));
    T = sortrows(T,{'representative_frequency', ...
        'consensus_cross_specificity_score', ...
        'mean_normalized_cross_weight'}, ...
        {'descend','descend','descend'});
    T.consensus_rank = (1:height(T))';
    T = movevars(T,'consensus_rank','Before',1);
end


function T = pcm_stock_pair_consensus(allRelations,nRepresentatives)
    if isempty(allRelations)
        T = table();
        return;
    end
    stockKeys = allRelations.stock_a_ticker + "--" + ...
        allRelations.stock_b_ticker;
    keys = unique(stockKeys,'stable');
    rows = cell(numel(keys),1);
    for q = 1:numel(keys)
        R = allRelations(stockKeys==keys(q),:);
        rec = struct();
        rec.stock_pair = keys(q);
        rec.stock_a_ticker = R.stock_a_ticker(1);
        rec.stock_a_sector = R.stock_a_sector(1);
        rec.stock_b_ticker = R.stock_b_ticker(1);
        rec.stock_b_sector = R.stock_b_sector(1);
        rec.representative_count = numel(unique(R.representative_id));
        rec.representative_frequency = rec.representative_count / ...
            max(nRepresentatives,1);
        rec.number_of_distinct_horizon_pairs = numel(unique( ...
            string(R.horizon_low_days)+"--"+string(R.horizon_high_days)));
        rec.horizon_pairs = strjoin(unique(string(R.horizon_low_days)+ ...
            "--"+string(R.horizon_high_days))',"|");
        rec.mean_normalized_cross_weight = mean( ...
            R.cross_weight_normalized_to_representative_top,'omitnan');
        rec.mean_cross_to_endpoint_max_ratio = mean( ...
            R.cross_to_endpoint_max_ratio,'omitnan');
        rec.mean_selected_cross_fraction_of_pair_maximum = mean( ...
            R.selected_cross_fraction_of_pair_maximum,'omitnan');
        rec.multiscale_pair_score = rec.representative_frequency * ...
            rec.mean_normalized_cross_weight;
        rec.ranks_present = strjoin(string(unique(R.model_rank))',"|");
        rows{q} = rec;
    end
    T = struct2table(vertcat(rows{:}));
    T = sortrows(T,{'representative_frequency','multiscale_pair_score'}, ...
        {'descend','descend'});
    T.stock_pair_rank = (1:height(T))';
    T = movevars(T,'stock_pair_rank','Before',1);
end


function pcm_plot_across_representative_consensus(T,summaryDir,cfg)
    if isempty(T), return; end
    n = min(20,height(T));
    R = T(1:n,:);
    f = pcm_new_figure([1700,max(760,45*n+250)]);
    ax = axes('Parent',f,'Position',[0.34 0.12 0.62 0.78]);
    values = [R.representative_frequency, ...
        R.mean_normalized_cross_weight, ...
        min(R.mean_cross_to_endpoint_max_ratio,2)/2];
    barh(ax,values(end:-1:1,:),'grouped');
    set(ax,'YTick',1:n,'YTickLabel', ...
        cellstr(R.relation_label(end:-1:1)),'FontSize',8);
    grid(ax,'on'); box(ax,'on');
    xlabel(ax,'scaled diagnostic value');
    legend(ax,{'representative frequency', ...
        'mean normalized cross weight', ...
        'cross/endpoint-within ratio divided by 2 (clipped)'}, ...
        'Location','southoutside','Orientation','horizontal');
    title(ax,'Across-representative cross relation recurrence and specificity');
    pcm_export_figure(f,fullfile(summaryDir, ...
        '01_cross_relation_recurrence_and_specificity.png'),cfg);
end


function pcm_plot_stock_pair_multiscale_summary(T,summaryDir,cfg)
    if isempty(T), return; end
    n = min(20,height(T));
    R = T(1:n,:);
    f = pcm_new_figure([1500,max(700,42*n+240)]);
    ax = axes('Parent',f,'Position',[0.28 0.12 0.68 0.78]);
    barh(ax,R.multiscale_pair_score(end:-1:1));
    set(ax,'YTick',1:n,'YTickLabel', ...
        cellstr(R.stock_pair(end:-1:1)),'FontSize',9);
    grid(ax,'on'); box(ax,'on');
    xlabel(ax,'representative frequency x mean normalized cross weight');
    title(ax,'Stable stock pairs across ranks, representatives, and horizons');
    pcm_export_figure(f,fullfile(summaryDir, ...
        '02_stock_pair_multiscale_recurrence.png'),cfg);
end


% ========================================================================
% Small utilities
% ========================================================================

function row = pcm_empty_representative_row()
    row = struct('representative_id','', ...
        'rank',NaN,'source_fit_file','', ...
        'number_of_top_relations',NaN, ...
        'mean_cross_to_endpoint_max_ratio',NaN, ...
        'median_cross_to_endpoint_max_ratio',NaN, ...
        'fraction_cross_exceeds_both_endpoint_within',NaN, ...
        'mean_selected_cross_fraction_of_pair_maximum',NaN, ...
        'stationarity_residual',NaN,'objective',NaN);
end


function row = pcm_make_representative_row(id,rankValue,fitFile,T,fit)
    row = pcm_empty_representative_row();
    row.representative_id = id;
    row.rank = rankValue;
    row.source_fit_file = fitFile;
    row.number_of_top_relations = height(T);
    row.mean_cross_to_endpoint_max_ratio = mean( ...
        T.cross_to_endpoint_max_ratio,'omitnan');
    row.median_cross_to_endpoint_max_ratio = median( ...
        T.cross_to_endpoint_max_ratio,'omitnan');
    row.fraction_cross_exceeds_both_endpoint_within = mean( ...
        double(T.cross_is_larger_than_both_endpoint_within),'omitnan');
    row.mean_selected_cross_fraction_of_pair_maximum = mean( ...
        T.selected_cross_fraction_of_pair_maximum,'omitnan');
    if isfield(fit,'stationarity_combined')
        row.stationarity_residual = fit.stationarity_combined;
    end
    if isfield(fit,'objective')
        row.objective = fit.objective;
    end
end


function T = pcm_vertcat_nonempty(C)
% Vertically concatenate nonempty tables after aligning their variables.
% This is required because r=3 and r=4 relation tables contain different
% numbers of mode-specific columns. Missing variables are filled with a
% type-compatible missing value (NaN for floating-point mode columns).
    nonempty = ~cellfun(@isempty,C);
    C = C(nonempty);
    if isempty(C)
        T = table();
        return;
    end

    % Union of variable names, preserving first occurrence order.
    allNames = string.empty(1,0);
    for q = 1:numel(C)
        names = string(C{q}.Properties.VariableNames);
        for j = 1:numel(names)
            if ~any(allNames==names(j))
                allNames(end+1) = names(j); %#ok<AGROW>
            end
        end
    end

    % Add every missing variable using an exemplar from another table.
    for q = 1:numel(C)
        nRows = height(C{q});
        current = string(C{q}.Properties.VariableNames);
        missingNames = allNames(~ismember(allNames,current));
        for j = 1:numel(missingNames)
            name = char(missingNames(j));
            exemplar = [];
            found = false;
            for z = 1:numel(C)
                if ismember(name,C{z}.Properties.VariableNames)
                    exemplar = C{z}.(name);
                    found = true;
                    break;
                end
            end
            if ~found
                error('Could not find exemplar for table variable %s.',name);
            end
            C{q}.(name) = pcm_missing_like(exemplar,nRows);
        end
        C{q} = C{q}(:,cellstr(allNames));
    end

    T = vertcat(C{:});
end


function value = pcm_missing_like(exemplar,nRows)
    sz = size(exemplar);
    if isempty(sz)
        sz = [0 1];
    end
    sz(1) = nRows;

    if isfloat(exemplar)
        value = nan(sz,'like',exemplar);
    elseif isinteger(exemplar)
        value = zeros(sz,'like',exemplar);
    elseif islogical(exemplar)
        value = false(sz);
    elseif isstring(exemplar)
        value = strings(sz);
        value(:) = missing;
    elseif iscell(exemplar)
        value = cell(sz);
    elseif iscategorical(exemplar)
        value = categorical(strings(sz));
    elseif isdatetime(exemplar)
        value = NaT(sz);
    elseif isduration(exemplar)
        value = seconds(nan(sz));
    elseif ischar(exemplar)
        value = repmat(' ',sz);
    else
        try
            value = repmat(missing,sz);
        catch
            error('Unsupported table variable class for alignment: %s', ...
                class(exemplar));
        end
    end
end


function A = pcm_symmetrize_nonnegative(A)
    A = max(0,0.5*(double(A)+double(A')));
end


function f = pcm_new_figure(position)
    f = figure('Visible','off','Color','w','Position',[50 50 position]);
end


function pcm_export_figure(f,filePath,cfg)
    [folder,~,~] = fileparts(filePath);
    pcm_ensure_folder(folder);
    try
        exportgraphics(f,filePath,'Resolution',cfg.figure_resolution);
    catch
        print(f,filePath,'-dpng',sprintf('-r%d',cfg.figure_resolution));
    end
    close(f);
end


function pcm_draw_sector_boundaries(ax,boundaries)
    hold(ax,'on');
    for q = 1:numel(boundaries)
        xline(ax,boundaries(q),'w:','LineWidth',0.7);
        yline(ax,boundaries(q),'w:','LineWidth',0.7);
    end
end


function pcm_ensure_folder(folder)
    if ~isfolder(folder)
        mkdir(folder);
    end
end


function safe = pcm_safe_filename(name)
    safe = regexprep(name,'[^A-Za-z0-9_\-]+','_');
    safe = regexprep(safe,'_+','_');
    safe = regexprep(safe,'^_|_$','');
    if isempty(safe), safe='relation'; end
end


function out = pcm_merge_struct(base,override)
    out = base;
    if isempty(override), return; end
    names = fieldnames(override);
    for q = 1:numel(names)
        out.(names{q}) = override.(names{q});
    end
end


function files = pcm_find_files_recursive(rootFolder,targetName)
    files = {};
    if ~isfolder(rootFolder), return; end
    listing = dir(rootFolder);
    for q = 1:numel(listing)
        name = listing(q).name;
        if strcmp(name,'.') || strcmp(name,'..'), continue; end
        fullPath = fullfile(rootFolder,name);
        if listing(q).isdir
            child = pcm_find_files_recursive(fullPath,targetName);
            files = [files;child]; %#ok<AGROW>
        elseif strcmp(name,targetName)
            files{end+1,1} = fullPath; %#ok<AGROW>
        end
    end
end


function folders = pcm_find_folders_recursive(rootFolder,targetName,maxDepth)
    folders = {};
    if maxDepth<0 || ~isfolder(rootFolder), return; end
    listing = dir(rootFolder);
    for q = 1:numel(listing)
        name = listing(q).name;
        if ~listing(q).isdir || strcmp(name,'.') || strcmp(name,'..')
            continue;
        end
        fullPath = fullfile(rootFolder,name);
        if strcmp(name,targetName)
            folders{end+1,1} = fullPath; %#ok<AGROW>
        elseif maxDepth>0
            child = pcm_find_folders_recursive(fullPath,targetName,maxDepth-1);
            folders = [folders;child]; %#ok<AGROW>
        end
    end
end


function pcm_write_readme(filePath,runRoot,dataFile,cfg,nFits)
    fid = fopen(filePath,'w');
    if fid<0, return; end
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid,'MOSAIC finance cross-vs-within and mode post-processing\n');
    fprintf(fid,'======================================================\n\n');
    fprintf(fid,'Source run: %s\n',runRoot);
    fprintf(fid,'Dataset   : %s\n',dataFile);
    fprintf(fid,'Representative fits processed: %d\n',nFits);
    fprintf(fid,'Top cross relations per fit   : %d\n\n', ...
        cfg.top_cross_relations);
    fprintf(fid,['For every selected cross-stock/cross-horizon relation, ', ...
        'the output compares the cross edge with the same stock pair''s ', ...
        'within-view weights at all horizons.\n']);
    fprintf(fid,['The full pair-specific K-by-K horizon coupling matrix is ', ...
        'saved, with its diagonal representing within-view weights and ', ...
        'its off-diagonal entries representing cross-horizon weights.\n']);
    fprintf(fid,['A^(m) matrices, B profiles, Gamma copy coefficients, complete ', ...
        'within-view graphs, top within edges, and mode decompositions ', ...
        'are also saved.\n']);
end
