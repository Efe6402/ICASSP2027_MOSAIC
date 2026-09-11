%% build_mosaic_finance_multihorizon_dataset.m
% =========================================================================
% Prepare a real multi-view S&P 500 dataset for MOSAIC.
%
% EXPECTED FOLDER STRUCTURE
%   <this script>/
%       build_mosaic_finance_multihorizon_dataset.m
%       stocks/
%           sp500-data-2016-2020_converted.csv
%           SP500-sectors.csv
%
% MODELING CHOICE
%   Nodes : stocks
%   Views : return horizons h = [1, 5, 25, 60, 100] trading days
%   Signal coordinates : common calendar ending dates, sampled every 5
%                        trading days
%
% For stock i, horizon h, and common ending date t_q,
%
%   R_i^(h)(q) = log(P_i(t_q)) - log(P_i(t_q-h)).
%
% Each stock is standardized separately within each horizon:
%
%   X_i^(h)(q) = (R_i^(h)(q) - mu_i,h) / sigma_i,h.
%
% The final MOSAIC inputs are
%
%   data.X{k} in R^(n x p), k = 1,...,K,
%
% where rows are stocks and columns are the same ending dates in every view.
%
% IMPORTANT
%   - No correlation adjacency is constructed here.
%   - No absolute correlation, thresholding, fixed density, MST, detoning,
%     or Marchenko-Pastur filtering is applied.
%   - GICS sectors are saved only for interpretation/diagnostics; they are
%     not used to build the signals.
%
% OUTPUT
%   A timestamped folder containing:
%     1) mosaic_finance_multihorizon_dataset.mat
%     2) selected_stocks_and_sectors.csv
%     3) common_ending_dates.csv
%     4) view_diagnostics.csv
%     5) same_stock_cross_horizon_similarity.csv
%     6) same_stock_cross_horizon_rms_distance.csv
%     7) diagnostic plots and a text report
%
% MATLAB REQUIREMENTS
%   Written as a script with local functions (MATLAB R2016b or newer).
% =========================================================================

clear; clc; close all;

%% ============================= CONFIGURATION ============================
cfg = struct();

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

cfg.stocksDir = fullfile(scriptDir, 'stocks');
cfg.priceCsv  = fullfile(cfg.stocksDir, 'sp500-data-2016-2020_converted.csv');
cfg.sectorCsv = fullfile(cfg.stocksDir, 'SP500-sectors.csv');

% A frozen ticker list is created on the first successful run. Later runs
% reuse it so that the stock universe and row ordering remain reproducible.
cfg.frozenTickerCsv = fullfile(cfg.stocksDir, 'mosaic_selected_tickers.csv');

cfg.numStocks   = 100;
cfg.horizons    = [1, 20, 45, 70, 120]; %[1, 5, 25, 60, 100]
cfg.endpointStep = 5;                 % trading days between signal columns
cfg.scaling      = 'row_zscore';      % each stock, separately per horizon
cfg.winsorize    = false;             % keep crisis/extreme observations
cfg.saveDPI      = 200;
cfg.heatmapQuantile = 0.995;

% Output folder
stamp = datestr(now, 'yyyymmdd_HHMMSS');
cfg.outputDir = fullfile(scriptDir, ...
    ['mosaic_finance_multihorizon_output_' stamp]);

assert(isfolder(cfg.stocksDir), ...
    'Cannot find the stocks folder: %s', cfg.stocksDir);
assert(isfile(cfg.priceCsv), ...
    'Cannot find the price CSV: %s', cfg.priceCsv);
assert(isfile(cfg.sectorCsv), ...
    'Cannot find the sector CSV: %s', cfg.sectorCsv);

if ~exist(cfg.outputDir, 'dir')
    mkdir(cfg.outputDir);
end

logFile = fullfile(cfg.outputDir, 'run_log.txt');
diary(logFile);
cleanupDiary = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('\n============================================================\n');
fprintf('MOSAIC FINANCE: MULTI-HORIZON SIGNAL PREPARATION\n');
fprintf('============================================================\n');
fprintf('Script directory : %s\n', scriptDir);
fprintf('Stocks directory : %s\n', cfg.stocksDir);
fprintf('Output directory : %s\n', cfg.outputDir);
fprintf('Horizons         : %s trading days\n', mat2str(cfg.horizons));
fprintf('Endpoint step    : %d trading days\n', cfg.endpointStep);
fprintf('Target stocks    : %d\n\n', cfg.numStocks);

%% =============================== LOAD CSV ===============================
fprintf('Reading price CSV...\n');
Traw = readtable(cfg.priceCsv, 'VariableNamingRule', 'preserve');
assert(width(Traw) >= 2, ...
    'The price CSV must contain a date column and at least one stock column.');

rawDate = Traw{:,1};
priceDates = parse_dates_strict(rawDate);
priceTable = Traw(:,2:end);

[Praw, numericKeep] = table_to_numeric_matrix_with_mask(priceTable);
tickersAll = string(priceTable.Properties.VariableNames);
tickersAll = tickersAll(numericKeep);

assert(size(Praw,1) == numel(priceDates), ...
    'Price rows and date rows are inconsistent.');
assert(size(Praw,2) == numel(tickersAll), ...
    'Price columns and ticker names are inconsistent.');

% Remove rows with unparsed dates, then sort chronologically.
goodDate = ~isnat(priceDates);
if nnz(~goodDate) > 0
    fprintf('Dropping %d rows with invalid dates.\n', nnz(~goodDate));
    priceDates = priceDates(goodDate);
    Praw = Praw(goodDate,:);
end

[priceDates, sortIdx] = sort(priceDates, 'ascend');
Praw = Praw(sortIdx,:);

% Remove duplicate dates while preserving the first chronological entry.
[priceDatesUnique, uniqueIdx] = unique(priceDates, 'stable');
if numel(priceDatesUnique) < numel(priceDates)
    fprintf('Removing %d duplicated date rows.\n', ...
        numel(priceDates) - numel(priceDatesUnique));
    priceDates = priceDatesUnique;
    Praw = Praw(uniqueIdx,:);
end

fprintf('Price date range  : %s to %s\n', ...
    datestr(priceDates(1), 'yyyy-mm-dd'), ...
    datestr(priceDates(end), 'yyyy-mm-dd'));
fprintf('Price matrix size : %d dates x %d candidate columns\n', ...
    size(Praw,1), size(Praw,2));

%% ======================= FIX THE STOCK UNIVERSE =========================
validStock = all(isfinite(Praw),1) & all(Praw > 0,1);
fprintf('Complete positive-price columns: %d / %d\n', ...
    nnz(validStock), numel(validStock));

if isfile(cfg.frozenTickerCsv)
    fprintf('Using frozen ticker list:\n  %s\n', cfg.frozenTickerCsv);
    F = readtable(cfg.frozenTickerCsv, 'VariableNamingRule', 'preserve');
    assert(width(F) >= 1, 'Frozen ticker CSV is empty.');
    frozenTickers = string(F{:,1});
    frozenTickers = frozenTickers(strlength(strtrim(frozenTickers)) > 0);

    [found, selectedIdx] = match_tickers(frozenTickers, tickersAll);
    if any(~found)
        missing = strjoin(frozenTickers(~found), ', ');
        error('Frozen tickers not found in price CSV: %s', char(missing));
    end
    if any(~validStock(selectedIdx))
        bad = strjoin(tickersAll(selectedIdx(~validStock(selectedIdx))), ', ');
        error('Some frozen tickers do not have complete positive prices: %s', char(bad));
    end
    if numel(selectedIdx) ~= cfg.numStocks
        error(['Frozen ticker file contains %d tickers, but cfg.numStocks=%d. ' ...
               'Update one of them deliberately.'], ...
               numel(selectedIdx), cfg.numStocks);
    end
else
    validIdx = find(validStock);
    if numel(validIdx) < cfg.numStocks
        error('Only %d valid stocks are available; %d were requested.', ...
            numel(validIdx), cfg.numStocks);
    end
    selectedIdx = validIdx(1:cfg.numStocks);
    frozenTickers = tickersAll(selectedIdx);

    fprintf('Creating frozen ticker list:\n  %s\n', cfg.frozenTickerCsv);
    writetable(table(frozenTickers(:), 'VariableNames', {'Ticker'}), ...
        cfg.frozenTickerCsv);
end

P = Praw(:,selectedIdx);
tickers = tickersAll(selectedIdx);
n = size(P,2);

fprintf('Selected stock count: n = %d\n', n);
fprintf('First selected tickers: %s\n', ...
    char(strjoin(tickers(1:min(15,n)), ', ')));

%% ============================= LOAD SECTORS =============================
[sectorLabels, sectorKnown] = load_sector_labels(cfg.sectorCsv, tickers);
[sectorOrder, sectorNames, sectorBoundaries] = ...
    make_sector_order(sectorLabels);

fprintf('Known GICS sectors: %d / %d stocks\n', nnz(sectorKnown), n);
fprintf('Distinct sector labels: %d\n', numel(sectorNames));
for s = 1:numel(sectorNames)
    fprintf('  %2d) %-25s %3d stocks\n', s, char(sectorNames(s)), ...
        nnz(sectorLabels == sectorNames(s)));
end

selectedStockTable = table((1:n)', tickers(:), sectorLabels(:), sectorKnown(:), ...
    'VariableNames', {'NodeIndex','Ticker','GICSSector','SectorKnown'});
writetable(selectedStockTable, ...
    fullfile(cfg.outputDir, 'selected_stocks_and_sectors.csv'));

%% =========================== DAILY LOG RETURNS ==========================
fprintf('\nComputing log prices and daily log returns...\n');
logP = log(P);
dailyLogReturns = diff(logP,1,1).';     % n x (T_price-1)
dailyReturnDates = priceDates(2:end);

assert(all(isfinite(dailyLogReturns(:))), ...
    'Daily log returns contain non-finite values.');

fprintf('Daily-return matrix: %d stocks x %d trading days\n', ...
    size(dailyLogReturns,1), size(dailyLogReturns,2));

%% ======================= COMMON ENDING DATE GRID ========================
horizons = cfg.horizons(:).';
K = numel(horizons);
Hmax = max(horizons);
numPriceDates = size(P,1);

if numPriceDates <= Hmax
    error('Only %d price dates are available; maximum horizon is %d.', ...
        numPriceDates, Hmax);
end

% Price row t can support an h-day return if t-h >= 1.
endPriceIdx = (Hmax + 1):cfg.endpointStep:numPriceDates;
commonDates = priceDates(endPriceIdx);
p = numel(endPriceIdx);

if p < 50
    warning(['Only %d common signal coordinates remain. Consider reducing ' ...
             'the maximum horizon or endpoint step.'], p);
end

fprintf('\nCommon aligned endpoints:\n');
fprintf('  First ending date : %s\n', datestr(commonDates(1), 'yyyy-mm-dd'));
fprintf('  Last ending date  : %s\n', datestr(commonDates(end), 'yyyy-mm-dd'));
fprintf('  Number of columns : p = %d\n', p);
fprintf('  Endpoint spacing  : %d trading days\n', cfg.endpointStep);

writetable(table((1:p)', endPriceIdx(:), commonDates(:), ...
    'VariableNames', {'CoordinateIndex','PriceRowIndex','EndingDate'}), ...
    fullfile(cfg.outputDir, 'common_ending_dates.csv'));

%% ===================== BUILD MULTI-HORIZON SIGNALS ======================
fprintf('\nBuilding multi-horizon node signals...\n');

X = cell(1,K);                    % standardized signals, n x p
Xraw = cell(1,K);                 % raw cumulative log returns, n x p
normMu = zeros(n,K);
normSigma = zeros(n,K);
viewNames = cell(1,K);

for k = 1:K
    h = horizons(k);
    viewNames{k} = sprintf('%d-day', h);

    % p x n, then transpose to n x p.
    Rh = logP(endPriceIdx,:) - logP(endPriceIdx-h,:);
    Rh = Rh.';

    assert(isequal(size(Rh), [n,p]), ...
        'Unexpected raw signal dimensions for horizon h=%d.', h);
    assert(all(isfinite(Rh(:))), ...
        'Non-finite raw horizon returns found for h=%d.', h);

    if cfg.winsorize
        error(['cfg.winsorize=true is intentionally not implemented in the ' ...
               'primary builder. Keep false, or add a separately documented ' ...
               'robustness pipeline.']);
    end

    mu_h = mean(Rh,2);
    sigma_h = std(Rh,0,2);

    nearConstant = sigma_h < 1e-12 | ~isfinite(sigma_h);
    if any(nearConstant)
        warning('%d near-constant rows at horizon h=%d; scale set to one.', ...
            nnz(nearConstant), h);
        sigma_h(nearConstant) = 1;
    end

    Xh = (Rh - mu_h) ./ sigma_h;

    assert(all(isfinite(Xh(:))), ...
        'Non-finite standardized signals found for h=%d.', h);

    Xraw{k} = Rh;
    X{k} = Xh;
    normMu(:,k) = mu_h;
    normSigma(:,k) = sigma_h;

    fprintf('  h=%3d days -> X{%d}: %d x %d\n', ...
        h, k, size(Xh,1), size(Xh,2));
end

%% ======================== PRELIMINARY DIAGNOSTICS =======================
fprintf('\nComputing preliminary diagnostics...\n');

viewId = (1:K)';
horizonDays = horizons(:);
viewNameColumn = string(viewNames(:));
nStocksColumn = repmat(n,K,1);
nCoordinatesColumn = repmat(p,K,1);
firstDateColumn = repmat(commonDates(1),K,1);
lastDateColumn = repmat(commonDates(end),K,1);

rawGlobalMean = zeros(K,1);
rawGlobalStd = zeros(K,1);
standardizedGlobalMean = zeros(K,1);
standardizedGlobalStd = zeros(K,1);
meanAbsRowMean = zeros(K,1);
meanRowStd = zeros(K,1);
medianRowNorm = zeros(K,1);
medianPairwiseDistance = zeros(K,1);
q99AbsStandardized = zeros(K,1);
maxAbsStandardized = zeros(K,1);
meanLag1Autocorr = zeros(K,1);
medianLag1Autocorr = zeros(K,1);
meanWithinSectorCorr = nan(K,1);
meanBetweenSectorCorr = nan(K,1);
withinMinusBetweenCorr = nan(K,1);

upperMask = triu(true(n),1);
knownPair = sectorKnown & sectorKnown.';
sameSectorPair = knownPair & (sectorLabels == sectorLabels.');
differentSectorPair = knownPair & (sectorLabels ~= sectorLabels.');
sameSectorPair = sameSectorPair & upperMask;
differentSectorPair = differentSectorPair & upperMask;

for k = 1:K
    Rh = Xraw{k};
    Xh = X{k};

    rawGlobalMean(k) = mean(Rh(:));
    rawGlobalStd(k) = std(Rh(:));
    standardizedGlobalMean(k) = mean(Xh(:));
    standardizedGlobalStd(k) = std(Xh(:));

    rowMeans = mean(Xh,2);
    rowStds = std(Xh,0,2);
    meanAbsRowMean(k) = mean(abs(rowMeans));
    meanRowStd(k) = mean(rowStds);
    medianRowNorm(k) = median(sqrt(sum(Xh.^2,2)));

    pairDistances = upper_pairwise_distances(Xh);
    medianPairwiseDistance(k) = median(pairDistances);

    absVals = abs(Xh(:));
    q99AbsStandardized(k) = quantile(absVals,0.99);
    maxAbsStandardized(k) = max(absVals);

    lag1 = rowwise_lag1_correlation(Xh);
    meanLag1Autocorr(k) = mean(lag1,'omitnan');
    medianLag1Autocorr(k) = median(lag1,'omitnan');

    % Since every row is standardized, this is the sample-correlation matrix.
    Ck = (Xh * Xh.') / max(p-1,1);
    Ck = max(min(Ck,1),-1);

    if any(sameSectorPair(:))
        meanWithinSectorCorr(k) = mean(Ck(sameSectorPair),'omitnan');
    end
    if any(differentSectorPair(:))
        meanBetweenSectorCorr(k) = mean(Ck(differentSectorPair),'omitnan');
    end
    withinMinusBetweenCorr(k) = ...
        meanWithinSectorCorr(k) - meanBetweenSectorCorr(k);
end

viewDiagnostics = table( ...
    viewId, horizonDays, viewNameColumn, nStocksColumn, nCoordinatesColumn, ...
    firstDateColumn, lastDateColumn, rawGlobalMean, rawGlobalStd, ...
    standardizedGlobalMean, standardizedGlobalStd, meanAbsRowMean, ...
    meanRowStd, medianRowNorm, medianPairwiseDistance, ...
    q99AbsStandardized, maxAbsStandardized, meanLag1Autocorr, ...
    medianLag1Autocorr, meanWithinSectorCorr, meanBetweenSectorCorr, ...
    withinMinusBetweenCorr, ...
    'VariableNames', { ...
    'view_id','horizon_days','view_name','n_stocks','n_coordinates', ...
    'first_ending_date','last_ending_date','raw_global_mean','raw_global_std', ...
    'standardized_global_mean','standardized_global_std','mean_abs_row_mean', ...
    'mean_row_std','median_row_norm','median_pairwise_distance', ...
    'q99_abs_standardized','max_abs_standardized','mean_lag1_autocorr', ...
    'median_lag1_autocorr','mean_within_sector_corr', ...
    'mean_between_sector_corr','within_minus_between_corr'});

writetable(viewDiagnostics, ...
    fullfile(cfg.outputDir, 'view_diagnostics.csv'));

disp(viewDiagnostics(:, {'horizon_days','n_coordinates', ...
    'mean_row_std','median_pairwise_distance','mean_lag1_autocorr', ...
    'mean_within_sector_corr'}));

%% ================ SAME-STOCK CROSS-HORIZON DIAGNOSTICS =================
sameStockSimilarity = eye(K);
sameStockRmsDistance = zeros(K);

for k = 1:K
    for ell = k:K
        Xk = X{k};
        Xell = X{ell};

        % Sample correlation between the two horizon signals, calculated
        % separately for each stock and then averaged across stocks.
        stockCorr = sum(Xk .* Xell, 2) / max(p-1,1);
        stockCorr = max(min(stockCorr,1),-1);
        avgCorr = mean(stockCorr,'omitnan');

        % Per-stock RMS difference over the aligned ending dates.
        stockRms = sqrt(mean((Xk - Xell).^2,2));
        avgRms = mean(stockRms,'omitnan');

        sameStockSimilarity(k,ell) = avgCorr;
        sameStockSimilarity(ell,k) = avgCorr;
        sameStockRmsDistance(k,ell) = avgRms;
        sameStockRmsDistance(ell,k) = avgRms;
    end
end

similarityTable = array2table(sameStockSimilarity, ...
    'VariableNames', matlab.lang.makeValidName(viewNames), ...
    'RowNames', viewNames);
distanceTable = array2table(sameStockRmsDistance, ...
    'VariableNames', matlab.lang.makeValidName(viewNames), ...
    'RowNames', viewNames);

writetable(similarityTable, ...
    fullfile(cfg.outputDir, 'same_stock_cross_horizon_similarity.csv'), ...
    'WriteRowNames', true);
writetable(distanceTable, ...
    fullfile(cfg.outputDir, 'same_stock_cross_horizon_rms_distance.csv'), ...
    'WriteRowNames', true);

fprintf('\nAverage same-stock cross-horizon correlation:\n');
disp(similarityTable);

%% ============================ BUILD DATA STRUCT =========================
data = struct();

data.description = [ ...
    'S&P 500 multi-horizon MOSAIC dataset. Nodes are stocks; views are ' ...
    'cumulative-return horizons; columns are common calendar ending dates.'];

data.X = X;                          % 1 x K cell, each n x p
data.X_raw = Xraw;                   % unstandardized cumulative log returns
data.horizons = horizons;
data.view_names = viewNames;
data.K = K;
data.n = n;
data.p = p;
data.endpoint_step = cfg.endpointStep;
data.common_dates = commonDates;
data.end_price_indices = endPriceIdx;

data.tickers = tickers;
data.node_names = tickers;
data.sector_labels = sectorLabels;
data.sector_known = sectorKnown;
data.sector_order = sectorOrder;
data.sector_names = sectorNames;
data.sector_boundaries = sectorBoundaries;

data.prices = P;
data.price_dates = priceDates;
data.daily_log_returns = dailyLogReturns;
data.daily_return_dates = dailyReturnDates;

data.normalization = cfg.scaling;
data.norm_mu = normMu;
data.norm_sigma = normSigma;

data.same_stock_cross_horizon_similarity = sameStockSimilarity;
data.same_stock_cross_horizon_rms_distance = sameStockRmsDistance;
data.view_diagnostics = viewDiagnostics;
data.config = cfg;

data.preprocessing_notes = { ...
    'No graph or correlation adjacency is supplied to MOSAIC.'; ...
    'Signals retain their sign; positive and negative co-movement are distinct.'; ...
    'Each stock is z-scored separately within each horizon.'; ...
    'All views use identical calendar ending dates.'; ...
    'Ending dates are sampled every five trading days to reduce redundancy.'; ...
    'Long-horizon columns remain overlapping and are not independent samples.'; ...
    'GICS sectors are used only for diagnostics and interpretation.'};

matPath = fullfile(cfg.outputDir, ...
    'mosaic_finance_multihorizon_dataset.mat');
save(matPath, 'data', '-v7.3');
fprintf('Saved dataset:\n  %s\n', matPath);

%% =============================== PLOTS ==================================
fprintf('Saving diagnostic plots...\n');

plot_cross_horizon_matrices(sameStockSimilarity, ...
    sameStockRmsDistance, viewNames, cfg.outputDir, cfg.saveDPI);

plot_view_diagnostics(viewDiagnostics, viewNames, ...
    cfg.outputDir, cfg.saveDPI);

plot_signal_heatmaps(X, horizons, commonDates, sectorOrder, ...
    sectorBoundaries, cfg.heatmapQuantile, cfg.outputDir, cfg.saveDPI);

%% ============================ TEXT REPORT ===============================
reportPath = fullfile(cfg.outputDir, 'preliminary_diagnostics_report.txt');
write_diagnostic_report(reportPath, cfg, data, viewDiagnostics, ...
    sameStockSimilarity, sameStockRmsDistance);

fprintf('\n============================================================\n');
fprintf('DONE\n');
fprintf('============================================================\n');
fprintf('Final data dimensions: n=%d stocks, K=%d views, p=%d columns\n', ...
    n, K, p);
fprintf('Horizons: %s\n', mat2str(horizons));
fprintf('Output folder:\n  %s\n', cfg.outputDir);
fprintf('Main dataset:\n  %s\n', matPath);
fprintf('\nUse data.X{k} as the MOSAIC signal matrix for view k.\n');

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

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

function d = upper_pairwise_distances(X)
% Euclidean distances between all distinct row pairs, without pdist.
    rowSq = sum(X.^2,2);
    D2 = rowSq + rowSq.' - 2*(X*X.');
    D2 = max(D2,0);
    mask = triu(true(size(X,1)),1);
    d = sqrt(D2(mask));
end

function rho = rowwise_lag1_correlation(X)
% Lag-one Pearson correlation for each row.
    if size(X,2) < 3
        rho = nan(size(X,1),1);
        return;
    end

    A = X(:,1:end-1);
    B = X(:,2:end);
    A = A - mean(A,2);
    B = B - mean(B,2);

    numerator = sum(A.*B,2);
    denominator = sqrt(sum(A.^2,2).*sum(B.^2,2));
    rho = numerator ./ denominator;
    rho(denominator < 1e-14) = NaN;
    rho = max(min(rho,1),-1);
end

function plot_cross_horizon_matrices(S, D, viewNames, outputDir, dpi)
    fig = figure('Visible','off','Color','w', ...
        'Position',[80 80 1300 560]);
    tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    ax1 = nexttile(tl);
    imagesc(ax1,S,[-1 1]);
    axis(ax1,'image');
    colorbar(ax1);
    set(ax1,'XTick',1:numel(viewNames),'XTickLabel',viewNames, ...
        'YTick',1:numel(viewNames),'YTickLabel',viewNames, ...
        'XTickLabelRotation',35);
    title(ax1,'Average same-stock cross-horizon correlation');
    add_matrix_text(ax1,S,'%.2f');

    ax2 = nexttile(tl);
    imagesc(ax2,D);
    axis(ax2,'image');
    colorbar(ax2);
    set(ax2,'XTick',1:numel(viewNames),'XTickLabel',viewNames, ...
        'YTick',1:numel(viewNames),'YTickLabel',viewNames, ...
        'XTickLabelRotation',35);
    title(ax2,'Average same-stock RMS difference');
    add_matrix_text(ax2,D,'%.2f');

    sgtitle(tl,'Cross-horizon signal diagnostics','FontWeight','bold');
    save_figure(fig, fullfile(outputDir, ...
        'cross_horizon_signal_diagnostics.png'), dpi);
    close(fig);
end

function add_matrix_text(ax,M,fmt)
    hold(ax,'on');
    for i = 1:size(M,1)
        for j = 1:size(M,2)
            text(ax,j,i,sprintf(fmt,M(i,j)), ...
                'HorizontalAlignment','center', ...
                'FontSize',9,'FontWeight','bold');
        end
    end
    hold(ax,'off');
end

function plot_view_diagnostics(T, viewNames, outputDir, dpi)
    x = 1:height(T);

    fig = figure('Visible','off','Color','w', ...
        'Position',[80 80 1450 850]);
    tl = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

    ax = nexttile(tl);
    bar(ax,x,T.median_pairwise_distance);
    grid(ax,'on'); box(ax,'on');
    set(ax,'XTick',x,'XTickLabel',viewNames,'XTickLabelRotation',30);
    ylabel(ax,'Euclidean distance');
    title(ax,'Median within-view pairwise distance');

    ax = nexttile(tl);
    bar(ax,x,T.mean_lag1_autocorr);
    grid(ax,'on'); box(ax,'on'); ylim(ax,[-1 1]);
    set(ax,'XTick',x,'XTickLabel',viewNames,'XTickLabelRotation',30);
    ylabel(ax,'Correlation');
    title(ax,'Mean lag-one signal autocorrelation');

    ax = nexttile(tl);
    plot(ax,x,T.mean_within_sector_corr,'-o','LineWidth',1.5); hold(ax,'on');
    plot(ax,x,T.mean_between_sector_corr,'-s','LineWidth',1.5);
    yline(ax,0,'--');
    grid(ax,'on'); box(ax,'on');
    set(ax,'XTick',x,'XTickLabel',viewNames,'XTickLabelRotation',30);
    ylabel(ax,'Mean correlation');
    title(ax,'Sector diagnostic (not used for learning)');
    legend(ax,{'Within sector','Between sector'},'Location','best');

    ax = nexttile(tl);
    plot(ax,x,T.q99_abs_standardized,'-o','LineWidth',1.5); hold(ax,'on');
    plot(ax,x,T.max_abs_standardized,'-s','LineWidth',1.5);
    grid(ax,'on'); box(ax,'on');
    set(ax,'XTick',x,'XTickLabel',viewNames,'XTickLabelRotation',30);
    ylabel(ax,'Absolute standardized value');
    title(ax,'Tail magnitude by horizon');
    legend(ax,{'99th percentile','Maximum'},'Location','best');

    sgtitle(tl,'Preliminary multi-horizon diagnostics','FontWeight','bold');
    save_figure(fig, fullfile(outputDir,'view_diagnostics.png'), dpi);
    close(fig);
end

function plot_signal_heatmaps(X, horizons, commonDates, sectorOrder, ...
        sectorBoundaries, q, outputDir, dpi)

    K = numel(X);
    allAbs = [];
    for k = 1:K
        vals = abs(X{k}(:));
        if numel(vals) > 50000
            idx = round(linspace(1,numel(vals),50000));
            vals = vals(idx);
        end
        allAbs = [allAbs; vals]; %#ok<AGROW>
    end
    cmax = quantile(allAbs,q);
    if ~isfinite(cmax) || cmax <= 0
        cmax = max(allAbs);
    end
    if ~isfinite(cmax) || cmax <= 0
        cmax = 1;
    end

    fig = figure('Visible','off','Color','w', ...
        'Position',[40 40 1900 1050]);
    tl = tiledlayout(ceil(K/2),2,'Padding','compact','TileSpacing','compact');
    cmap = blue_white_red(256);

    for k = 1:K
        ax = nexttile(tl);
        imagesc(ax,X{k}(sectorOrder,:),[-cmax cmax]);
        colormap(ax,cmap);
        hold(ax,'on');
        for b = sectorBoundaries
            yline(ax,b,'k-','LineWidth',0.4);
        end
        hold(ax,'off');
        title(ax,sprintf('%d-day return view',horizons(k)), ...
            'FontWeight','bold');
        ylabel(ax,'Stocks (sector-sorted)');
        xlabel(ax,'Common ending-date coordinate');
        set(ax,'YDir','normal');

        tickPos = unique(round(linspace(1,numel(commonDates),5)));
        tickLabels = cellstr(datestr(commonDates(tickPos),'yyyy-mm-dd'));
        set(ax,'XTick',tickPos,'XTickLabel',tickLabels, ...
            'XTickLabelRotation',30);
    end

    cb = colorbar;
    cb.Layout.Tile = 'east';
    ylabel(cb,'Standardized cumulative log return');
    sgtitle(tl,sprintf(['MOSAIC node-signal matrices | rows sorted by sector | ' ...
        'color limit = %.3g quantile'],q),'FontWeight','bold');

    save_figure(fig, fullfile(outputDir,'signal_heatmaps_by_horizon.png'), dpi);
    close(fig);
end

function C = blue_white_red(m)
    if nargin < 1
        m = 256;
    end
    half = floor(m/2);
    blueToWhite = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];
    redCount = m-half;
    whiteToRed = [ones(redCount,1), linspace(1,0,redCount)', ...
        linspace(1,0,redCount)'];
    C = [blueToWhite; whiteToRed];
end

function write_diagnostic_report(path, cfg, data, T, S, D)
    fid = fopen(path,'w');
    if fid < 0
        warning('Could not create report: %s', path);
        return;
    end
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid,'MOSAIC FINANCE MULTI-HORIZON DATASET REPORT\n');
    fprintf(fid,'============================================================\n\n');
    fprintf(fid,'Nodes                         : %d stocks\n',data.n);
    fprintf(fid,'Views                         : %d return horizons\n',data.K);
    fprintf(fid,'Horizons (trading days)       : %s\n',mat2str(data.horizons));
    fprintf(fid,'Endpoint spacing              : %d trading days\n',data.endpoint_step);
    fprintf(fid,'Signal coordinates per view   : %d\n',data.p);
    fprintf(fid,'Common date range             : %s to %s\n', ...
        datestr(data.common_dates(1),'yyyy-mm-dd'), ...
        datestr(data.common_dates(end),'yyyy-mm-dd'));
    fprintf(fid,'Scaling                       : %s\n',cfg.scaling);
    fprintf(fid,'Winsorization                 : %d\n\n',cfg.winsorize);

    fprintf(fid,'SIGNAL DEFINITION\n');
    fprintf(fid,'-----------------\n');
    fprintf(fid,['For horizon h and common ending price row t_q:\n' ...
        '  R_i^(h)(q) = log(P_i(t_q)) - log(P_i(t_q-h)).\n' ...
        'Each stock is standardized separately within each horizon.\n\n']);

    fprintf(fid,'IMPORTANT INTERPRETATION\n');
    fprintf(fid,'------------------------\n');
    fprintf(fid,['All views use the same ending date at coordinate q. Therefore,\n' ...
        'cross-view distances compare returns of different horizons ending\n' ...
        'on the same calendar dates. No date reordering is used.\n\n']);

    fprintf(fid,'VIEW DIAGNOSTICS\n');
    fprintf(fid,'----------------\n');
    for k = 1:height(T)
        fprintf(fid,['h=%3d | p=%d | mean row std=%.6f | median pair dist=%.4f | ' ...
            'mean lag1=%.4f | within-sector corr=%.4f | between-sector corr=%.4f\n'], ...
            T.horizon_days(k),T.n_coordinates(k),T.mean_row_std(k), ...
            T.median_pairwise_distance(k),T.mean_lag1_autocorr(k), ...
            T.mean_within_sector_corr(k),T.mean_between_sector_corr(k));
    end

    fprintf(fid,'\nAVERAGE SAME-STOCK CROSS-HORIZON CORRELATION\n');
    fprintf(fid,'--------------------------------------------\n');
    write_numeric_matrix(fid,S,data.view_names);

    fprintf(fid,'\nAVERAGE SAME-STOCK CROSS-HORIZON RMS DIFFERENCE\n');
    fprintf(fid,'-----------------------------------------------\n');
    write_numeric_matrix(fid,D,data.view_names);

    fprintf(fid,'\nSCIENTIFIC NOTES\n');
    fprintf(fid,'----------------\n');
    fprintf(fid,['1. The five-day endpoint spacing reduces redundancy but does not\n' ...
        '   make long-horizon return columns independent.\n' ...
        '2. Hyperparameter validation should later use chronological blocks,\n' ...
        '   not a random column split, because cumulative-return windows overlap.\n' ...
        '3. GICS sectors and all preliminary correlation diagnostics are external\n' ...
        '   interpretation tools; they are not supplied to MOSAIC.\n' ...
        '4. The primary signals remain signed. Absolute/squared returns would\n' ...
        '   define a different volatility-network experiment.\n']);
end

function write_numeric_matrix(fid,M,names)
    fprintf(fid,'%14s','');
    for j = 1:numel(names)
        fprintf(fid,'%14s',names{j});
    end
    fprintf(fid,'\n');
    for i = 1:size(M,1)
        fprintf(fid,'%14s',names{i});
        for j = 1:size(M,2)
            fprintf(fid,'%14.4f',M(i,j));
        end
        fprintf(fid,'\n');
    end
end

function save_figure(fig,path,dpi)
    folder = fileparts(path);
    if ~isempty(folder) && ~exist(folder,'dir')
        mkdir(folder);
    end
    try
        exportgraphics(fig,path,'Resolution',dpi);
    catch
        print(fig,path,'-dpng',sprintf('-r%d',dpi));
    end
end
