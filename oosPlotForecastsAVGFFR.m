%% plot paths of forecast-relevant rates

clear
close all
fclose all;

%#ok<*DATNM>
%#ok<*DATST>

%% load em toolboxes
warning('off','MATLAB:handle_graphics:exceptions:SceneNode')

path(pathdef)

addpath matlabtoolbox/emtools/
addpath matlabtoolbox/emtexbox/
addpath matlabtoolbox/emgibbsbox/
addpath matlabtoolbox/emeconometrics/
addpath matlabtoolbox/emstatespace/

%% setup

% resultsdir = '~/jam/lager2/QUANTICO/quantico2023logscoresXL';
resultsdir = '../matfilesShadowrateVAR/lagerFREDblock';
doDateTicks = true;
doCumulate  = false;
doOutcomes  = false;

datalabel = 'fredblockMD20-2022-09'; %#ok<*NASGU>
% datalabel = 'fredblockMD20exYield-2022-09'; %#ok<*NASGU>

%% get LSAP dates
[lsapDates, lsapLabels] = getLSAPdates();
lsapDates = datenum(lsapDates);
ndxDropTaper = ~strcmpi(lsapLabels, 'taper begins');
lsapLabels   = lsapLabels(ndxDropTaper);
lsapDates    = lsapDates(ndxDropTaper);

%% list of models
m = 0;
% STANDARD
m = m + 1;
models(m).datalabel    = datalabel; %#ok<*SAGROW>
models(m).resultlabel  = 'standardVAR-RATSbvarshrinkage-p12';
models(m).prettylabel  = 'Standard VAR';
models(m).shortlabel   = 'StandardVAR';
models(m).fcstType     = 'fcstY';
models(m).hasShadow     = false;

% % SHADOW-RATE
% m = m + 1;
% models(m).datalabel    = datalabel; %#ok<*SAGROW>
% models(m).resultlabel  = 'ELBsampling-RATSbvarshrinkage-p12';
% models(m).prettylabel  = 'Simple shadow-rate VAR';
% models(m).shortlabel   = 'ShadowRate';
% models(m).fcstType     = 'fcstY';
% models(m).hasShadow     = true;

% BLOCK-HYBRID
m = m + 1;
models(m).datalabel    = datalabel; %#ok<*SAGROW>
models(m).resultlabel  = 'ELBblockhybrid-RATSbvarshrinkage-p12';
models(m).prettylabel  = 'Hybrid shadow-rate VAR';
models(m).shortlabel   = 'BlockHybridVAR';
models(m).fcstType     = 'fcstY';
models(m).hasShadow     = false;


%% parameters

m0  = 2;

fontsize = 16;

titlename = strcat('oosPlotForecastsFFR-', datalabel);
if doCumulate
    titlename = strcat('CUMULATED-', titlename);
end
if isdesktop
    wrap = [];
else
    initwrap
end

initwrap

%#ok<*UNRCH>
%% load data

varlist = {'ydates', 'Tjumpoffs', 'N', 'tcode', 'ncode', 'cumcode', 'MCMCdraws', ...
    'fcstNhorizons', 'fcstYrealized', 'fcstYhaterror', 'fcstYmederror', 'fcstYcrps', 'fcstYlogscore', ...
    'fcstYhat', 'fcstYmedian', ...
    };

mat0 = matfile(fullfile(resultsdir, sprintf('%s-%s.mat', models(m0).datalabel, models(m0).resultlabel)));

ydates    = mat0.ydates;

Tjumpoffs = mat0.Tjumpoffs;
Njumpoffs = length(Tjumpoffs);

dates     = ydates(Tjumpoffs);

Nhorizons  = mat0.fcstNhorizons;

ncode   = mat0.ncode;
Ylabels = fredMDprettylabel(ncode);
N       = length(Ylabels);
Ylabels = strrep(Ylabels, '_', '');

tcode   = mat0.tcode;


% if isdesktop
%     maturities  = 24;
% else
%     maturities  = [3 6 12 24];
% end

maturities  = 12; % [12 24]; % [3 6 12 24]; % 12;

Nmaturities = length(maturities);

ndxFFR  = find(strcmpi(ncode, 'FEDFUNDS'));
setShadowYields
Nshadowrates = length(ndxSHADOWRATE);

%% collect data
Nmodels = length(models);
avgFFR  = NaN(Njumpoffs, Nmaturities, Nshadowrates, Nmodels);

fcstY     = NaN(Njumpoffs, Nmaturities, N, Nmodels);
realizedY = NaN(Njumpoffs, Nmaturities, N, Nmodels);

for m = 1 : Nmodels
    thismat = load(fullfile(resultsdir, sprintf('%s-%s.mat', models(m).datalabel, models(m).resultlabel)));


    for ndxT = 1 : Njumpoffs

        for hh = 1 : Nmaturities
            thismaturity = maturities(hh);

            % loop over shadowrate variables
            for nn = 1 : Nshadowrates
                if models(m).hasShadow
                    thisYhat = thismat.fcstShadowYhat(nn,1:thismaturity,ndxT);
                else
                    thisYhat = thismat.(sprintf('%shat', models(m).fcstType))(ndxSHADOWRATE(nn),1:thismaturity,ndxT);
                end
                avgFFR(ndxT,hh,nn,m) = mean(thisYhat);
            end

            % loop over all outcomes
            for nn = 1 : N

                if (doCumulate && (tcode(nn) == 5))
                    fcstY(ndxT,hh,nn,m)     = mean(thismat.(sprintf('%shat', models(m).fcstType))(nn,1:thismaturity,ndxT));
                    realizedY(ndxT,hh,nn,m) = mean(thismat.fcstYrealized(nn,1:thismaturity,ndxT));
                else
                    fcstY(ndxT,hh,nn,m)     = thismat.(sprintf('%shat', models(m).fcstType))(nn,thismaturity,ndxT);
                    realizedY(ndxT,hh,nn,m) = thismat.(sprintf('%srealized', models(m).fcstType))(nn,thismaturity,ndxT);
                end

            end


        end


    end
end

%% fix averaging
if doOutcomes
    fcstY(:,:,tcode == 5,:)     = fcstY(:,:,tcode == 5,:) ./ maturities;
    realizedY(:,:,tcode == 5,:) = realizedY(:,:,tcode == 5,:) ./ maturities;
end


%% plot


for doPreCOVID = [true false]

    if doPreCOVID
        sam = dates < datenum(2018,1,1);
        tickdates = datenum(2009:2:2017,1,1);
        sampleLabel = '-preCOVID';
    else
        sam = true(size(dates));
        tickdates = datenum(2010:2:2023,1,1);
        sampleLabel = '';
    end

    h = NaN(length(models),1);
    LineStyleOrder = {':','-.', '-'};
    % newcolors = [0 0.4470 0.7410
    %     1 0 0
    %     0 0 0
    %     0.9290 0.6940 0.1250
    %     0.4940 0.1840 0.5560
    %     0.3010 0.7450 0.9330
    %     0.6350 0.0780 0.1840];
    LineColors = Colors4Plots([8 1 7]);
    LineWidths ={4, 4, 3};

    for nn = 1 % : length(ndxSHADOWRATE)



        for hh = 1 : Nmaturities

            thisfig = figure;
            hold on

            for mm = 1 : length(models)
                h(mm) = plot(dates, squeeze(avgFFR(:,hh,nn,mm)), 'linewidth', LineWidths{mm}, ...
                    'color', LineColors{mm}, 'linestyle', LineStyleOrder{mm});
                if mm == 3
                    uistack(h(mm),'bottom')
                end
            end
            ylim([-5 5])
            
            set(gca, 'fontsize', fontsize)
            xticks(tickdates)
            xtickdates(dates(sam), 'keepticks')
            
            % plotOrigin
            wrapthisfigure(thisfig, sprintf('avg%s-%s-%s-h%d%s', ncode{ndxSHADOWRATE(nn)}, models(m0).datalabel, 'TRIPLE', maturities(hh), sampleLabel), wrap, [], [], [], [], true)

            hl = legend(h, {models.prettylabel}, 'location', 'northwest', 'box', 'on', 'AutoUpdate', 'off');
            wrapthisfigure(thisfig, sprintf('avg%s-%s-%s-h%d%s-WITHLEGEND', ncode{ndxSHADOWRATE(nn)}, models(m0).datalabel, 'TRIPLE', maturities(hh), sampleLabel), wrap, [], [], [], [], false)

            % if doPreCOVID
            %     ylim([-4 4])
            % else
            %     ylim([-5 10])
            % end
            % hLSAP = xline(lsapDates, '--', lsapLabels, 'fontsize', fontsize, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center');
            % set(hl, 'location', 'southeast', 'box','on');
            % wrapthisfigure(thisfig, sprintf('avg%s-%s-%s-h%d%s-LSAP-WITHLEGEND', ncode{ndxSHADOWRATE(nn)}, models(m0).datalabel, 'TRIPLE', maturities(hh), sampleLabel), wrap, [], [], [], [], true)
            % 
            % title(sprintf('%s (h=%d)', Ylabels{ndxSHADOWRATE(nn)}, maturities(hh)))
            % wrapthisfigure(thisfig, sprintf('avg%s-%s-%s-h%d%s-WITHTITLE', ncode{ndxSHADOWRATE(nn)}, models(m0).datalabel, 'TRIPLE', maturities(hh), sampleLabel), wrap, [], [], [], [], false)
        end
    end

    if doOutcomes
        for nn = 1 : N

            for hh = 1 : Nmaturities
                thisfig = figure;

                hold on
                for mm = 1 : length(models)
                    h(mm) = plot(dates, squeeze(fcstY(:,hh,nn,mm)), 'linewidth', 3, 'linestyle', LineStyleOrder{mm});
                end
                colororder(newcolors)
                hr = plot(dates, squeeze(realizedY(:,hh,nn,:)), '-', 'color', [0.4660 0.6740 0.1880], 'linewidth', 1); % note: redundant, should be same line for every model
                set(gca, 'fontsize', fontsize)

                xticks(tickdates)
                xtickdates(dates(sam), 'keepticks')

                wrapthisfigure(thisfig, sprintf('fcst%s-%s-%s-h%d%s', ncode{nn}, models(m0).datalabel, 'TRIPLE', maturities(hh), sampleLabel), wrap, [], [], [], [], true)
                legend([h; hr(1)], {models.prettylabel, 'realized'}, 'location', 'southeast', 'box', 'on')
                wrapthisfigure(thisfig, sprintf('fcst%s-%s-%s-h%d%s-WITHLEGEND', ncode{nn}, models(m0).datalabel, 'TRIPLE', maturities(hh), sampleLabel), wrap, [], [], [], [], true)

                title(sprintf('%s (h=%d)', Ylabels{nn}, maturities(hh)))
                wrapthisfigure(thisfig, sprintf('fcst%s-%s-%s-h%d%s-WITHTITLE', ncode{nn}, models(m0).datalabel, 'TRIPLE', maturities(hh), sampleLabel), wrap, [], [], [], [], false)



            end
        end
    end
end

%% finish script
dockAllFigures
finishwrap
finishscript



