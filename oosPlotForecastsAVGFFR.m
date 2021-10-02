%% plot paths of forecast-relevant rates

clear
close all
fclose all;

%% load em toolboxes
warning('off','MATLAB:handle_graphics:exceptions:SceneNode')

path(pathdef)

addpath matlabtoolbox/emtools/
addpath matlabtoolbox/emtexbox/
addpath matlabtoolbox/emgibbsbox/
addpath matlabtoolbox/emeconometrics/
addpath matlabtoolbox/emstatespace/

%% setup

resultsdir = '~/jam/lager/quanticoELBmatfiles2021cum';

doDateTicks = true;
doCumulate  = false;
doOutcomes  = false;

datalabel = 'fredMD20baa-2021-07';
% datalabel = 'fredMD20baaExYieldBAA-2021-07';




%% list of models
% STANDARD
m = 1;
models(m).datalabel    = datalabel; %#ok<*SAGROW>
models(m).resultlabel  = 'standardVAR-tightBVARshrinkage-p12';
models(m).prettylabel  = 'Standard VAR';
models(m).shortlabel   = 'StandardVAR';
models(m).fcstType     = 'fcstYcum';
models(m).hasShadow     = false;

% SIMPLE SHADOW-RATE VAR
m = 2;
models(m).datalabel    = datalabel; %#ok<*SAGROW>
models(m).resultlabel  = 'simpleshadowrateVAR-tightBVARshrinkage-p12';
models(m).prettylabel  = 'Simple shadow-rate VAR';
models(m).shortlabel   = 'SimpleShadowRateVAR';
models(m).fcstType     = 'fcstYcum';
models(m).hasShadow     = true;

% HYBRID SHADOW-RATE VAR
m = 3;
models(m).datalabel    = datalabel; %#ok<*SAGROW>
models(m).resultlabel  = 'hybridshadowrateVAR-tightBVARshrinkage-p12';
models(m).prettylabel  = 'Hybrid shadow-rate VAR';
models(m).shortlabel   = 'HybridShadowRateVAR';
models(m).fcstType     = 'fcstYcum';
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

maturities  = 12; % [3 6 12 24]; 

Nmaturities = length(maturities);

ndxFFR  = find(strcmpi(ncode, 'FEDFUNDS'));
setShadowYields
Nshadowrates = length(ndxSHADOWRATE);

%% some parameters
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
                
                fcstY(ndxT,hh,nn,m)     = thismat.(sprintf('%shat', models(m).fcstType))(nn,thismaturity,ndxT);
                realizedY(ndxT,hh,nn,m) = thismat.(sprintf('%srealized', models(m).fcstType))(nn,thismaturity,ndxT);
               
            end
            
            
        end
        
        
    end
end

%% plot

sam = dates < datenum(2018,1,1);

h = NaN(length(models),1);
LineStyleOrder = {'-.','-', ':'};
newcolors = [0 0.4470 0.7410
    1 0 0
    0 0 0
    0.9290 0.6940 0.1250
    0.4940 0.1840 0.5560
    0.3010 0.7450 0.9330
    0.6350 0.0780 0.1840];


for nn = 1 : length(ndxSHADOWRATE)
    
    
    
    for hh = 1 : Nmaturities
        
        thisfig = figure;
        hold on
        
        for mm = 1 : length(models)
            h(mm) = plot(dates, squeeze(avgFFR(:,hh,nn,mm)), 'linewidth', 3, 'linestyle', LineStyleOrder{mm});
        end
        if nn == 1
            ylim([-10 4])
        else
            ylim([-2 2])
        end
        colororder(newcolors)
        set(gca, 'fontsize', fontsize)
        xtickdates(dates)
        
        plotOrigin
        xtickdates(dates(sam))
        wrapthisfigure(thisfig, sprintf('avg%s-%s-%s-h%d-preCOVID', ncode{ndxSHADOWRATE(nn)}, models(m0).datalabel, 'TRIPLE', maturities(hh)), wrap, [], [], [], [], true)
        
        legend(h, {models.prettylabel}, 'location', 'best', 'box', 'on')
        wrapthisfigure(thisfig, sprintf('avg%s-%s-%s-h%d-preCOVID-WITHLEGEND', ncode{ndxSHADOWRATE(nn)}, models(m0).datalabel, 'TRIPLE', maturities(hh)), wrap, [], [], [], [], true)
        
        title(sprintf('%s (h=%d)', Ylabels{ndxSHADOWRATE(nn)}, maturities(hh)))
        wrapthisfigure(thisfig, sprintf('avg%s-%s-%s-h%d-preCOVID-WITHTITLE', ncode{ndxSHADOWRATE(nn)}, models(m0).datalabel, 'TRIPLE', maturities(hh)), wrap, [], [], [], [], false)
        
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
            xtickdates(dates)
            
            xtickdates(dates(sam))
            wrapthisfigure(thisfig, sprintf('fcst%s-%s-%s-h%d-preCOVID', ncode{nn}, models(m0).datalabel, 'TRIPLE', maturities(hh)), wrap, [], [], [], [], true)
            legend([h; hr(1)], {models.prettylabel, 'realized'}, 'location', 'best', 'box', 'on')
            wrapthisfigure(thisfig, sprintf('fcst%s-%s-%s-h%d-preCOVID-WITHLEGEND', ncode{nn}, models(m0).datalabel, 'TRIPLE', maturities(hh)), wrap, [], [], [], [], true)
            
            title(sprintf('%s (h=%d)', Ylabels{nn}, maturities(hh)))
            wrapthisfigure(thisfig, sprintf('fcst%s-%s-%s-h%d-preCOVID-WITHTITLE', ncode{nn}, models(m0).datalabel, 'TRIPLE', maturities(hh)), wrap, [], [], [], [], false)
            
            
            
        end
    end
end

%% finish script
dockAllFigures
finishwrap
finishscript



