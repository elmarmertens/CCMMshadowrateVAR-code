%% compare OOS results from two estimates
% load quantico*.mat files and assess OOS performance

clear
close all
fclose all;

%#ok<*DATNM>
%#ok<*DATST>
%#ok<*DATIC>

%% load em toolboxes
warning('off','MATLAB:handle_graphics:exceptions:SceneNode')

path(pathdef)

addpath matlabtoolbox/emtools/
addpath matlabtoolbox/emtexbox/
addpath matlabtoolbox/emgibbsbox/
addpath matlabtoolbox/emeconometrics/
addpath matlabtoolbox/emstatespace/

%% setup

resultsdir = pwd;

doDateTicks        = true;

YLIM = [-8 8];

Nhorizons = 24;

for DATALABEL = {'fredsxMD20-2022-09'}
    
    datalabel = DATALABEL{:};
    
    %% list of models
    % STANDARD
    m = 1;
    models(m).datalabel    = datalabel; %#ok<*SAGROW>
    models(m).resultlabel  = 'standardVARAR1SV-RATSbvarshrinkage-p12';
    models(m).prettylabel  = 'Standard VAR';
    models(m).shortlabel   = 'Standard';
    models(m).fcstType     = 'fcstY';
    
    m = 2;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = 'ELBblocknonstructuralAR1SV-RATSbvarshrinkage-p12';
    models(m).prettylabel  = 'Restricted non-structural shadow-rate VAR';
    models(m).shortlabel   = 'BlockNonstructuralShadowRateVAR';
    models(m).fcstType     = 'fcstY';
    
    
    
    %% parameters
    
    m1  = 1;
    
    m0  = 2;
    shadowshortlabel = 'blocknonstructuralshadowrateAR1SV';
    
    fontsize = 16;
    
    if isdesktop
        wrap = [];
        titlename = sprintf('oosPlotForecasts2023chartsSHADOWFFR-%s-%s', datalabel, shadowshortlabel);
        initwrap
    else
        titlename = sprintf('oosPlotForecasts2023chartsSHADOWFFR-%s-%s', datalabel, shadowshortlabel);
        initwrap
    end
    
    %#ok<*UNRCH>
    %% load data
    
    mat0 = matfile(fullfile(resultsdir, sprintf('%s-%s.mat', models(m0).datalabel, models(m0).resultlabel)));
    mat1 = matfile(fullfile(resultsdir, sprintf('%s-%s.mat', models(m1).datalabel, models(m1).resultlabel)));
    
    ydates    = mat0.ydates;
    
    Tjumpoffs = mat0.Tjumpoffs;
    
    %% cut eval sample if desired
    ndxJumpoff = true(size(Tjumpoffs)); % legacy variable
    
    Tjumpoffs   = Tjumpoffs(ndxJumpoff);
    
    
    dates    = ydates(Tjumpoffs);
    
    %% some parameters
    % Nhorizons  = mat0.fcstNhorizons;
    
    ncode = mat0.ncode;
    
    Ylabels = fredMDprettylabel(ncode);
    N       = length(Ylabels);
    
    Ylabels = strrep(Ylabels, '_', '');
    
    
    %% choose dates to plot
    datechoice = [datenum(2020,3,1), datenum(2020,9,1)];
    
    %% loop over forecast origins
    for thisdate = datechoice
        
        thisT = find(dates == thisdate);
        thisY = find(ismember(ncode, 'FEDFUNDS'));
        
        if doDateTicks
            [thisYear, thisMonth] = datevec(dates(thisT));
            forecastDates = datenum(thisYear, thisMonth + (1 : Nhorizons), 1);
        else
            forecastDates = 1:Nhorizons;
        end
        
        %% realized
        realized = squeeze(mat0.fcstYrealized(thisY,:,thisT)); % actuals from shadow-rate model are censored at ELB
        
        %% load draws
        
        ndxCI = find(ismember(mat0.setQuantiles, normcdf([-1 1]) * 100));
        
        
        fcstYmedian0 = mat0.(sprintf('%smedian', models(m0).fcstType));
        med0 = fcstYmedian0(thisY,1:Nhorizons,thisT);
        fcstYhat0 = mat0.(sprintf('%shat', models(m0).fcstType));
        mid0 = fcstYhat0(thisY,1:Nhorizons,thisT);
        
        fcstYquantiles0 = mat0.(sprintf('%squantiles', models(m0).fcstType));
        tail0 = squeeze(fcstYquantiles0(thisY,1:Nhorizons,ndxCI,thisT));
        
        n = 1;
        
        shadowmid   = mat0.fcstShadowYmedian(n,1:Nhorizons,thisT);
        these       = mat0.fcstShadowYquantiles;
        shadowtails = squeeze(these(n,1:Nhorizons,ndxCI,thisT));
        
        linearmid   = mat1.fcstYmedian(thisY,1:Nhorizons,thisT);
        these       = mat1.fcstYquantiles;
        lineartails = squeeze(these(thisY,1:Nhorizons,ndxCI,thisT));
        
        
        
        %% plot w/linear
        thisfig = figure;
        ax = gca;
        set(ax, 'fontsize', fontsize)
        
        hold on
        
        plotCIaltcolor(shadowmid, shadowtails, forecastDates, [], 'w--', 'linewidth', 3);
        hs = plot(forecastDates, shadowmid, '-',  'color', [0.5843    0.8157    0.9882], 'linewidth', 4);
        plot(forecastDates, shadowmid, 'w--',  'linewidth', 4);
        
        % linear
        hlin = plot(forecastDates, linearmid, 'k-.',  'linewidth', 4);
        plot(forecastDates, lineartails, 'k-.',  'linewidth', 2);
        
        %     hrealized = plot(forecastDates, realized, 'd', 'color', [0 .5 0], 'linewidth', 3);
        
        h1 = plot(forecastDates, med0, 'b-', 'linewidth', 4);
        h1tail = plot(forecastDates, tail0, 'b-', 'linewidth', 2); %#ok<NASGU>
        
        if doDateTicks
            xticks(forecastDates([1, 6 : 6 : end]))
            xlim(forecastDates([1 end]));
            datetick('x', 'yyyy:mm', 'keeplimits', 'keepticks')
        else
            xticks(forecastDates([1, 3 : 3 : end]))
            xlim(forecastDates([1 end]));
        end
        
        ylim(YLIM)
        
        wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity1-%s', ...
            ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
            wrap, [], [], [], [], true)
        ht = title(datestr(dates(thisT), 'yyyy:mm'));
        wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity1-%s-WITHDATE', ...
            ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
            wrap, [], [], [], [], true)
        box off
        hl = legend([hs h1  hlin], 'shadow rate density', 'actual rate (median)', 'linear VAR', 'location', 'northwest'); %#ok<NASGU>
        wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity1-%s-WITHDATELEGEND', ...
            ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
            wrap)
        delete(ht)
        wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity1-%s-WITHLEGEND', ...
            ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
            wrap, [], [], [], [], true)
        
        
        
        %% tabulate
        tabname = sprintf('%s-%s-%s-predictivedensity1lin-%s', ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm'));
        %             call tabulateFigure(wrap, thisPlotLabel, n, ndxT, thisdateT, theseFcstdates, prettylabel1, prettylabel2, prettylabel3, fcstMid1, fcstMid2, fcstMid3, fcstTails1, fcstTails2, fcstTails3, ncode, ndxCI, setQuantiles)
        
        labels4table = {sprintf('%s median', 'SR'), sprintf('%s %4.2f\\%%', 'SR', mat0.setQuantiles(1,ndxCI(1))), sprintf('%s %4.2f\\%%', 'SR', mat0.setQuantiles(1,ndxCI(2))), ...
            sprintf('%s median', 'LIN'), sprintf('%s %4.2f\\%%', 'LIN', mat0.setQuantiles(1,ndxCI(1))), sprintf('%s %4.2f\\%%', 'LIN', mat0.setQuantiles(1,ndxCI(2)))};
        
        data4table = [shadowmid(:) shadowtails linearmid(:) lineartails];
        writedatatable(wrap, tabname, forecastDates, data4table, labels4table, 'yyyy:mm');
        writedatatable2tex(wrap, tabname, forecastDates, data4table, labels4table, 'yyyy:mm');
        
    end
    
    
    %% finish script
    if ~isdesktop
        close all
    end
    finishwrap
    
end % datalabel

dockAllFigures
finishscript
