%% load em toolboxes
warning('off','MATLAB:handle_graphics:exceptions:SceneNode')

path(pathdef)

addpath matlabtoolbox/emtools/
addpath matlabtoolbox/emtexbox/
addpath matlabtoolbox/emgibbsbox/
addpath matlabtoolbox/emeconometrics/
addpath matlabtoolbox/emstatespace/

%#ok<*DATNM>
%#ok<*DATST>
%#ok<*UNRCH>

%% prep
clear; close all; clc;

fontsize = 16;
doTitle  = false;

ELBbound = .25;
ELBcolor = Colors4Plots(8);



switch ELBbound
    case .25
        resultsdir = '~/jam/lager/quantico2023/';
        ELBtag    = '';
        ELBlegend = '25 bp';
        YLIM         = [-8 4]; % set common yaxis limits for all plots
        YLIMqrt      = [-15 5]; % set common yaxis limits for all plots
    case .125
        resultsdir = '~/jam/lager/quantico2023logscoresELB125/';
        ELBtag    = '-ELB125';
        ELBlegend = '12.5 bp';
        YLIM      = [-8 4]; % set common yaxis limits for all plots
        YLIMqrt   = YLIM;
    case .5
        resultsdir = '~/jam/lager/quantico2023logscoresELB500/';
        ELBtag    = '-ELB500';
        ELBlegend = '50 bp';
        YLIM      = [-8 4]; % set common yaxis limits for all plots
        YLIMqrt   = YLIM;
    otherwise
        error('ELBbound value of %5.2f not recognized', ELBbound)
end

%% patch in wuxia
wuxia = importdata('WUXIASHADOWRATE.csv');
wuxiaDates = wuxia.data(:,1);
wuxiaRate  = wuxia.data(:,2);

%% patch in krippner
krippner = importdata('KRIPPNERSHADOWRATE.csv');
krippnerDates = krippner.data(:,1);
krippnerRate  = krippner.data(:,2);

%% define datasets
dd = 0;

dd = dd + 1;
DATA(dd).datalabel     = 'fredMD20';
DATA(dd).prettylabel   = 'all';
DATA(dd).linestyle     = '-';
DATA(dd).marker        = 'none';

dd = dd + 1;
DATA(dd).datalabel     = 'fredMD20exYield';
DATA(dd).prettylabel   = 'only FFR';
DATA(dd).linestyle     = '-';
DATA(dd).marker        = 'none';

% dd = dd + 1;
% DATA(dd).datalabel     = 'fredMD15plus6M';
% DATA(dd).prettylabel = 'FFR, 6M';
% DATA(dd).linestyle     = '--';
% DATA(dd).marker        = 'none';
% 
% dd = dd + 1;
% DATA(dd).datalabel     = 'fredMD15plus1Y';
% DATA(dd).prettylabel = 'FFR, 1Y';
% DATA(dd).linestyle     = '--';
% DATA(dd).marker        = 'none';

dd = dd + 1;
DATA(dd).datalabel     = 'fredMD15plusShortYields';
DATA(dd).prettylabel   = 'FFR, 6M, 1Y';
DATA(dd).linestyle     = '--';
DATA(dd).marker        = 'none';

% dd = dd + 1;
% DATA(dd).datalabel     = 'fredMD15plus5Y';
% DATA(dd).prettylabel   = 'FFR, 5Y';
% DATA(dd).linestyle     = ':';
% DATA(dd).marker        = 'none';
% 
% dd = dd + 1;
% DATA(dd).datalabel     = 'fredMD15plus10Y';
% DATA(dd).prettylabel   = 'FFR, 10Y';
% DATA(dd).linestyle     = ':';
% DATA(dd).marker        = 'none';
% 
% dd = dd + 1;
% DATA(dd).datalabel     = 'fredMD15plusLongYields';
% DATA(dd).prettylabel = 'FFR, 5Y, 10Y';
% DATA(dd).linestyle     = ':';
% DATA(dd).marker        = 'none';
% 
% dd = dd + 1;
% DATA(dd).datalabel     = 'fredMD15plusBAA';
% DATA(dd).prettylabel   = 'FFR, BAA';
% DATA(dd).linestyle     = ':';
% DATA(dd).marker        = 'none';

dd = dd + 1;
DATA(dd).datalabel     = 'fredMD15plusLongYieldsBAA';
DATA(dd).prettylabel   = 'FFR, 5Y, 10Y, BAA';
DATA(dd).linestyle     = ':';
DATA(dd).marker        = 'none';


%% allocate memory
T      = 762;
Ntails = 4;

[shadowrateQRTmid, shadowrateFINALmid]     = deal(NaN(T, length(DATA)));
[shadowrateQRTtails, shadowrateFINALtails] = deal(NaN(T, length(DATA), Ntails));

thisShadow = 1;

%% loop over models stuff

initwrap

for MODELTYPE = {'ELBsampling', 'ELBblockhybrid'}

    modeltype = MODELTYPE{:};

    switch modeltype
        case 'ELBsampling'
            modellabel = 'simple shadow-rate VAR';
        case 'ELBblockhybrid'
            modellabel = 'hybrid shadow-rate VAR';
    end

    %% collect shadow rates from different data sets
    for dd =  1 : length(DATA)

        datalabel = strcat(DATA(dd).datalabel, '-2022-09');

        mat = matfile(fullfile(resultsdir, sprintf('%s-%s%s-RATSbvarshrinkage-p12.mat', datalabel, modeltype, ELBtag)));

        shadowshortlabel = 'shadowrate';

        Nvin = size(mat.shadowrateQRTmid,3);

        firstQRTobs = mat.Tjumpoffs(1,1);
        ydates      = mat.ydates;
        if ~isequal(T, length(ydates))
            error('mismatch betwenn T and ydates');
        end
        p           = mat.p;

        shadowrateQRTmid(:,dd)     = mat.shadowrateQRTmid(:,thisShadow);
        %         shadowrateQRTtails      = mat.shadowrateQRTtails;
        shadowrateFINALmid(:,dd)   = mat.shadowrateVintagesMid(:,thisShadow,end);
        %         shadowrateVintagesTails = mat.shadowrateVintagesTails;

        ndxSHADOWRATE = mat.ndxSHADOWRATE;

        %% patch in actual rate data
        FREDmd = importdata(sprintf('%s.csv', datalabel),',');
        checkdiff(ydates, FREDmd.data(3:end,1));

        data           = FREDmd.data(3:end,2:end);
        ActualRateData = data(:,mat.ndxSHADOWRATE(1,thisShadow));

        % patch QRT
        ndxActual                      = isnan(shadowrateQRTmid(:,dd));
        shadowrateQRTmid(ndxActual,dd) = ActualRateData(ndxActual);
        %         for nn = 1 : size(shadowrateQRTtails,3)
        %             this                       = shadowrateQRTtails(:,:,nn);
        %             this(ndxActual)            = ActualRateData(ndxActual);
        %             shadowrateQRTtails(:,:,nn) = this;
        %         end

        % patch FINAL
        ndxActual                        = isnan(shadowrateFINALmid(:,dd));
        shadowrateFINALmid(ndxActual,dd) = ActualRateData(ndxActual);

        %         shadowrateVintagesMid(:,:,end) = this;
        %         for nn = 1 : size(shadowrateVintagesTails,3)
        %             this                       = shadowrateVintagesTails(:,:,nn,end);
        %             this(ndxActual)            = ActualRateData(ndxActual);
        %             shadowrateVintagesTails(:,:,nn,end) = this;
        %         end

    end % done collecting shadow rates


    %% plot QRT rates
    thisfig = figure;
    hold on
    set(gca, 'fontsize', fontsize)
    set(gca,'linestyleorder',{':','-.','--'})


    hshadow      = plot(ydates, shadowrateQRTmid, 'linewidth', 2);
    for dd = 1 : length(DATA)
        hshadow(dd).LineStyle = DATA(dd).linestyle;
        hshadow(dd).Marker    = DATA(dd).marker;
    end

    xtickdates(ydates([firstQRTobs end]))
    ylim(YLIMqrt)
    hELB = yline(ELBbound, '--', 'color', ELBcolor); %#ok<NASGU>

    if doTitle
        title(sprintf('Quasi-real-time estimates from %s', modellabel))
    end
    wrapthisfigure(thisfig, sprintf('%s%d-QRTp%d-groups-%s%s', shadowshortlabel, thisShadow, p, modeltype, ELBtag), wrap)
    hl = legend(hshadow, {DATA.prettylabel}, 'location', 'southoutside');
    hl.NumColumns = 3;
    wrapthisfigure(thisfig, sprintf('%s%d-QRTp%d-groups-%s%s-WITHLEGEND', shadowshortlabel, thisShadow, p, modeltype, ELBtag), wrap)

    %% plot FINAL rates
    thisfig = figure;
    hold on
    set(gca, 'fontsize', fontsize)
    set(gca,'linestyleorder',{':','-.','--'})


    hshadow      = plot(ydates, shadowrateFINALmid, 'linewidth', 2);
    for dd = 1 : length(DATA)
        hshadow(dd).LineStyle = DATA(dd).linestyle;
        hshadow(dd).Marker    = DATA(dd).marker;
    end

    xtickdates(ydates([firstQRTobs end]))
    ylim(YLIM)
    hELB = yline(ELBbound, '--', 'color', ELBcolor);

    if doTitle
        title(sprintf('Full-sample estimates from %s', modellabel))
    end
    wrapthisfigure(thisfig, sprintf('%s%d-FINALp%d-groups-%s%s', shadowshortlabel, thisShadow, p, modeltype, ELBtag), wrap)
    hl = legend(hshadow, {DATA.prettylabel}, 'location', 'southoutside');
    hl.NumColumns = 3;
    wrapthisfigure(thisfig, sprintf('%s%d-FINALp%d-groups-%s%s-WITHLEGEND', shadowshortlabel, thisShadow, p, modeltype, ELBtag), wrap)



end

%% wrap up
dockAllFigures
finishwrap
finishscript
