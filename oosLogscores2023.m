%% compare OOS results from two estimates
% load quantico*.mat files and assess OOS performance

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

%#ok<*UNRCH>

%% setup

resultsdir = '~/jam/lager/quantico2023logscoresXL/';

datalabel = 'fredMD20baa-2022-09';

Nmodels     = 4;

for scoreType = {'X', 'YIELDS', 'ALL'}


    %% one big wrapper

    titlename = sprintf('logscores-%s', datalabel);
    switch scoreType{:}
        case 'X'
            titlename = strcat('X', titlename);
            scoreLabel = 'Xlogscores';
        case 'YIELDS'
            titlename = strcat('YIELDS', titlename);
            scoreLabel = 'YIELDSlogscores';
        otherwise
            scoreLabel = 'ALLlogscores';
            titlename = strcat('ALL', titlename);
    end
    wrap = [];
    initwrap



    %% BASELINE p=12

    m = 0;

    % STANDARD
    m = m +  1;
    models(m).datalabel    = datalabel; %#ok<*SAGROW>
    models(m).resultlabel  = 'standardVAR-RATSbvarshrinkage-p12';
    models(m).prettylabel  = 'standard linear VAR';
    models(m).shortlabel   = 'Standard';
    models(m).fcstLogscore = 'fcstYmvlogscoreELB';

    % SHADOW-RATE
    m = m +  1;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = 'ELBsampling-RATSbvarshrinkage-p12';
    models(m).prettylabel  = 'simple shadow-rate VAR';
    models(m).shortlabel   = 'ShadowRateVAR';
    models(m).fcstLogscore = 'fcstYmvlogscore';

    % BLOCK-HYBRID RATE
    m = m +  1;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = 'ELBblockhybrid-RATSbvarshrinkage-p12';
    models(m).prettylabel  = 'hybrid shadow-rate VAR';
    models(m).shortlabel   = 'BlockHybridVAR';
    models(m).fcstLogscore = 'fcstYmvlogscore';

    % % BLOCK-HYBRID w/switching A
    m = m +  1;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = 'ELBblockhybridAelbPrec10000-RATSbvarshrinkage-p12';
    models(m).prettylabel  = 'hybrid shadow-rate VAR (A switching)';
    models(m).shortlabel   = 'BlockHybridVAR-Aelb';
    models(m).fcstLogscore = 'fcstYmvlogscore';

    %% alternative with p=3

    % STANDARD
    m = m +  1;
    models(m).datalabel    = datalabel; %#ok<*SAGROW>
    models(m).resultlabel  = 'standardVAR-RATSbvarshrinkage-p3';
    models(m).prettylabel  = 'standard linear VAR (p=3)';
    models(m).shortlabel   = 'Standard-p3';
    models(m).fcstLogscore = 'fcstYmvlogscoreELB';

    % SHADOW-RATE
    m = m +  1;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = 'ELBsampling-RATSbvarshrinkage-p3';
    models(m).prettylabel  = 'simple shadow-rate VAR  (p=3)';
    models(m).shortlabel   = 'ShadowRateVAR-p3';
    models(m).fcstLogscore = 'fcstYmvlogscore';

    % BLOCK-HYBRID RATE
    m = m +  1;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = 'ELBblockhybrid-RATSbvarshrinkage-p3';
    models(m).prettylabel  = 'hybrid shadow-rate VAR (p=3)';
    models(m).shortlabel   = 'BlockHybridVAR-p3';
    models(m).fcstLogscore = 'fcstYmvlogscore';

    % % BLOCK-HYBRID w/switching A
    m = m +  1;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = 'ELBblockhybridAelbPrec10000-RATSbvarshrinkage-p3';
    models(m).prettylabel  = 'hybrid shadow-rate VAR (A switching, p=3)';
    models(m).shortlabel   = 'BlockHybridVAR-Aelb-p3';
    models(m).fcstLogscore = 'fcstYmvlogscore';


    %% alternative with p=6
    % STANDARD
    m = m +  1;
    models(m).datalabel    = datalabel; %#ok<*SAGROW>
    models(m).resultlabel  = 'standardVAR-RATSbvarshrinkage-p6';
    models(m).prettylabel  = 'standard linear VAR (p=6)';
    models(m).shortlabel   = 'Standard-p6';
    models(m).fcstLogscore = 'fcstYmvlogscoreELB';

    % SHADOW-RATE
    m = m +  1;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = 'ELBsampling-RATSbvarshrinkage-p6';
    models(m).prettylabel  = 'simple shadow-rate VAR  (p=6)';
    models(m).shortlabel   = 'ShadowRateVAR-p6';
    models(m).fcstLogscore = 'fcstYmvlogscore';

    % BLOCK-HYBRID RATE
    m = m +  1;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = 'ELBblockhybrid-RATSbvarshrinkage-p6';
    models(m).prettylabel  = 'hybrid shadow-rate VAR (p=6)';
    models(m).shortlabel   = 'BlockHybridVAR-p6';
    models(m).fcstLogscore = 'fcstYmvlogscore';

    % % BLOCK-HYBRID w/switching A
    m = m +  1;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = 'ELBblockhybridAelbPrec10000-RATSbvarshrinkage-p6';
    models(m).prettylabel  = 'hybrid shadow-rate VAR (A switching, p=6)';
    models(m).shortlabel   = 'BlockHybridVAR-Aelb-p6';
    models(m).fcstLogscore  = 'fcstYmvlogscore';

    %% patch score metrics
    switch scoreType{:}
        case 'X'
            for m = 1 : length(models)
                models(m).fcstLogscore = 'fcstYmvlogscoreX';
            end
        case 'YIELDS'
            for m = 1 : length(models)
                models(m).fcstLogscore = 'fcstYmvlogscoreI';
            end
    end

    %% load benchmark model
    m0        = 1;
    oos0      = matfile(fullfile(resultsdir, sprintf('%s-%s.mat', models(m0).datalabel, models(m0).resultlabel)));
    ydates    = oos0.ydates;
    Tjumpoffs = oos0.Tjumpoffs;
    Njumpoffs = length(Tjumpoffs);

    dates      = ydates(Tjumpoffs);
    Nhorizons  = oos0.fcstNhorizons;

    Ylabels = fredMDprettylabel(oos0.ncode);
    N       = length(Ylabels);

    Ylabels = strrep(Ylabels, '_', '');


    LOGSCORES       = NaN(Njumpoffs, length(models));
    LOGSCORES(:,m0) = oos0.(models(m0).fcstLogscore);

    %% loop over models

    for m = 2 : length(models)


        %% load data
        oos1 = matfile(fullfile(resultsdir, sprintf('%s-%s.mat', models(m).datalabel, models(m).resultlabel)));
        if oos0.MCMCdraws ~= oos1.MCMCdraws
            warning('unequal numbers of MCMCdraws, 0 has %d, 1 has %d', oos0.MCMCdraws, oos1.MCMCdraws)
        end


        if oos0.ydates ~= oos1.ydates
            error('oos estimates based on different samples')
        end
        if ~isequal(oos0.Tjumpoffs, oos1.Tjumpoffs)
            error('oos jumpoffs differ')
        end


        LOGSCORES(:,m) = oos1.(models(m).fcstLogscore);

    end

    %% cumulate scores
    CUMscores  = cumsum(LOGSCORES,1);

    %% define model groups
    ndxExAelb = 1 : length(models);
    ndxExAelb = ndxExAelb(~ismember(ndxExAelb, Nmodels:Nmodels:length(models)));
    GROUPS      = {1:length(models), ndxExAelb, 1:Nmodels, 1:Nmodels-1, Nmodels+(1:Nmodels), Nmodels+(1:Nmodels-1), ...
        2*Nmodels+(1:Nmodels), 2*Nmodels+(1:Nmodels-1)};
    GROUPLABELS = {'all', 'all ex Aelb', 'p12', 'p12 ex Aelb',  'p3', 'p3 ex Aelb', 'p6', 'p6 ex Aelb'};
    GROUPNMODELS = [4 3 4 3 4 3 4 3];

    %% first look
    
    LineTypes = {'-', '--', ':'};
    Colors    = Colors4Plots;
    Colors    = cell2mat(Colors(1:Nmodels));
    preCovid  = dates <= datenum(2020,1,1);

    for gg = 1 : length(GROUPS)

        gndx   = GROUPS{gg};
        grplbl = GROUPLABELS{gg};
        m0     = gndx(1);
        thisNmodel = GROUPNMODELS(gg);

        theseColors = Colors(1:thisNmodel,:); 

        theseSCORES  = CUMscores(:,gndx) - CUMscores(:,m0);
        thisfig = figure;
        h = plot(dates, theseSCORES, 'LineWidth', 3);
        set(gca, 'ColorOrder', theseColors, 'LineStyleOrder', LineTypes);
        %     set(gca, 'FontSize', 12)
        hl = legend(h, {models(gndx).shortlabel}, 'location', 'southoutside');
        hl.NumColumns = 3;
        xtickdates(dates)
        title(sprintf('Models: %s', grplbl))
        wrapthisfigure(thisfig, sprintf('%s-%s-%s', scoreLabel, datalabel, grplbl),wrap);
        
        thisfig = figure;
        h = plot(dates, theseSCORES, 'LineWidth', 3);
        set(gca, 'ColorOrder', theseColors, 'LineStyleOrder', LineTypes);
        %     set(gca, 'FontSize', 12)
        hl = legend(h, {models(gndx).shortlabel}, 'location', 'southoutside');
        hl.NumColumns = 3;
        xtickdates(dates(preCovid))
        title(sprintf('Models: %s -- pre-COVID', grplbl))
        wrapthisfigure(thisfig, sprintf('%s-%s-%s-preCOVID', scoreLabel, datalabel, grplbl),wrap);

    end

    dockAllFigures
    finishwrap

end
%% finish script
finishscript


