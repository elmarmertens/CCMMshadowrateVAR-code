%#ok<*NOSEL>
%#ok<*DISPLAYPROG>
%#ok<*UNRCH>
%#ok<*ASGLU>
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

rng(01012023)

%% Initial operations
clear; close all; clc;

DATALABELS = {'fredMD20VXO-2022-09', 'fredMD20VXOexYield-2022-09', 'fredMD20EBP-2022-09', 'fredMD20EBPexYield-2022-09'};

initwrap
%% set parameters for VAR and MCMC
for dd = 1 : length(DATALABELS)
    datalabel  = DATALABELS{dd};
    jumpDate            = datenum(2022,08,01);
    check_stationarity  = 0;                  % Truncate nonstationary draws? (1=yes)

    % resultsdir           = 'mcmcMatfiles';
    resultsdir          = pwd;

    irfSCALES           = [1 10];
    irfDATES            = [datenum(2007,1,1) datenum(2009,1,1) datenum([2010 2012 2014],12,1)];

    p                   = 12;                    % Number of lags on dependent variables
    irfNdraws           = 1e3;                   % per MCMC node
    irfHorizon          = 25;

    % SED-PARAMETERS-HERE
    resultsdir          = 'irfVXOEBP';

    doRATSprior         = true;

    samStart            = [];                 % truncate start of sample if desired (leave empty if otherwise)

    doELBsampling       = true;
    ELBbound            = 0.25;

    np = 12;

    %% load data
    % load CSV file
    dum=importdata(sprintf('%s.csv', datalabel),',');


    ydates=dum.data(3:end,1);
    % Variable names
    ncode=dum.textdata(1,2:end);
    % Transformation codes (data are already transformed)
    tcode  =dum.data(1,2:end);
    cumcode=logical(dum.data(2,2:end));
    cumcode(tcode == 5) = 1;
    % Data
    data=dum.data(3:end,2:end);

    setShadowYields

    ndxYIELDS = union(ndxSHADOWRATE, ndxOTHERYIELDS);
    Nyields   = length(ndxYIELDS);

    Nshadowrates = length(ndxSHADOWRATE);
    Tdata = length(ydates);

    Ylabels = fredMDprettylabel(ncode);

    %% process settings
    N     = size(data,2);
    K     = N * p + 1; % number of regressors per equation
    Nstates = K + p * Nshadowrates;

    ndxSHADOWRATELAGS = cat(2, false, repmat(ismember(1:N, ndxSHADOWRATE), 1, p)); % prepend by false for CONST

    % truncate start of sample (if desired)
    if ~isempty(samStart)
        ndx    = ydates >= samStart;
        data   = data(ndx,:);
        ydates = ydates(ndx);
        Tdata  = length(ydates);
    end

    % define oos jump offs

    ELBdummy = data(:,ndxSHADOWRATE) <= ELBbound;

    startELB       = find(any(ELBdummy,2), 1);
    elbT0          = startELB - 1 - p;
    % elbT0: first obs prior to missing obs, this is the jump off for the state space
    % note: startELB is counted against the available obs in sample, which include
    % p additional obs compared to the VAR



    %% some parameters

    fontsize = 12;

    TID   = parid;

    thisT = find(ydates == jumpDate);
    T     = thisT - p;

    setQuantiles    = [.5, 2.5, 5, normcdf(-1) * 100, 25 , 75,  (1 - normcdf(-1)) * 100, 95, 97.5, 99.5];
    Nquantiles      = length(setQuantiles);
    ndxCI68         = ismember(setQuantiles, [normcdf(-1) * 100, 100 - normcdf(-1) * 100]);
    ndxCI90         = ismember(setQuantiles, [5 95]);
    ndxCI           = ndxCI68 | ndxCI90;

    rndStream   = getDefaultStream;

    actualrateBlock = ~ismember(1:N, ndxYIELDS);

    %% collect MCMC results
    matfilename=sprintf('%s-p%d-jumpoff%s', datalabel, p, datestr(jumpDate, 'yyyymmm'));
    mcmc = matfile(fullfile(resultsdir, sprintf('mcmcBlockHybrid-%s', matfilename)));
    MCMCdraws = mcmc.MCMCdraws;
    dates     = mcmc.ydates;
    dates     = dates(p+1:end);


    %% plot SV
    nn = 1;
    SV        = mcmc.sqrtht_all;
    SVmedian1 = median(SV(:,:,nn),1);
    SVtails1  = prctile(SV(:,:,nn),setQuantiles, 1);

    thisfig = figure;
    hold on
    plot(dates, SVmedian1, 'k-', 'LineWidth', 2)
    plot(dates, SVtails1(ndxCI68,:), 'k:', 'LineWidth', 1)
    xtickdates(dates)
    title(sprintf('SV of %s shock\n (%s)', Ylabels{nn}, datalabel))
    wrapthisfigure(thisfig, sprintf('SV-%s-%s', Ylabels{nn}, matfilename), wrap)

    exportdata = transpose([SVmedian1; SVtails1(ndxCI68,:)]);
    writedatatable(wrap, sprintf('SV-%s-%s', Ylabels{nn}, matfilename), dates, exportdata, {'median', '15%', '68%'})

end % datalabel
%% finish
dockAllFigures
finishwrap
finishscript

