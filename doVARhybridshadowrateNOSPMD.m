%% Hybrid shadow-rate VAR of „Forecasting with Shadow-Rate VARs“
% by Carriero, Clark, Marcellino and Mertens (2021)
% The working paper and supplementary appendices are available here: https://doi.org/10.26509/frbc-wp-202109
%
% Estimation at a given jumpoffDate, and comparison against missing-rate sampling

% This is a simplified version of doVARhybridshadowrate, that works
% without calling spmd and computes only the VAR with default options
% (default optin is shadow-rate sampling instead of missing-rate sampling)

%#ok<*NOSEL>
%#ok<*DISPLAYPROG>
%#ok<*UNRCH>
%#ok<*NASGU>

%% load em toolboxes
warning('off','MATLAB:handle_graphics:exceptions:SceneNode')

path(pathdef)

addpath matlabtoolbox/emtools/
addpath matlabtoolbox/emtexbox/
addpath matlabtoolbox/emgibbsbox/
addpath matlabtoolbox/emeconometrics/
addpath matlabtoolbox/emstatespace/

%% Initial operations
clear; close all; clc;

Nworker     = 1;
rndStreams  = initRandStreams(Nworker, [], 0);


%% set parameters for VAR and MCMC
datalabel           = 'fredMD20baaExYieldBAA-2021-07';
p                   = 12;                  % Number of lags on dependent variables
check_stationarity  = 0;                  % Truncate nonstationary draws? (1=yes)

MCMCdraws           = 1e3;                % Final number of MCMC draws after burn in

jumpoffDate         = datenum(2021,6,1);

% quick options:
MCMCdraws           = 1e2;                % Final number of MCMC draws after burn in
jumpoffDate         = datenum(2017,6,1);



doPlotData          = false;
doVMA               = false;
doTightPrior        = true;

samStart            = [];                 % truncate start of sample if desired (leave empty if otherwise)

ELBbound            = 0.25;

np =12;

%% SED-PARAMETERS-HERE

%% load data
% load CSV file
dum=importdata(sprintf('%s.csv', datalabel),',');


ydates=dum.data(3:end,1);
% Variable names
ncode=dum.textdata(1,2:end);
% Transformation codes (data are already transformed)
tcode  =dum.data(1,2:end);
cumcode=logical(dum.data(2,2:end));
% Data
data=dum.data(3:end,2:end);

setShadowYields

ndxYIELDS = union(ndxSHADOWRATE, ndxOTHERYIELDS);
Nyields   = length(ndxYIELDS);

Nshadowrates = length(ndxSHADOWRATE);
Tdata = length(ydates);

Ylabels = fredMDprettylabel(ncode);

%% process settings
N = size(data,2);
Kbvar = N * p + 1; % number of regressors per equation
K = Kbvar;


% truncate start of sample (if desired)
if ~isempty(samStart)
    ndx  = ydates >= samStart;
    data   = data(ndx,:);
    ydates = ydates(ndx);
    Tdata  = length(ydates);
end



ELBdummy = data(:,ndxSHADOWRATE) <= ELBbound;

startELB       = find(ELBdummy, 1);
elbT0          = startELB - 1 - p;
% elbT0: first obs prior to missing obs, this is the jump off for the state space
% note: startELB is counted against the available obs in sample, which include
% p additional obs compared to the VAR



%% mean for Minnesota prior
setMinnesotaMean

%% Forecast parameters

fcstNhorizons     = 24;  % number of steps forecasted (1:fcstNhorizon)
fcstNdraws        = 10 * MCMCdraws; % draws sampled from predictive density


%% start latexwrapper to collect results
labelELBsampling = 'HybridshadowrateVARvsMissingratesampling';

if doTightPrior
    labelELBsampling = strcat(labelELBsampling, '-tightBVARshrinkage');
end

titlename=sprintf('%s-%s-p%d', datalabel, labelELBsampling, p);
if ~isempty(samStart)
    titlename = strcat(titlename,'-', datestr(samStart, 'yyyymmm'));
end

wrap = [];
initwrap

fontsize = 16;
%% plot input data
if doPlotData
    for n = 1 : N
        this = figure;
        plot(ydates, data(:,n))
        xtickdates(ydates)
        wrapthisfigure(this, sprintf('data%s', ncode{n}), wrap)
    end
end

%% define blocks
actualrateBlock = ~ismember(1:N, ndxSHADOWRATE);

%% loop over QRT estimates

thisT = find(ydates == jumpoffDate);
T     = thisT - p;

%% MCMC sampler


    
    
    
            doELBsampling       = true;
            
    
    doPAIactual = false;
    
    [PAI_all, PHI_all, invA_all, SV_all, shadowrate_all, missingrate_all, ...
        ydraws, yhat, ...
        shadowratedraws, shadowratehat ] ...
        = mcmcVARhybridshadowrate(thisT, MCMCdraws, p, np, data, ydates, ...
        actualrateBlock, ...
        minnesotaPriorMean, doTightPrior, ...
        ndxSHADOWRATE, ndxOTHERYIELDS, doELBsampling, true, ELBbound, elbT0, ...
        check_stationarity, ...
        fcstNdraws, fcstNhorizons, rndStreams{labindex}, true);
    
    
    
    
    % simulate VMA
    
    Kbvar       = 1 + N * p;
    % prepare companion and impulse matrix
    ndxY         = 1+(1:N);
    
    comp0                                = zeros(1 + N * p);
    comp0(1,1)                           = 1; % intercept
    comp0(1 + N + 1 : end,1+(1:N*(p-1))) = eye(N*(p-1));
    
    vma0         = zeros(1+N*p,N);
    vma0(1,:)    = 1;
    vma0(ndxY,:) = eye(N);
    VMAdraws = NaN(N,N,fcstNhorizons,MCMCdraws);
    
    for m = 1 : MCMCdraws
        
        thisPAI      = squeeze(PAI_all(m,:,:));
        
        comp = comp0;
        comp(ndxY,:) = thisPAI(1:Kbvar,:)';
        %     offelbMaxlambdas(m) = max(abs(eig(comp(2:end,2:end))));
        comppow = vma0;
        for h = 1 : fcstNhorizons
            comppow           = comp * comppow;
            VMAdraws(:,:,h,m) = comppow(ndxY,1:N);
        end
    end
    



%% collect PAI and SV composites
% this is legacy code from the version using SPMD composites
PAIshadow  = PAI_all;
SVshadow  = SV_all;

%% plot distribution of sampled shadow rates

% shadowrate moments
shadowrates     = permute(shadowrate_all, [3 2 1]);
shadowrateMid   = median(shadowrates,3);
shadowrateTails = prctile(shadowrates, [5 25 75 95], 3);


%% plot shadowrate
n = 1;

ELBjumpoff = p + elbT0;
thesedates = ydates(ELBjumpoff+1:thisT);
ndxELB     = ELBdummy(ELBjumpoff+1:thisT);
thisfig = figure;
set(gca, 'fontsize', fontsize)
hold on
h0 = plotCI(shadowrateMid(:,n), squeeze(shadowrateTails(:,n,[1 4])), thesedates, [], 'k-', 'linewidth', 3);
plothorzline(ELBbound, [], 'k--')

% ylim([-2.5 2.5])
xtickdates(thesedates)
% ylim([-2.5 2.5])
srlim = ylim;
wrapthisfigure(thisfig, sprintf('hybridshadowrate%dQR-%s', n, datalabel), wrap)



%% wrap up
dockAllFigures
finishwrap
finishscript
