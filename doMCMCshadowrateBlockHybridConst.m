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
addpath matlabtoolbox/empbsbox/

rng(01012023)

%% Initial operations
clear; close all; clc;


%% set parameters for VAR and MCMC
datalabel           = 'fredblockMD20EBP-2022-09';
jumpDate            = datenum(2022,08,01);
check_stationarity  = 0;                  % Truncate nonstationary draws? (1=yes)

p                   = 12;                    % Number of lags on dependent variables
MCMCdraws           = 1e3;                   % Final number of MCMC draws after burn in
ELBbound            = 0.25;

% SED-PARAMETERS-HERE

%% process ELB setting
switch ELBbound
    case .25
        ELBtag    = '';
    case .125
        ELBtag    = '-ELB125';
    case .5
        ELBtag    = '-ELB500';
    otherwise
        error('ELBbound value of %5.2f not recognized', ELBbound)
end

fcstNhorizons       = [];
fcstNdraws          = [];


doRATSprior         = true;

samStart            = [];                 % truncate start of sample if desired (leave empty if otherwise)

doELBsampling       = true;
doPAIactual         = false;

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


ELBdummy = data(:,ndxSHADOWRATE) <= ELBbound;

startELB       = find(any(ELBdummy,2), 1);
elbT0          = startELB - 1 - p;
% elbT0: first obs prior to missing obs, this is the jump off for the state space
% note: startELB is counted against the available obs in sample, which include
% p additional obs compared to the VAR



%% some parameters
setMinnesotaMean

fontsize = 12;

TID   = parid;

thisT = find(ydates == jumpDate);
T     = thisT - p;

setQuantiles    = [.5, 2.5, 5, normcdf(-1) * 100, 25 , 75,  (1 - normcdf(-1)) * 100, 95, 97.5, 99.5];
Nquantiles      = length(setQuantiles);
ndxCI68         = ismember(setQuantiles, [normcdf(-1) * 100, 100 - normcdf(-1) * 100]);
ndxCI90         = ismember(setQuantiles, [5 95]);
ndxCI           = ndxCI68 | ndxCI90;

%% MCMC sampler

rndStream   = getDefaultStream;

actualrateBlock = ~ismember(1:N, ndxYIELDS);

[PAI_all, PHI_all, invA_all, sqrtht_all, shadowrate_all, ...
    ] = mcmcVARshadowrateBlockHybridConst(thisT, MCMCdraws, ...
    p, np, data, ydates, ...
    actualrateBlock, ...
    minnesotaPriorMean, doRATSprior, ...
    ndxSHADOWRATE, ndxOTHERYIELDS, doELBsampling, false, ELBbound, elbT0, check_stationarity, ...
    [], [], ...
    [], ... % yrealized
    fcstNdraws, fcstNhorizons, rndStream, true);

%% store MCMC
titlename=sprintf('%s%s-p%d-jumpoff%s', datalabel, ELBtag, p, datestr(jumpDate, 'yyyymmm'));

allw = whos;
ndx = contains({allw.class}, 'Figure');
if any(ndx)
    clear(allw(ndx).name)
end
save(sprintf('mcmcBlockHybridConst-%s.mat', titlename), '-v7.3')

%% finish
finishscript
