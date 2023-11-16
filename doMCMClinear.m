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


%% set parameters for VAR and MCMC
datalabel           = 'fredblockMD20EBP-2022-09';
jumpDate            = datenum(2022,08,01);
check_stationarity  = 0;                  % Truncate nonstationary draws? (1=yes)

p                   = 12;                    % Number of lags on dependent variables
MCMCdraws           = 1e3;                   % Final number of MCMC draws after burn in

% SED-PARAMETERS-HERE

%% process ELB setting
fcstNhorizons       = [];
fcstNdraws          = [];


doRATSprior         = true;

samStart            = [];                 % truncate start of sample if desired (leave empty if otherwise)


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


% truncate start of sample (if desired)
if ~isempty(samStart)
    ndx    = ydates >= samStart;
    data   = data(ndx,:);
    ydates = ydates(ndx);
    Tdata  = length(ydates);
end



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

ELBbound = .25; % legacy argument, not relevant for GIRF
[PAI_all, PHI_all, invA_all, sqrtht_all, ...
    ] = mcmcVAR(thisT, MCMCdraws,...
    p, np, data, ydates, ...
    minnesotaPriorMean, doRATSprior, ...
    ndxYIELDS, ELBbound, check_stationarity, ...
    [], ... % yrealized
    fcstNdraws, fcstNhorizons, rndStream, true);

%% store MCMC
titlename=sprintf('%s-p%d-jumpoff%s', datalabel, p, datestr(jumpDate, 'yyyymmm'));

allw = whos;
ndx = contains({allw.class}, 'Figure');
if any(ndx)
    clear(allw(ndx).name)
end
save(sprintf('mcmcLinear-%s.mat', titlename), '-v7.3')

%% finish
finishscript


%% define forecast simulation as function
function ydraws = simVAR(N, fcstX0, ndxfcstY, cumcode, np, fcstA, fcstB, invA, nushocks, irfHorizon, irfNdraws)

ydraws      = NaN(N,irfHorizon,irfNdraws);
theseShocks = zeros(N, irfHorizon+1); % padded with zeros for use with ltitr

for nn = 1 : irfNdraws
    theseShocks(:,1:irfHorizon)  = invA * nushocks(:,:,nn);
    xdraws         = ltitr(fcstA, fcstB, theseShocks', fcstX0); % faster forecast simulation using ltitr
    ydraws(:,:,nn) = xdraws(2:end,ndxfcstY)';
end

ydraws(cumcode, :,:)         = cumsum(ydraws(cumcode,:,:), 2) / np;

end % function simVARshadowrate