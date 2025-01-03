%% non-structural VAR of <<Forecasting with Shadow-Rate VARs>>
% by Carriero, Clark, Marcellino and Mertens (2021)
% The working paper and supplementary appendices are available here: https://doi.org/10.26509/frbc-wp-202109
%
% Recursive estimation of quasi-real-time forecasts

%#ok<*NOSEL>
%#ok<*DISPLAYPROG>
%#ok<*UNRCH>
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

%% Initial operations
clear; close all; clc;

tic % start clocking time

% NOTE: to utilize multiple cores, "parpool" must be started prior to executing the script
%       (or MATLAB must be enabled to launch parpool as needed)
%        otherwise, only a single core gets used (and Nstreams will return 1 after executing the next line)
Nstreams    = max(1,getparpoolsize);
rndStreams  = parallel.pool.Constant(RandStream('Threefry'));

%% set parameters for VAR and MCMC

datalabel           = 'fredsxMD20exYield-2022-09';
doQuarterly         = false;
doRATSprior         = true;
MCMCdraws           = 1e3;               % Final number of MCMC draws after burn in
fcstNdraws          = 10 * MCMCdraws;    % draws sampled from predictive density
doELBsampling       = true;              % shadowrate sampling if true, otherwise treats fedfunds as missing data

ELBbound    = 0.25;

if doQuarterly
    p  = 4;
    np = 4; % number of periods per year, used for calibrating priors
    datalabel = strcat(datalabel, '-quarterly');
else
    p  = 12;
    np = 12;
end

% SED-PARAMETERS-HERE
       


doStoreXL           = false; %#ok<*NASGU>
check_stationarity  = 0;                  % Truncate nonstationary draws? (1=yes)
Compute_diagnostics = false;              % compute Inefficiency Factors and Potential


doLoMem             = true; % do not store memory intensive stuff, just do oos forecasts

doPlotData          = false;

samStart            = [];                 % truncate start of sample if desired (leave empty if otherwise)



if doLoMem
    doStoreXL = false;
end

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
N = size(data,2);

setShadowYields

Nyields   = length(ndxYIELDS);

Nshadowrates = length(ndxSHADOWRATE);
Tdata = length(ydates);

Ylabels = fredMDprettylabel(ncode);

%% process settings
Kbvar = N * p + Nshadowrates * p + 1; % number of regressors per equation
K     = Kbvar;

modellabel = 'ELBnonstructuralAR1SV';


if ELBbound ~= 0.25
    modellabel = strcat(modellabel, sprintf('-ELB%d', ELBbound * 1000));
end

if doRATSprior
    modellabel = strcat(modellabel, '-RATSbvarshrinkage');
end




% truncate start of sample (if desired)
if ~isempty(samStart)
    ndx  = ydates >= samStart;
    data   = data(ndx,:);
    ydates = ydates(ndx);
    Tdata  = length(ydates);
end

% define oos jump offs
Tjumpoffs   = find(ydates > datenum(2008,12,1));

Njumpoffs = length(Tjumpoffs);

ELBdummy = data(:,ndxSHADOWRATE) <= ELBbound;

startELB    = find(any(ELBdummy,2), 1);
elbT0       = startELB - 1 - p;
% elbT0: first obs prior to missing obs, this is the jump off for the state space
% note: startELB is counted against the available obs in sample, which include
% p additional obs compared to the VAR

% other settings

setQuantiles = [.5, 2.5, 5, normcdf(-1) * 100, 25 , 75,  (1 - normcdf(-1)) * 100, 95, 97.5, 99.5];
Nquantiles   = length(setQuantiles);
fractiles = [normcdf(-1) * 100, 100 - normcdf(-1) * 100];
ndxCI68   = ismember(setQuantiles, fractiles);
ndxCI90   = ismember(setQuantiles, [5 95]);
ndxCI     = ndxCI68 | ndxCI90;

%% mean for Minnesota prior
setMinnesotaMean


%% allocate memory for tracking random states
randomStates      = NaN(17, Njumpoffs);

%% allocate QRT memory for shadowrate draws
shadowrateVintagesMid   = NaN(length(ydates), Nshadowrates, Njumpoffs);
shadowrateVintagesTails = NaN(length(ydates), Nshadowrates, 4, Njumpoffs);

missingrateVintagesMid   = NaN(length(ydates), Nshadowrates, Njumpoffs);
missingrateVintagesTails = NaN(length(ydates), Nshadowrates, 4, Njumpoffs);

%% allocate memory for out-of-sample forecasts

fcstNhorizons     = 48;  % number of steps forecasted (1:fcstNhorizon)

% fcstYdraws        = NaN(N,fcstNhorizons,fcstNdraws,Njumpoffs);
fcstYrealized     = NaN(N,fcstNhorizons,Njumpoffs);
% fcstYhatRB        = NaN(N,fcstNhorizons,Njumpoffs); % predictive mean (linear RB)

% linear forecasts
fcstYhat          = NaN(N,fcstNhorizons,Njumpoffs); % predictive mean
fcstYmedian       = NaN(N,fcstNhorizons,Njumpoffs); % predictive median
fcstYhaterror     = NaN(N,fcstNhorizons,Njumpoffs);
fcstYmederror     = NaN(N,fcstNhorizons,Njumpoffs);
fcstYcrps         = NaN(N,fcstNhorizons,Njumpoffs);
fcstYquantiles    = NaN(N,fcstNhorizons,Nquantiles, Njumpoffs);

fcstYmvlogscoreDraws = NaN(fcstNdraws,Njumpoffs); % one-step ahead only
fcstYmvlogscore      = NaN(1,Njumpoffs); % one-step ahead only

fcstYmvlogscoreXdraws = NaN(fcstNdraws,Njumpoffs); % one-step ahead only
fcstYmvlogscoreX      = NaN(1,Njumpoffs); % one-step ahead only
fcstYmvlogscoreIdraws = NaN(fcstNdraws,Njumpoffs); % one-step ahead only
fcstYmvlogscoreI      = NaN(1,Njumpoffs); % one-step ahead only


% cumulated forecasts
fcstYcumrealized     = NaN(N,fcstNhorizons,Njumpoffs);
fcstYcumhat          = NaN(N,fcstNhorizons,Njumpoffs); % predictive mean
fcstYcummedian       = NaN(N,fcstNhorizons,Njumpoffs); % predictive median
fcstYcumhaterror     = NaN(N,fcstNhorizons,Njumpoffs);
fcstYcummederror     = NaN(N,fcstNhorizons,Njumpoffs);
fcstYcumcrps         = NaN(N,fcstNhorizons,Njumpoffs);
fcstYcumquantiles    = NaN(N,fcstNhorizons,Nquantiles, Njumpoffs);



% censored forecasts
% fcstYcensorhat          = NaN(N,fcstNhorizons,Njumpoffs); % predictive mean
% fcstYcensormedian       = NaN(N,fcstNhorizons,Njumpoffs); % predictive median
% fcstYcensorhaterror     = NaN(N,fcstNhorizons,Njumpoffs);
% fcstYcensormederror     = NaN(N,fcstNhorizons,Njumpoffs);
% fcstYcensorlogscore     = NaN(N,fcstNhorizons,Njumpoffs);
% fcstYcensorcrps         = NaN(N,fcstNhorizons,Njumpoffs);
% fcstYcensorquantiles    = NaN(N,fcstNhorizons,Nquantiles, Njumpoffs);

% shadow rate forecasts
fcstShadowYhat        = NaN(Nyields, fcstNhorizons, Njumpoffs);
fcstShadowYmedian     = NaN(Nyields, fcstNhorizons, Njumpoffs);
fcstShadowYquantiles  = NaN(Nyields,fcstNhorizons,Nquantiles, Njumpoffs);



[PAImedian, PAImean, PAIstdev]    = deal(NaN(K, N, Njumpoffs));
PAIquantiles                      = NaN(K, N, Nquantiles, Njumpoffs);
[hRHOmedian, hRHOmean, hRHOstdev] = deal(NaN(N, Njumpoffs));
hRHOquantiles                     = NaN(N, Nquantiles, Njumpoffs);
[hBARmedian, hBARmean, hBARstdev] = deal(NaN(N, Njumpoffs));
hBARquantiles                     = NaN(N, Nquantiles, Njumpoffs);

shadowratePSRF = NaN(Nshadowrates,Njumpoffs);
stackAccept    = NaN(MCMCdraws, Njumpoffs);


%% allocate memory for MCMC output (ex forecast)
if ~doLoMem
    drawsPAI       = NaN(MCMCdraws, K, N, Njumpoffs);
    drawsPHI       = NaN(MCMCdraws, N*(N-1)/2+N, Njumpoffs);
    drawsINVA      = NaN(MCMCdraws, N, N, Njumpoffs);
    drawsSQRTHT    = NaN(MCMCdraws, Tdata, N, Njumpoffs);
    
    %% allocate memory for IRF and sum of FFR coeffs
    VMAmid  = NaN(N,N,fcstNhorizons,Njumpoffs);
    VMAtail = NaN(N,N,fcstNhorizons,Nquantiles,Njumpoffs);
    
    if ~isempty(ndxSHADOWRATE)
        sumFFRmid  = NaN(N,Njumpoffs);
        sumFFRtail = NaN(N,Nquantiles,Njumpoffs);
    end
    
end

%% start latexwrapper to collect results
titlename=sprintf('%s-%s-p%d', datalabel, modellabel, p);
if ~isempty(samStart)
    titlename = strcat(titlename, '-', datestr(samStart, 'yyyymmm'));
end
initwrap
% wrap = [];

%% plot input data
if doPlotData
    for n = 1 : N
        this = figure;
        plot(ydates, data(:,n))
        xtickdates(ydates)
        wrapthisfigure(this, sprintf('data%s', ncode{n}), wrap)
    end
end

%% loop over QRT estimates

% progressbar(0)
parfor ndxT = 1 : Njumpoffs % parfor
    
    
    TID   = parid;
    thisT = Tjumpoffs(ndxT);
    T     = thisT - p;
    
    thisStream           = rndStreams.Value;
    thisStream.Substream = ndxT;
    
    
    fprintf('loop %d, thisT %d, with TID %d\n', ndxT, thisT, TID)
    
    %% collect realized values (without cumulation)
    thisdata = data; % to avoid parfor warning
    yrealized = NaN(N, fcstNhorizons);
    for h = 1 : fcstNhorizons
        if thisT + h <= Tdata
            yrealized(:,h) = thisdata(thisT+h,:)';
        end
    end
    
    % set Funds Rate equal to ELB when at ELB
    % (Note: ELB may be set higher than actual funds rate readings, e.g. 25bp)
    if ~isempty(ELBbound)
        yieldsrealized             = yrealized(ndxSHADOWRATE,:);
        ndx                        = yieldsrealized < ELBbound;
        yieldsrealized(ndx)        = ELBbound;
        yrealized(ndxSHADOWRATE,:) = yieldsrealized;
    end
    
    [PAI_all, hRHO_all, hBAR_all, PHI_all, invA_all, sqrtht_all, shadowrate_all, missingrate_all, ...
                ydraws, yhat, ...
                shadowratedraws, shadowratehat, ...
                logscoredraws, ...
                logscoreXdraws, logscoreIdraws ...
                ] = deal([]); % to avoid parfor warning
    %% MCMC sampler
    
    mcmcOK = false;
    while ~mcmcOK
        try % catch crashes and continue
            % launch mcmc sampler
            [PAI_all, hRHO_all, hBAR_all, PHI_all, invA_all, sqrtht_all, shadowrate_all, missingrate_all, ...
                ydraws, yhat, ...
                shadowratedraws, shadowratehat, ...
                logscoredraws, ...
                logscoreXdraws, logscoreIdraws ...
                ] = mcmcVARshadowrateNonstructuralAR1SV(thisT, MCMCdraws, p, np, data, ydates, ...
                [], ...
                minnesotaPriorMean, doRATSprior, ...
                ndxSHADOWRATE, ndxOTHERYIELDS, doELBsampling, false, ELBbound, elbT0, ...
                check_stationarity, ...
                [], cumcode, ... % IRF1scale
                yrealized, ...
                fcstNdraws, fcstNhorizons, thisStream, false); %#ok<PFBNS>
            mcmcOK = true;
        catch ME
            fprintf('Crash at TID %d, thisT %d\n', TID, thisT)
            fprintf('Error message: %s\n', ME.message)
            continue
        end
    end

    randomStates(:,ndxT) = thisStream.State;
	
    %% Convergence diagnostics
    if Compute_diagnostics
        % display('computing convergence diagnostics..')
        Diagnostics(sqrtht_all,invA_all,PAI_all,PHI_all,N,K,MCMCdraws);
    end
    
    %% collect sampled shadow rates
    
    thisELBdummy = ELBdummy; % to avoid parfor warning about variable slicing
    
    % Convergence Diagnostics for shadowrate Draws
    for s = 1 : Nshadowrates
        shadowratePSRF(s,ndxT) = DiagnosticsShadowrate(shadowrate_all(:,s,thisELBdummy(startELB:thisT,s)),s);
    end
    
    
    
    % shadowrate_all is Ndraws x Nshadowrates x Nobs
    % first: permute into Nobs, Nshadowrates, Ndraws
    shadowrate_all = permute(shadowrate_all, [3 2 1]);
    % now compute moments
    shadowrateMid   = median(shadowrate_all,3);
    shadowrateTails = prctile(shadowrate_all, [5 25 75 95], 3);

    
    %% compute out-of-sample forecasts
    
    % a word on parfor strategy:
    % to make matlab better see the intended use of sliced variabes, use
    % local temp variables and then copy those into the slices at end of
    % loop
    
    
    % cumulated forecasts
    ycumrealized            = yrealized;
    ycumdraws               = ydraws;
    ycumhat                 = yhat;
    ycumrealized(cumcode,:) = cumsum(ycumrealized(cumcode,:),2); % ./ (1:fcstNhorizons);
    ycumdraws(cumcode,:,:)  = cumsum(ycumdraws(cumcode,:,:),2); % ./ (1:fcstNhorizons);
    ycumhat(cumcode,:)      = cumsum(ycumhat(cumcode,:),2); % ./ (1:fcstNhorizons);
    
    % CRPS
    yCRPS = NaN(N,fcstNhorizons);
    for h = 1 : fcstNhorizons
        for n = 1 : N % loop over elements of Y
            yCRPS(n,h) = crpsDraws(yrealized(n,h), ydraws(n,h,:));
        end
    end
    
    ycumCRPS = NaN(N,fcstNhorizons);
    for h = 1 : fcstNhorizons
        for n = 1 : N % loop over elements of Y
            ycumCRPS(n,h) = crpsDraws(ycumrealized(n,h), ycumdraws(n,h,:));
        end
    end
    
    
    
    
    %% collect PAI moments
    PAImedian(:,:,ndxT)       = squeeze(median(PAI_all,1));
    PAImean(:,:,ndxT)         = squeeze(mean(PAI_all,1));
    PAIstdev(:,:,ndxT)        = squeeze(std(PAI_all,1,1));
    PAIquantiles(:,:,:,ndxT)  = permute(prctile(PAI_all,setQuantiles,1), [2 3 1]);

    %% collect RHO moments
    hRHOmedian(:,ndxT)       = median(hRHO_all,1);
    hRHOmean(:,ndxT)         = mean(hRHO_all,1);
    hRHOstdev(:,ndxT)        = std(hRHO_all,1);
    hRHOquantiles(:,:,ndxT)  = transpose(prctile(hRHO_all,setQuantiles,1));
    
    %% collect BAR moments
    hBARmedian(:,ndxT)       = median(hBAR_all,1);
    hBARmean(:,ndxT)         = mean(hBAR_all,1);
    hBARstdev(:,ndxT)        = std(hBAR_all,1);
    hBARquantiles(:,:,ndxT)  = transpose(prctile(hBAR_all,setQuantiles,1));
    %% copy results into sliced variables
    
    fcstYrealized(:,:,ndxT) = yrealized;
    % fcstYhatRB(:,:,ndxT)    = yhatRB;
    
    % predictive likelihood scores
    fcstYmvlogscoreDraws(:,ndxT)  = logscoredraws;
    maxlogscoredraw               = max(logscoredraws);
    fcstYmvlogscore(:,ndxT)       = log(mean(exp(logscoredraws - maxlogscoredraw))) + maxlogscoredraw;
    
    fcstYmvlogscoreXdraws(:,ndxT)   = logscoreXdraws;
    maxlogscoredraw                 = max(logscoreXdraws);
    fcstYmvlogscoreX(:,ndxT)        = log(mean(exp(logscoreXdraws - maxlogscoredraw))) + maxlogscoredraw;
    
    fcstYmvlogscoreIdraws(:,ndxT)   = logscoreIdraws;
    maxlogscoredraw                 = max(logscoreIdraws);
    fcstYmvlogscoreI(:,ndxT)        = log(mean(exp(logscoreIdraws - maxlogscoredraw))) + maxlogscoredraw;
    
    % forecast
    ymed = median(ydraws,3);
    fcstYhat(:,:,ndxT)          = yhat;
    fcstYmedian(:,:,ndxT)       = ymed;
    fcstYhaterror(:,:,ndxT)     = yrealized - yhat;
    fcstYmederror(:,:,ndxT)     = yrealized - ymed;
    fcstYcrps(:,:,ndxT)         = yCRPS;
    fcstYquantiles(:,:,:,ndxT)  = prctile(ydraws, setQuantiles, 3);
    
    % cumulated forecast
    fcstYcumrealized(:,:,ndxT)    = ycumrealized;
    ymed = median(ycumdraws,3);
    fcstYcumhat(:,:,ndxT)          = ycumhat;
    fcstYcummedian(:,:,ndxT)       = ymed;
    fcstYcumhaterror(:,:,ndxT)     = ycumrealized - ycumhat;
    fcstYcummederror(:,:,ndxT)     = ycumrealized - ymed;
    fcstYcumcrps(:,:,ndxT)         = ycumCRPS;
    fcstYcumquantiles(:,:,:,ndxT)  = prctile(ycumdraws, setQuantiles, 3);
    
    % censored forecast
    %     ymed = median(ycensordraws,3);
    %     fcstYcensorhat(:,:,ndxT)          = ycensorhat;
    %     fcstYcensormedian(:,:,ndxT)       = ymed;
    %     fcstYcensorhaterror(:,:,ndxT)     = yrealized - ycensorhat;
    %     fcstYcensormederror(:,:,ndxT)     = yrealized - ymed;
    %     fcstYcensorcrps(:,:,ndxT)         = ycensorCRPS;
    %     fcstYcensorquantiles(:,:,:,ndxT)  = prctile(ycensordraws, setQuantiles, 3);
    
    % shadow rate forecast
    fcstShadowYhat(:,:,ndxT)          = shadowratehat;
    fcstShadowYmedian(:,:,ndxT)       = median(shadowratedraws, 3);
    fcstShadowYquantiles(:,:,:,ndxT)  = prctile(shadowratedraws, setQuantiles, 3);
    
        
    if ~doLoMem
        drawsPAI(:,:,:,ndxT)  = PAI_all;
        drawsPHI(:,:,ndxT)    = PHI_all;
        drawsINVA(:,:,:,ndxT) = invA_all;
        % prepare dummy to make parfor work
        dummy                      = NaN(MCMCdraws, Tdata, N);
        dummy(:, p+1:thisT, :)     = sqrtht_all;
        drawsSQRTHT(:, :, :, ndxT) = dummy;
    end
    
    %% store shadowrate results
    jumpoff = p+elbT0;
    
    % need to work with dummy to get around Matlab's parfor rules
    dummy                               = NaN(length(ydates),Nshadowrates);
    dummy(jumpoff+1:thisT,:)            = shadowrateMid;
    shadowrateVintagesMid(:, :, ndxT)   = dummy;
    
    dummy                               = NaN(length(ydates),Nshadowrates,4);
    dummy(jumpoff+1:thisT,:,:)          = shadowrateTails;
    shadowrateVintagesTails(:,:,:,ndxT) = dummy;
    
   
end

%% plot QRT shadow rate

Nvin = size(shadowrateVintagesMid,3);

firstQRTobs = find(ydates <= datenum(2008,12,1), 1, 'last'); % find(~isnan(shadowrateVintagesMid(:,1,1)),1, 'last'); %  suffcient to check with first yield

%% collect realtime
shadowrateQRTmid   = NaN(length(ydates),Nshadowrates);
shadowrateQRTtails = NaN(length(ydates),Nshadowrates,4);
for v = 1 : Nvin
    ndx = find(~isnan(shadowrateVintagesMid(:,1,v)),1, 'last'); % note: sufficient to check only for first yield
    if ~isnan(shadowrateQRTmid(ndx))
        error houston
    end
    shadowrateQRTmid(ndx,:)     = shadowrateVintagesMid(ndx,:,v);
    shadowrateQRTtails(ndx,:,:) = shadowrateVintagesTails(ndx,:,:,v);
end

%% shadowrate: qrt vs final
for n = 1 : Nshadowrates
    thisfig = figure;
    hold on
    plot(ydates, shadowrateQRTmid(:,n), 'k-', 'linewidth', 2)
    plot(ydates, squeeze(shadowrateQRTtails(:,n,:)), 'k--', 'linewidth', 1)
    
    plot(ydates, shadowrateVintagesMid(:,n,end), 'r-', 'linewidth', 2)
    plot(ydates, squeeze(shadowrateVintagesTails(:,n,:,end)), 'r--', 'linewidth', 1)
    
    %     if n == 1 % should be fedfunds
    %         ylim([-8 3])
    %     end
    xtickdates(ydates([firstQRTobs end]))
    if exist('p', 'var')
        title(sprintf('Shadowrate VAR(%d)', p))
    end
    wrapthisfigure(thisfig, sprintf('shadowrate%dQRTp%d', n, p), wrap)
end

%% collect computer system info

thisArch = computer('arch');
COMPUTERmatlab  = ver;
if ismac
    [~, COMPUTERsys]   = system('sysctl -a | grep machdep.cpu ', '-echo');
    [~, COMPUTERbrand] = system('sysctl -a | grep machdep.cpu | grep brand_string ', '-echo');
    if contains(COMPUTERbrand, 'M1') || contains(COMPUTERbrand, 'M2') || contains(COMPUTERbrand, 'M3') || contains(COMPUTERbrand, 'M4')
        COMPUTERbrand = 'AppleSilicon';
    else
        COMPUTERbrand = 'MacOSIntel';
    end
elseif isunix
    [~, COMPUTERsys] = system('cat /proc/cpuinfo ', '-echo');
    COMPUTERbrand = 'IntelUbuntu';
else % ispc
    COMPUTERsys   = 'Intel(R) Xeon(R) Gold 6248R CPU @ 3.0 GHz';
    COMPUTERbrand = 'WindowsXeon';
end

COMPUTERnstreams  = Nstreams;
COMPUTERtotaltime = toc; % stop clocking time

%% store qrt summary
matfilename = sprintf('%s-%s-p%d', datalabel, modellabel, p);
if ~isempty(samStart)
    matfilename = strcat(matfilename, '-', datestr(samStart, 'yyyymmm'));
end
varlist = {'data', 'ydates', 'p', 'Tjumpoffs', 'N', ...
    'ncode', 'tcode', 'cumcode', ...
    'fcst*', 'fcstNhorizons', ...
    'PAI*', 'hRHO*', 'hBAR*', ...
    'shadowrate*', 'missingrate*', ...
    'ndxSHADOWRATE', 'ndxYIELDS', 'ndxOTHERYIELDS', 'ELBbound', 'ELBdummy',...
    'datalabel', 'modellabel', ...
    'doQuarterly', ...
    'setQuantiles', ...
    'MCMCdraws', ...
    'randomStates', ...
    'COMPUTER*'};
if ~doLoMem
    if ~isempty(ndxSHADOWRATE)
        varlist = cat(2, varlist, 'sumFFR*');
    end
    varlist = cat(2, varlist, 'VMA*');
end

if doStoreXL
    matfilename = sprintf('%s-%s-p%d-draws', datalabel, modellabel, p);
    save(matfilename, varlist{:}, 'draws*', '-v7.3');
end

clear *_all
save(matfilename, varlist{:}, '-v7.3');

%% wrap up
dockAllFigures
finishwrap

