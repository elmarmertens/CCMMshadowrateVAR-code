%% estimate shadow rate / missing-data rate for restricted non-structrual model
% restricted non-structrual model is a.k.a. block-nonstructural model

%#ok<*DATNM>
%#ok<*DATST>
%#ok<*NOSEL>
%#ok<*DISPLAYPROG>
%#ok<*UNRCH>
%#ok<*DLABINDEX>

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

Nworker     = 2;

if getparpoolsize < Nworker
    parpool(Nworker)
end

rndStreams  = parallel.pool.Constant(RandStream('Threefry'));

%% set parameters for VAR and MCMC
datalabel           = 'fredsxMD20-2022-09';
p                   = 12;                  % Number of lags on dependent variables
check_stationarity  = 0;                  % Truncate nonstationary draws? (1=yes)

MCMCdraws           = 1e3;                % Final number of MCMC draws after burn in

jumpoffDate         = datenum(2022,8,1);

doPlotData          = false;
doVMA               = false;
doRATSprior        = true;

samStart            = [];                 % truncate start of sample if desired (leave empty if otherwise)

ELBbound            = 0.25;

np = 12;

setQuantiles = [.5, 2.5, 5, normcdf(-1) * 100, 25 , 75,  (1 - normcdf(-1)) * 100, 95, 97.5, 99.5];
Nquantiles    = length(setQuantiles);
fractiles = [normcdf(-1) * 100, 100 - normcdf(-1) * 100];
ndxCI68   = ismember(setQuantiles, fractiles);
ndxCI90   = ismember(setQuantiles, [5 95]);
ndxCI     = ndxCI68 | ndxCI90;

%% SED-PARAMETERS-HERE
        datalabel='fredsxMD20exYield-2022-09'; 
        p=12; 
        ELBbound=0.25; 
		MCMCdraws=1e3; 
		fcstNdraws= 10 * MCMCdraws; 
        

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

startELB       = find(any(ELBdummy,2), 1);
elbT0          = startELB - 1 - p;
% elbT0: first obs prior to missing obs, this is the jump off for the state space
% note: startELB is counted against the available obs in sample, which include
% p additional obs compared to the VAR



%% mean for Minnesota prior
setMinnesotaMean

%% Forecast parameters (obsolete)

fcstNhorizons     = 24;  % number of steps forecasted (1:fcstNhorizon)
fcstNdraws        = 10 * MCMCdraws; % draws sampled from predictive density


%% start latexwrapper to collect results
modellabel = 'blocknonstructuralShadowVsMissingRateSampling';

if doRATSprior
    modellabel = strcat(modellabel, '-RATSbvarshrinkage');
end

titlename=sprintf('%s%s-%s-p%d', datalabel, ELBtag, modellabel, p);
if ~isempty(samStart)
    titlename = strcat(titlename,'-', datestr(samStart, 'yyyymmm'));
end

wrap = [];
initwrap

fontsize = 16;
%% define blocks
actualrateBlock = ~ismember(1:N, ndxSHADOWRATE);

thisT = find(ydates == jumpoffDate);
T     = thisT - p;

%% MCMC sampler

shadowrate_all  = Composite(Nworker); %#ok<NASGU>
missingrate_all = Composite(Nworker); %#ok<NASGU>
% VMAdraws        = Composite(Nworker); %#ok<NASGU>
SV_all          = Composite(Nworker); %#ok<NASGU>

spmd(Nworker)



    switch labindex


        case 1 % shadow-rate sampling

            doELBsampling       = true;

        case 2 % missing-rate sampling VAR

            doELBsampling       = false;
    end

    doPAIactual = false;

    thisStream           = rndStreams.Value;
    thisStream.Substream = labindex;

    [PAI_all, ~, ~, PHI_all, invA_all, SV_all, shadowrate_all, missingrate_all ...
        ] = mcmcVARshadowrateBlockNonstructuralAR1SV(thisT, MCMCdraws, p, np, data, ydates, ...
        actualrateBlock, ...
        minnesotaPriorMean, doRATSprior, ...
        ndxSHADOWRATE, ndxOTHERYIELDS, doELBsampling, true, ELBbound, elbT0, check_stationarity, ...
        [], [], [], ...
        fcstNdraws, fcstNhorizons, thisStream, true);






end



%% collect PAI and SV composites
PAIshadow  = PAI_all{1};
PAImissing = PAI_all{2};

SVshadow  = SV_all{1};
SVmissing = SV_all{2};

%% plot distribution of sampled shadow rates

% shadowrate moments
shadowrates     = permute(shadowrate_all{1}, [3 2 1]);
shadowrateMid   = median(shadowrates,3);
shadowrateTails = prctile(shadowrates, [5 25 75 95], 3);

% missingrate moments
missingrates     = permute(missingrate_all{1}, [3 2 1]);
missingrateMid   = median(missingrates,3);
missingrateTails = prctile(missingrates, [5 25 75 95], 3);

missingrates2     = permute(missingrate_all{2}, [3 2 1]);
missingrate2Mid   = median(missingrates2,3);
missingrate2Tails = prctile(missingrates2, [5 25 75 95], 3);


%% plot shadowrate
for n = 1 : Nshadowrates

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
    h3 = plot(thesedates, missingrateMid(:,n), 'r-.', 'linewidth', 3);
    h3b = plot(thesedates, squeeze(missingrateTails(:,n,[1 4])), 'r-', 'linewidth', 2);
    legend([h0 h3], 'shadow-rate sampling', 'missing-data sampling', 'location', 'best')
    % ylim([-8 4])
    srlim = ylim;
    wrapthisfigure(thisfig, sprintf('blocknonstructuralshadowrate%d-vs-missingrate-%s%s', n, datalabel, ELBtag), wrap)


    % missing-data draw comparison
    thisfig = figure;
    set(gca, 'fontsize', fontsize)

    hold on
    h3 = plot(thesedates, missingrate2Mid(:,n), 'b-.', 'linewidth', 3);
    plot(thesedates, squeeze(missingrate2Tails(:,n,[1 4])), 'b-', 'linewidth', 2);
    plothorzline(ELBbound, [], 'k--')
    xtickdates(thesedates)
    ylim(srlim)
    wrapthisfigure(thisfig, sprintf('blocknonstructuralshadowrate%d-vs-missingrate-MissingDataSampling-%s%s-step1', n, datalabel, ELBtag), wrap)
    h1 = plot(thesedates, missingrateMid(:,n), 'r-', 'linewidth', 3);
    plot(thesedates, squeeze(missingrateTails(:,n,[1 4])), 'r-', 'linewidth', 2)
    legend([h1 h3], 'shadow-rate VAR', 'missing-data VAR', 'location', 'best')
    wrapthisfigure(thisfig, sprintf('blocknonstructuralshadowrate%d-vs-missingrate-MissingDataSampling-%s%s', n, datalabel, ELBtag), wrap)


end

%% 2-panel version
medianPAIshadow  = squeeze(median(PAIshadow, 1));
medianPAImissing = squeeze(median(PAImissing, 1));


ndxLag1     = 1 + (1:N);
ndxLagOther = N + 2 : K;

for n = 1  :  N

    thisfig = figure;
    subplot(1,2,1)
    set(gca, 'fontsize', 16)
    hold on
    hlag = plot(medianPAIshadow(2:end,n), medianPAImissing(2:end,n), 'rx', 'linewidth', 2);
    plot(medianPAIshadow(1,n), medianPAImissing(1,n), 'bx', 'linewidth', 2)
    hint = plot(medianPAIshadow(1,n), medianPAImissing(1,n), 'bo', 'linewidth', 2);

    legend([hlag hint], 'Lag coefficients', 'Intercept', 'location', 'best')
    xlabel('Shadow-rate VAR coeffs ')
    ylabel('Missing-rate VAR coeffs')
    title('all coefficients')

    subplot(1,2,2)
    set(gca, 'fontsize', 16)
    hold on
    hlagX = plot(medianPAIshadow(ndxLagOther,n), medianPAImissing(ndxLagOther,n), 'kx', 'linewidth', 2);
    hlag1 = plot(medianPAIshadow(ndxLag1,n), medianPAImissing(ndxLag1,n), 'rx', 'linewidth', 2);
    %     plot(medianPAIawayelb(1,n), medianPAIatelb(1,n), 'bx', 'linewidth', 2)
    %     hint = plot(medianPAIawayelb(1,n), medianPAIatelb(1,n), 'bo', 'linewidth', 2);
    legend([hlag1 hlagX], 'Lag 1 coefficients', 'other lags', 'location', 'best')
    xlabel('Shadow-rate VAR coeffs ')
    ylabel('Missing-rate VAR coeffs')
    title('all (w/o intercept)')




    %     sgtitle(sprintf('%s equation (medians)', Ylabels{n}))


    wrapthisfigure(thisfig, sprintf('blocknonstructural-%s%s-PAI-%s', datalabel, ELBtag, ncode{n}), wrap)

end

%% scatter of intercepts

thisfig = figure;
set(gca, 'fontsize', 16)
hold on
plot(medianPAIshadow(1,:), medianPAImissing(1,:), 'bx', 'linewidth', 2);
xlabel('Shadow-rate VAR coeffs ')
ylabel('Missing-rate VAR coeffs')
title('all intercepts')

wrapthisfigure(thisfig, sprintf('blocknonstructural-%s%s-interceptPAI', datalabel, ELBtag), wrap)


%% scatter of all

scatter = ols(medianPAIshadow(:), [ones(length(medianPAIshadow(:)), 1) medianPAImissing(:)]);
% prt(scatter)

thisfig = figure;
set(gca, 'fontsize', 16)
hold on
hlag = plot(vec(medianPAIshadow(2:end,:)), vec(medianPAImissing(2:end,:)), 'kx', 'linewidth', 2);
hint = plot(medianPAIshadow(1,:), medianPAImissing(1,:), 'rx', 'linewidth', 2);
xlabel('Shadow-rate VAR')
ylabel('Missing-rate VAR')


h45 = plot(xlim, ylim, 'k--');
wrapthisfigure(thisfig, sprintf('blocknonstructural-%s%s-missing-vs-shadow-allPAI', datalabel, ELBtag), wrap)

legend([hlag, hint, h45], 'Lag coefficients', 'Intercepts', '45 degrees', ...
    'box', 'off', 'location', 'northwest')
wrapthisfigure(thisfig, sprintf('blocknonstructural-%s%s-missing-vs-shadow-allPAI-WITHLEGEND', datalabel, ELBtag), wrap)

title(sprintf('R^2 = %4.2f', scatter.rsqr))
wrapthisfigure(thisfig, sprintf('blocknonstructural-%s%s-missing-vs-shadow-allPAI-withR2', datalabel, ELBtag), wrap)

stdabsresid = abs(scatter.resid); % / sqrt(scatter.sige);
[~, sortndx] = sort(stdabsresid);

cutoff = stdabsresid(sortndx) > 2;

[PAIrow, PAIcol] = ind2sub(size(medianPAIshadow), sortndx(cutoff));

unique(Ylabels(PAIcol))

%% lag scatter only
thisfig = figure;
set(gca, 'fontsize', 16)
hold on
hlag = plot(vec(medianPAIshadow(2:end,:)), vec(medianPAImissing(2:end,:)), 'kx', 'linewidth', 2);
xlabel('Shadow-rate VAR')
ylabel('Missing-rate VAR')


h45 = plot(xlim, ylim, 'k--');
wrapthisfigure(thisfig, sprintf('blocknonstructural-%s%s-missing-vs-shadow-PAIlag', datalabel, ELBtag), wrap)

% legend([hlag, hint, h45], 'Lag coefficients', 'Intercepts', '45 degrees', ...
%     'box', 'off', 'location', 'northwest')
wrapthisfigure(thisfig, sprintf('blocknonstructural-%s%s-missing-vs-shadow-PAIlag-WITHLEGEND', datalabel, ELBtag), wrap)

title(sprintf('R^2 = %4.2f', scatter.rsqr))
wrapthisfigure(thisfig, sprintf('blocknonstructural-%s%s-missing-vs-shadow-PAIlag-withR2', datalabel, ELBtag), wrap)


%% standardized PAI changes
PAI0    = squeeze(mean(PAImissing, 1));
PAI0se  = squeeze(std(PAImissing, 1, 1));
PAImean = squeeze(mean(PAIshadow, 1));
PAIdevs = (PAImean - PAI0) ./ PAI0se;

maxRED = 2.5;

% plot all
thisfig = figure;

hb = surf(1:N,1:K,abs(PAIdevs));

xlim([1 N])
xticks(1:N)
xticklabels(Ylabels);

yticks([1, 1 + 3 * N])
yticklabels({'intercepts', 'lags'});


% caxis([0 maxRED])


% color by height
shading interp
colorbar
colormap turbo
set(gca, 'fontsize', 12)
view(-20,67)
wrapthisfigure(thisfig, sprintf('blocknonstructural-%s%s-ELBPAIall', datalabel, ELBtag), wrap);






%% compare SV

shadowSVmid   = squeeze(median(SVshadow,1));
shadowSVtails = squeeze(prctile(SVshadow, setQuantiles,1));
missingSVmid = squeeze(median(SVmissing,1));
missingSVtails = squeeze(prctile(SVmissing, setQuantiles,1));

thesedates = ydates(ELBjumpoff:thisT);

for n = 1 : N

    thisfig = figure;
    hold on
    set(gca, 'fontsize', 16)
    hShadow  = plot(thesedates, shadowSVmid(ELBjumpoff-p:end,n), 'r-', 'linewidth', 2);
    hMissing = plot(thesedates, missingSVmid(ELBjumpoff-p:end,n), 'k--', 'linewidth', 2);

    xtickdates(thesedates)
    wrapthisfigure(thisfig, sprintf('blocknonstructural-%s%s-missing-vs-shadow-SVmid-%s', datalabel, ELBtag,  ncode{n}), wrap)
    legend([hShadow hMissing], 'Shadow-rate VAR', 'Missing-rate VAR', 'box', 'on', 'location', 'best')
    wrapthisfigure(thisfig, sprintf('blocknonstructural-%s%s-missing-vs-shadow-SVmid-%s-WITHLEGEND', datalabel, ELBtag, ncode{n}), wrap)

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

COMPUTERnworker   = Nworker;
COMPUTERnstreams  = Nworker; % for sake of compatibility with goVAR* files
COMPUTERtotaltime = toc; % stop clocking time

%% store output
% clean all figures
allw = whos;
ndx = contains({allw.class}, 'Figure') | contains({allw.class}, 'Composite') ...
    | contains({allw.class}, 'graphics');
% clear(allw(ndx).name)
% ndx = contains({allw.class}, 'Composite');
% clear(allw(ndx).name)
save(sprintf('do-%s.mat', titlename), allw(~ndx).name)

%% wrap up
dockAllFigures
finishwrap
finishscript
