%% plot shadow rate / missing-data rate

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

%% Initial operations
clear; close all; clc;

%% set parameters for VAR and MCMC
datalabel           = 'fredMD20exYield-2022-09';
p                   = 12;                  % Number of lags on dependent variables

jumpoffDate         = datenum(2022,8,1);

doRATSprior        = true;

samStart            = [];                 % truncate start of sample if desired (leave empty if otherwise)

ELBbound            = 0.5;

np = 12;

%% process ELB setting
switch ELBbound
    case .25
        resultsdir = '~/jam/lager/quantico2023logscoresXL/';
        ELBtag    = '';
        ELBlegend = '25 bp';
        switch datalabel
            case 'fredMD20exYield-2022-09';
                YLIM      = [-14 4]; 
            case 'fredMD20-2022-09';
                YLIM      = [-8 4]; 
        end
    case .125
        resultsdir = '~/jam/lager/quantico2023logscoresELB125/';
        ELBtag    = '-ELB125';
        ELBlegend = '12.5 bp';
        YLIM      = [-8 4]; % set common yaxis limits for all plots
    case .5
        resultsdir = '~/jam/lager/quantico2023logscoresELB500/';
        ELBtag    = '-ELB500';
        ELBlegend = '50 bp';
        YLIM      = [-14 4]; % set common yaxis limits for all plots
    otherwise
        error('ELBbound value of %5.2f not recognized', ELBbound)
end


%% load results

modellabel = 'ShadowVsMissingRateSampling';

if doRATSprior
    modellabel = strcat(modellabel, '-RATSbvarshrinkage');
end

titlename=sprintf('%s%s-%s-p%d', datalabel, ELBtag, modellabel, p);
if ~isempty(samStart)
    titlename = strcat(titlename,'-', datestr(samStart, 'yyyymmm'));
end

load(fullfile(resultsdir, sprintf('do-%s.mat', titlename)))

% patch
if ELBbound == .5
    ELBtag    = '-ELB500';
end    

%% start latexwrapper to collect results

initwrap

fontsize = 16;


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
    ylim(YLIM)
    wrapthisfigure(thisfig, sprintf('shadowrate%d-vs-missingrate-%s%s', n, datalabel, ELBtag), wrap)


    % missing-data draw comparison
    thisfig = figure;
    set(gca, 'fontsize', fontsize)

    hold on
    h3 = plot(thesedates, missingrate2Mid(:,n), 'b-.', 'linewidth', 3);
    plot(thesedates, squeeze(missingrate2Tails(:,n,[1 4])), 'b-', 'linewidth', 2);
    plothorzline(ELBbound, [], 'k--')
    xtickdates(thesedates)
    ylim(YLIM)
    wrapthisfigure(thisfig, sprintf('shadowrate%d-vs-missingrate-MissingDataSampling-%s%s-step1', n, datalabel, ELBtag), wrap)
    h1 = plot(thesedates, missingrateMid(:,n), 'r-', 'linewidth', 3);
    plot(thesedates, squeeze(missingrateTails(:,n,[1 4])), 'r-', 'linewidth', 2)
    legend([h1 h3], 'shadow-rate VAR', 'missing-data VAR', 'location', 'best')
    wrapthisfigure(thisfig, sprintf('shadowrate%d-vs-missingrate-MissingDataSampling-%s%s', n, datalabel, ELBtag), wrap)


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

    wrapthisfigure(thisfig, sprintf('simple-%s%s-PAI-%s', datalabel, ELBtag, ncode{n}), wrap)

end

%% scatter of intercepts

thisfig = figure;
set(gca, 'fontsize', 16)
hold on
plot(medianPAIshadow(1,:), medianPAImissing(1,:), 'bx', 'linewidth', 2);
xlabel('Shadow-rate VAR coeffs ')
ylabel('Missing-rate VAR coeffs')
title('all intercepts')

wrapthisfigure(thisfig, sprintf('simple-%s%s-interceptPAI', datalabel, ELBtag), wrap)


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
wrapthisfigure(thisfig, sprintf('simple-%s%s-missing-vs-shadow-allPAI', datalabel, ELBtag), wrap)

legend([hlag, hint, h45], 'Lag coefficients', 'Intercepts', '45 degrees', ...
    'box', 'off', 'location', 'northwest')
wrapthisfigure(thisfig, sprintf('simple-%s%s-missing-vs-shadow-allPAI-WITHLEGEND', datalabel, ELBtag), wrap)

title(sprintf('R^2 = %4.2f', scatter.rsqr))
wrapthisfigure(thisfig, sprintf('simple-%s%s-missing-vs-shadow-allPAI-withR2', datalabel, ELBtag), wrap)

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
wrapthisfigure(thisfig, sprintf('simple-%s%s-missing-vs-shadow-PAIlag', datalabel, ELBtag), wrap)

% legend([hlag, hint, h45], 'Lag coefficients', 'Intercepts', '45 degrees', ...
%     'box', 'off', 'location', 'northwest')
wrapthisfigure(thisfig, sprintf('simple-%s%s-missing-vs-shadow-PAIlag-WITHLEGEND', datalabel, ELBtag), wrap)

title(sprintf('R^2 = %4.2f', scatter.rsqr))
wrapthisfigure(thisfig, sprintf('simple-%s%s-missing-vs-shadow-PAIlag-withR2', datalabel, ELBtag), wrap)


%% standardized PAI changes
% PAI0    = squeeze(mean(PAImissing, 1));
% PAI0se  = squeeze(std(PAImissing, 1, 1));
% PAImean = squeeze(mean(PAIshadow, 1));
% PAIdevs = (PAImean - PAI0) ./ PAI0se;
% 
% maxRED = 2.5;

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
wrapthisfigure(thisfig, sprintf('simple-%s%s-ELBPAIall', datalabel, ELBtag), wrap);





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
    wrapthisfigure(thisfig, sprintf('simple-%s%s-missing-vs-shadow-SVmid-%s', datalabel, ELBtag,  ncode{n}), wrap)
    legend([hShadow hMissing], 'Shadow-rate VAR', 'Missing-rate VAR', 'box', 'on', 'location', 'best')
    wrapthisfigure(thisfig, sprintf('simple-%s%s-missing-vs-shadow-SVmid-%s-WITHLEGEND', datalabel, ELBtag, ncode{n}), wrap)




end



%% wrap up
dockAllFigures
finishwrap
finishscript
