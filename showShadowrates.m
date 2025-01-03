%% collect shadow rate series as plotted in Figure 2 of the paper and store as a csv table

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

%% prep
clear; close all; clc;

fontsize = 14;


% labels for Panel(a) exYield
datalabel    = 'fredsxMD20exYield-2022-09';
modeltype1   = 'ELBnonstructuralAR1SV';
prettylabel1 = 'Non-structural shadow-rate VAR';
modeltype2   = 'ELBshadowrateGeneralAR1SV';
prettylabel2 = 'General shadow-rate VAR';

% labels for Panel(b) withYield
datalabel3   = 'fredsxMD20-2022-09';
modeltype3   = 'ELBblocknonstructuralAR1SV';
prettylabel3 = 'Non-structural shadow-rate VAR (w/restrictions)';

YLIM = [-8 6];

ELBbound = .25;
ELBcolor = colors4plots(8);
ELBtag   = '';
ELBlegend = '25 bp';

resultsdir = pwd;

%% get LSAP dates
[lsapDates, lsapLabels] = getLSAPdates();
lsapDates = datenum(lsapDates);
ndxDropTaper = ~strcmpi(lsapLabels, 'taper begins');
lsapLabels   = lsapLabels(ndxDropTaper);
lsapDates    = lsapDates(ndxDropTaper);


tailNDX = [1 4]; % [1 4] for 90% or [2 3] for IQR

%% load stuff
%#ok<*UNRCH>

titlename = sprintf('showShadowrate%s', ELBtag);
initwrap


mat1 = matfile(fullfile(resultsdir, sprintf('%s-%s%s-RATSbvarshrinkage-p12.mat', datalabel, modeltype1, ELBtag)));
mat2 = matfile(fullfile(resultsdir, sprintf('%s-%s%s-RATSbvarshrinkage-p12.mat', datalabel, modeltype2, ELBtag)));
mat3 = matfile(fullfile(resultsdir, sprintf('%s-%s%s-RATSbvarshrinkage-p12.mat', datalabel3, modeltype3, ELBtag)));

shadowshortlabel = 'shadowrate';

ydates      = mat1.ydates;
if ydates ~= mat2.ydates
    error('date mismatch')
end
T           = length(ydates);
p           = mat1.p;

shadowrateVintagesMid1   = mat1.shadowrateVintagesMid;
shadowrateVintagesTails1 = mat1.shadowrateVintagesTails;

shadowrateVintagesMid2   = mat2.shadowrateVintagesMid;
shadowrateVintagesTails2 = mat2.shadowrateVintagesTails;

shadowrateVintagesMid3   = mat3.shadowrateVintagesMid;
shadowrateVintagesTails3 = mat3.shadowrateVintagesTails;

ndxSHADOWRATE = mat1.ndxSHADOWRATE;
if ~isequal(ndxSHADOWRATE, mat2.ndxSHADOWRATE)
    error('ndxSHADOWRATE mismatch')
end

%% patch in actual rate data
FREDmd = importdata(sprintf('%s.csv', datalabel),',');
checkdiff(ydates, FREDmd.data(3:end,1));

data           = FREDmd.data(3:end,2:end);
ActualRateData = data(:,ndxSHADOWRATE);

% patch mat1
ndxActual      = isnan(shadowrateVintagesMid1(:,:,end));
this = shadowrateVintagesMid1(:,:,end);
this(ndxActual) = ActualRateData(ndxActual);
shadowrateVintagesMid1(:,:,end) = this;
for nn = 1 : size(shadowrateVintagesTails1,3)
    this                       = shadowrateVintagesTails1(:,:,nn,end);
    this(ndxActual)            = ActualRateData(ndxActual);
    shadowrateVintagesTails1(:,:,nn,end) = this;
end

% patch mat2
this = shadowrateVintagesMid2(:,:,end);
this(ndxActual) = ActualRateData(ndxActual);
shadowrateVintagesMid2(:,:,end) = this;
for nn = 1 : size(shadowrateVintagesTails2,3)
    this                       = shadowrateVintagesTails2(:,:,nn,end);
    this(ndxActual)            = ActualRateData(ndxActual);
    shadowrateVintagesTails2(:,:,nn,end) = this;
end

% patch mat3
this = shadowrateVintagesMid3(:,:,end);
this(ndxActual) = ActualRateData(ndxActual);
shadowrateVintagesMid3(:,:,end) = this;
for nn = 1 : size(shadowrateVintagesTails3,3)
    this                       = shadowrateVintagesTails3(:,:,nn,end);
    this(ndxActual)            = ActualRateData(ndxActual);
    shadowrateVintagesTails3(:,:,nn,end) = this;
end

%% patch in wuxia
wuxia = importdata('WUXIASHADOWRATE.csv');
wuxiaDates = wuxia.data(:,1);
wuxiaRate  = wuxia.data(:,2);

%% patch in krippner
krippner = importdata('KRIPPNERSHADOWRATE.csv');
krippnerDates = krippner.data(:,1);
krippnerRate  = krippner.data(:,2);

%% plot comparison exYield
n = 1; % FFR

thisfig = figure;
set(gca, 'fontsize', fontsize)

hold on
h1 = plotCI(shadowrateVintagesMid1(:,n,end), squeeze(shadowrateVintagesTails1(:,n,tailNDX,end)), ydates, [], 'k-', 'linewidth', 3);

h2      = plot(ydates, shadowrateVintagesMid2(:,n,end), 'r-.', 'linewidth', 3);
h2tails = plot(ydates, squeeze(shadowrateVintagesTails2(:,n,tailNDX,end)), 'r-.', 'linewidth', 2);

% replot mid
plot(ydates, shadowrateVintagesMid1(:,n,end), 'k-', 'linewidth', 3);

xtickdates([datenum(2006,1,1) ydates(end)])
ylim(YLIM)
hELB = yline(ELBbound, ':', 'color', ELBcolor);
hl = legend([h1(1) h2 hELB], prettylabel1, prettylabel2, ELBlegend, ...
    'location', 'northwest', 'box', 'off', 'AutoUpdate', 'off');
thiswrapname = sprintf('%s%d-p%d-%s-%s-vs-%s%s', shadowshortlabel, n, p, datalabel, modeltype1,modeltype2, ELBtag);
wrapthisfigure(thisfig, thiswrapname, wrap, [], [], [], [], true)
hLSAP = xline(lsapDates, '--', lsapLabels, 'fontsize', fontsize, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center'); %#ok<NASGU>
set(hl, 'location', 'southeast', 'box','on');
wrapthisfigure(thisfig, strcat(thiswrapname, '-LSAP'), wrap)

%% plot comparison withYield against WuXia Krippner
n = 1;
thisfig = figure;
set(gca, 'fontsize', fontsize)

hold on
hfinal = plotCI(shadowrateVintagesMid3(:,n,end), squeeze(shadowrateVintagesTails3(:,n,tailNDX,end)), ydates, [], 'k-', 'linewidth', 3);
xtickdates([datenum(2006,1,1) ydates(end)])
hwx = plot(wuxiaDates, wuxiaRate, '-.', 'color', colors4plots('lightblue'), 'linewidth', 3);
hkrip = plot(krippnerDates, krippnerRate, ':', 'color', colors4plots('darkblue'), 'linewidth', 3);
ylim(YLIM)
hELB = yline(ELBbound, ':', 'color', ELBcolor);
hl = legend([hfinal hwx hkrip hELB], prettylabel3, 'Wu-Xia', 'Krippner', ELBlegend, ...
    'location', 'southeast', 'box', 'off', 'AutoUpdate', 'off');
thiswrapname = sprintf('%s%d-p%d-%s-%s%s-wuxiakrippner', shadowshortlabel, n, p, datalabel3, modeltype3, ELBtag);
wrapthisfigure(thisfig,  thiswrapname, wrap)
hLSAP = xline(lsapDates, '--', lsapLabels, 'fontsize', fontsize, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center');
set(hl, 'location', 'southeast', 'box','on');
wrapthisfigure(thisfig, strcat(thiswrapname, '-LSAP'), wrap)


%% tabulate values
tabname      = 'CCMMshadowrates';
data4table   = [shadowrateVintagesMid1(:,n,end), squeeze(shadowrateVintagesTails1(:,n,tailNDX,end)), ...
    shadowrateVintagesMid2(:,n,end), squeeze(shadowrateVintagesTails2(:,n,tailNDX,end)), ...
    shadowrateVintagesMid3(:,n,end), squeeze(shadowrateVintagesTails3(:,n,tailNDX,end))];
labels4table = {sprintf('Median: %s (ex Yields)', prettylabel1), sprintf('5%%: %s (ex Yields)', prettylabel1), sprintf('95%%: %s (ex Yields)', prettylabel1), ...
    sprintf('Median: %s (ex Yields)', prettylabel2), sprintf('5%%: %s (ex Yields)', prettylabel2), sprintf('95%%: %s (ex Yields)', prettylabel2), ...
    sprintf('Median: %s (w/Yields)', prettylabel3), sprintf('5%%: %s (w/Yields)', prettylabel3), sprintf('95%%: %s (w/Yields)', prettylabel3)};

sam          = ydates >= datenum(2006,1,1);
writedatatable(wrap, tabname, ydates(sam), data4table(sam,:), labels4table, 'yyyy-mm');

%% wrap up
dockAllFigures
finishwrap
