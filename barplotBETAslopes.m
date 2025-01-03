%% plot BETA coefficients of general model at different jumpoffs

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
%#ok<*NANMEAN>
%#ok<*DATNM>
%#ok<*DATST>

%% setup

doAR1SV = true;

datalabel  = 'fredsxMD20exYield-2022-09';
p = 12;

if doAR1SV
    resultsdir = pwd;
    resultlabel  = sprintf('ELBshadowrateGeneralAR1SV-RATSbvarshrinkage-p%d', p);
    titlename = sprintf('barplotBETA-AR1SV-%s-p%d', datalabel, p);
else

    resultsdir = pwd;
    resultlabel  = sprintf('ELBshadowrateGeneral-RATSbvarshrinkage-p%d', p);
    titlename = sprintf('barplotBETA-%s-p%d', datalabel, p);
end


%% one wrapper per lag choice

wrap = [];
initwrap



%% load data

oos = matfile(fullfile(resultsdir, sprintf('%s-%s.mat', datalabel, resultlabel)));

ydates    = oos.ydates;
Tjumpoffs = oos.Tjumpoffs;
dates     = ydates(Tjumpoffs);

ncode = oos.ncode;
Ylabels = fredMDshortlabel(oos.ncode);
Ylabels = strrep(Ylabels, '_', '');
N       = length(Ylabels);

ndxMacro = ~ismember(1:N,oos.ndxYIELDS);
Nmacro   = sum(ndxMacro);

setQuantiles = oos.setQuantiles;
ndxTails68   = find(ismember(setQuantiles, normcdf11)); % find to enable indexing of matfile objects
ndxTails     = find(ismember(setQuantiles, [5 95])); % find to enable indexing of matfile objects

%% barplot for given jumpoff
for thisT = length(Tjumpoffs) % [ 1 108 length(Tjumpoffs)]
    betaColor  = colors4plots("blue");
    errorColor = colors4plots("lightblue");

    theseBETAtails   = oos.BETAquantiles(:,:,ndxTails,thisT);
    theseBETAtails68 = oos.BETAquantiles(:,:,ndxTails68,thisT);
    mid = oos.BETAmean(:, :, thisT);
    hi  = theseBETAtails(:,:,2) - mid;
    lo  = mid - theseBETAtails(:,:,1);
    if hi < 0
        error houston
    end
    if lo < 0
        error houston
    end
    hi68 = theseBETAtails68(:,:,2) - mid;
    lo68 = mid - theseBETAtails68(:,:,1);
    if lo68 < 0
        error houston
    end
    if hi68 < 0
        error houston
    end
    thisfig = figure;
    bar(1:Nmacro, mid, 'FaceColor', betaColor, 'EdgeColor','flat');
    hold on
    errorbar(1:Nmacro, mid, hi, lo, 'color', errorColor, 'linestyle', 'none', 'linewidth', 2);
    % errorbar(1:Nmacro, mid, hi68, lo68, 'color', errorColor68, 'linestyle', 'none', 'linewidth', 2);
    hold off
    xticklabels(strrep(Ylabels(ndxMacro), '\', ''));
    set(gca, "FontSize", 16)
    set(gca, "box", "off")
    grid on
    % title(sprintf('BETA coefficients at jumpoff %s', datestr(dates(thisT))))
    wrapthisfigure(thisfig, sprintf('%s-%s-%s', 'BETA', resultlabel, datestr(dates(thisT), 'yyyymm')), wrap)
end

%% finish
dockAllFigures
finishwrap
finishscript
