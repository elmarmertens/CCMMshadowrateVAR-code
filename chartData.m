%% produce charts of data

%#ok<*NOSEL>
%#ok<*DISPLAYPROG>
%#ok<*UNRCH>

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


datalabel           = 'fredMD18wuxia-2021-07';
doQuarterly         = false;

samStart = [];

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


Tdata = length(ydates);

Ylabels = fredMDprettylabel(ncode);
Ylabelslong = Ylabels;
Ylabelslong(ismember(tcode, [2 5])) = strcat(Ylabelslong(ismember(tcode, [2 5])), ' (APR growth)');

%% process settings
N = size(data,2);

% truncate start of sample (if desired)
if ~isempty(samStart)
    ndx  = ydates >= samStart;
    data   = data(ndx,:);
    ydates = ydates(ndx);
    Tdata  = length(ydates);
end

yearStart = year(ydates(1));
yearStop  = year(ydates(end)) + 9;

if doQuarterly
    fcstdates = genrQdates(yearStart, yearStop);
else
    fcstdates = genrMdates(yearStart, yearStop, 1);
end
ndx = find(fcstdates == ydates(1));
if isempty(ndx)
    error houston
end
fcstdates = fcstdates(ndx:end);


wrap = [];
% initwrap

fontsize = 12;

%% plot input data
for n = 1 : N
    
    % determne SW outliers (based on pre COVID sample)
    %     sam = ydates < datenum(2020,1,1);
    sam = true(size(ydates)); 
    dev = abs(data(:,n) - median(data(sam,n)));
    iqr = range(prctile(data(sam,n), [25 75]));
    outndx = dev > 5 * iqr;
    
    this = figure;
    hold on
    plot(ydates, data(:,n), 'k-', 'linewidth', 2)
    plot(ydates(outndx), data(outndx,n), 'rd', 'linewidth', 2)
    ax = gca;
    set(ax, 'fontsize', fontsize)
    box(ax, 'off')
    xtickdates(ydates)
    
    wrapthisfigure(this, sprintf('data%s', ncode{n}), wrap)
    title(Ylabelslong{n})
    wrapthisfigure(this, sprintf('data%s-WITHTITLE', ncode{n}), wrap)
end



%% wrap up
dockAllFigures
finishwrap





