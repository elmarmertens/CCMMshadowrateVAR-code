%% produce charts of data

%#ok<*NOSEL>
%#ok<*DISPLAYPROG>
%#ok<*UNRCH>
%#ok<*DATNM>

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


datalabel           = 'fredMD20-2022-09';
doQuarterly         = false;

samStart            = []; % datenum(1988,12,1);                 % truncate start of sample if desired (leave empty if otherwise)


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



initwrap
% wrap = [];
fontsize = 12;

%% code beamer frames
fid = fopen('foo.tex', 'wt');


for n = 1 : N
    fprintf(fid, '%s\n', '\begin{frame}');
    fprintf(fid, '%s\n', '\setlength{\skipper}{.25\baselineskip}');
    fprintf(fid, '\\frametitle{%s}\n', upper(Ylabels{n}));
    switch tcode(n)
        case 1
            % % no transformation
        case 2
            fprintf(fid, '%s\n', '\framesubtitle{\fmath{\Delta x_t}}');
        case 4
            fprintf(fid, '%s\n', '\framesubtitle{\fmath{\log(x_t)}}');
        case 5
            fprintf(fid, '%s\n', '\framesubtitle{\fmath{\Delta\log(x_t) \cdot 1200}}');
        otherwise
            error houston
    end
    fprintf(fid, '%s\n', '\vspace{-.75\baselineskip}');
    
    fprintf(fid, '%s\n', '\begin{center}');
    
    fprintf(fid, '\\includegraphics[width=\\textwidth]{data%s}\n', ncode{n});
    
    
    fprintf(fid, '%s\n', '\end{center}');
    
    fprintf(fid, '%s\n\n\n', '\end{frame}');
end

fclose(fid);
type('foo.tex')

%% identify intrest rates
ndxFFR  = strcmp(ncode, 'FEDFUNDS');
ndxGS5  = strcmp(ncode, 'GS5');
ndxGS10 = strcmp(ncode, 'GS10');
ndxTB6M = strcmp(ncode, 'TB6MS');
ndxGS1  = strcmp(ncode, 'GS1');

ndxYIELDS = ndxFFR | ndxTB6M | ndxGS1 | ndxGS5 | ndxGS10;

%% plot termstructure
sam = ydates >= datenum(2008,1,1);

this = figure;
hold on
hFFR = plot(ydates(sam), data(sam,ndxFFR),  'linewidth', 2);
% h6M  = plot(ydates(sam), data(sam,ndxTB6M),  'linewidth', 2);
% h1Y  = plot(ydates(sam), data(sam,ndxGS1),  'linewidth', 2);
h5Y = plot(ydates(sam), data(sam,ndxGS5), 'linewidth', 2);
h10Y = plot(ydates(sam), data(sam,ndxGS10), 'linewidth', 2);
plothorzline(0.25, [], 'k:', 'linewidth', 2)
ax = gca;
set(ax, 'fontsize', fontsize)
box(ax, 'off')
xtickdates(ydates(sam))
wrapthisfigure(this, sprintf('data%s', 'TERMSTRUCTURE'), wrap)
legend([hFFR h5Y h10Y], 'FFR', '5-year', '10-year');
wrapthisfigure(this, sprintf('data%s-WITHLEGEND', 'TERMSTRUCTURE'), wrap)
% title('Interest rates')
% wrapthisfigure(this, sprintf('data%s-WITHLEGENDTITLE', 'TERMSTRUCTURE'), wrap)
h6M  = plot(ydates(sam), data(sam,ndxTB6M),  '--', 'linewidth', 2);
h1Y  = plot(ydates(sam), data(sam,ndxGS1), '-.',  'linewidth', 2);
hl = legend([hFFR h6M h1Y h5Y h10Y], 'FFR', '6-month', '1-year', '5-year', '10-year');
wrapthisfigure(this, sprintf('data%s-WITHLEGEND', 'TERMSTRUCTUREALL'), wrap)
delete(hl)
wrapthisfigure(this, sprintf('data%s', 'TERMSTRUCTUREALL'), wrap)

% this = figure;
% hold on
% hFFR = plot(ydates(sam), data(sam,ndxFFR),  'linewidth', 2);
% % h6M  = plot(ydates(sam), data(sam,ndxTB6M),  'linewidth', 2);
% % h1Y  = plot(ydates(sam), data(sam,ndxGS1),  'linewidth', 2);
% h5Y = plot(ydates(sam), data(sam,ndxGS5), 'linewidth', 2);
% h10Y = plot(ydates(sam), data(sam,ndxGS10), 'linewidth', 2);
% plothorzline(0.25, [], 'k:', 'linewidth', 2)
% ax = gca;
% set(ax, 'fontsize', fontsize)
% box(ax, 'off')
% xtickdates(ydates(sam))
% wrapthisfigure(this, sprintf('data%s', 'TERMSTRUCTURE'), wrap)
% legend([hFFR h5Y h10Y], 'FFR', '5-year', '10-year')
% wrapthisfigure(this, sprintf('data%s-WITHLEGEND', 'TERMSTRUCTURE'), wrap)
% title('Interest rates')
% wrapthisfigure(this, sprintf('data%s-WITHLEGENDTITLE', 'TERMSTRUCTURE'), wrap)

%% plot full sample termstructure
sam = true(size(ydates));

this = figure;
hold on
hFFR = plot(ydates(sam), data(sam,ndxFFR), 'linewidth', 2);
% h6M  = plot(ydates(sam), data(sam,ndxTB6M),  'linewidth', 2);
% h1Y  = plot(ydates(sam), data(sam,ndxGS1),  'linewidth', 2);
h5Y  = plot(ydates(sam), data(sam,ndxGS5),  'linewidth', 2);
h10Y = plot(ydates(sam), data(sam,ndxGS10), 'linewidth', 2);
plothorzline(0.25, [], 'k:', 'linewidth', 2)
ax = gca;
set(ax, 'fontsize', fontsize)
box(ax, 'off')
xtickdates(ydates(sam))
wrapthisfigure(this, sprintf('alldata%s', 'TERMSTRUCTURE'), wrap)
legend([hFFR h5Y h10Y], 'FFR', '5-year', '10-year')
wrapthisfigure(this, sprintf('alldata%s-WITHLEGEND', 'TERMSTRUCTURE'), wrap)
title('Interest rates')
wrapthisfigure(this, sprintf('alldata%s-WITHLEGENDTITLE', 'TERMSTRUCTURE'), wrap)


%% RECENT
sam = ydates >= datenum(2019,1,1); 
this = figure;
hold on
hFFR = plot(ydates(sam), data(sam,ndxFFR),  'linewidth', 2);
h5Y = plot(ydates(sam), data(sam,ndxGS5),  'linewidth', 2);
h10Y = plot(ydates(sam), data(sam,ndxGS10), 'linewidth', 2);
plothorzline(0.25, [], 'k:', 'linewidth', 2)
ax = gca;
set(ax, 'fontsize', fontsize)
box(ax, 'off')
thesedates = ydates(sam);
xlim(thesedates([1 end]))
% xticks([thesedates(1:6:end-1); thesedates(end)])
xticks(thesedates(1:6:end))
datetick('x', 'yyyy:mm', 'keeplimits', 'keepticks')
wrapthisfigure(this, sprintf('recentdata%s', 'TERMSTRUCTURE'), wrap)
legend([hFFR h5Y h10Y], 'FFR', '5-year', '10-year')
wrapthisfigure(this, sprintf('recentdata%s-WITHLEGEND', 'TERMSTRUCTURE'), wrap)
title('Interest rates')
wrapthisfigure(this, sprintf('recentdata%s-WITHLEGENDTITLE', 'TERMSTRUCTURE'), wrap)

%% report averages
sam = ydates >= datenum(2020,4,1) & ydates <= datenum(2020,12,1);
mean(data(sam, ndxYIELDS))

sam = ydates >= datenum(2021,1,1);
mean(data(sam, ndxYIELDS))


%% wrap up
dockAllFigures
finishwrap





