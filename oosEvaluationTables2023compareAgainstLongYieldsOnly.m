%% compare OOS results from two estimates

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

resultsdir = '../matfilesShadowrateVAR/lagerFREDblock';

doCharts = false;
doBold   = true; % matter only if doCharts = false;

datalabel0   = 'fredMD14longyields';
prettylabel0 = 'Linear VAR without short-term yields';

% datalabel0   = 'fredblockMD20exYield';
% prettylabel0 = 'all variables ex yields';

%% define alternative datasets
dd = 0;

dd = dd + 1;
ALTDATA(dd).dataset     = 'fredblockMD20';
ALTDATA(dd).prettylabel = 'all variables';

dd = dd + 1;
ALTDATA(dd).dataset     = 'fredblockMD20exYield';
ALTDATA(dd).prettylabel = 'all variables (ex yields)';

% dd = dd + 1;
% ALTDATA(dd).dataset     = 'fredMD15plus6M';
% ALTDATA(dd).prettylabel = 'FFR plus 6M (and 14 others)';
% 
% dd = dd + 1;
% ALTDATA(dd).dataset     = 'fredMD15plus1Y';
% ALTDATA(dd).prettylabel = 'FFR plus 1Y (and 14 others)';
% 
% dd = dd + 1;
% ALTDATA(dd).dataset     = 'fredMD15plus5Y';
% ALTDATA(dd).prettylabel = 'FFR plus 5Y (and 14 others)';
% 
% dd = dd + 1;
% ALTDATA(dd).dataset     = 'fredMD15plus10Y';
% ALTDATA(dd).prettylabel = 'FFR plus 10Y (and 14 others)';
% 
% dd = dd + 1;
% ALTDATA(dd).dataset     = 'fredMD15plusShortYields';
% ALTDATA(dd).prettylabel = 'FFR plus 6M and 1Y (and 14 others)';
% 
% dd = dd + 1;
% ALTDATA(dd).dataset     = 'fredMD15plusLongYields';
% ALTDATA(dd).prettylabel = 'FFR plus 5Y and 10Y (and 14 others)';
% 
% dd = dd + 1;
% ALTDATA(dd).dataset     = 'fredMD15plusBAA';
% ALTDATA(dd).prettylabel = 'FFR plus BAA (and 14 others)';
% 
% dd = dd + 1;
% ALTDATA(dd).dataset     = 'fredMD15plusLongYieldsBAA';
% ALTDATA(dd).prettylabel = 'FFR plus 5Y, 10Y, and BAA (and 14 others)';

%% define samples
s = 0;

s = s + 1;
samples(s).evalStart  = datenum(2009,1,1);
samples(s).evalStop   = datenum(2022,8,1);

s = s + 1;
samples(s).evalStart = datenum(2009,1,1);
samples(s).evalStop  = datenum(2017,12,1);

% s = s + 1;
% samples(s).evalStart = datenum(2018,1,1);
% samples(s).evalStop  = datenum(2022,8,1);
% 
% s = s + 1;
% samples(s).evalStart  = datenum(2009,1,1);
% samples(s).evalStop   = datenum(2015,12,1);

%% define models

m = 0;
% BASELINE p=12
% % STANDARD
% m = m +  1;
% models(m).resultlabel  = 'standardVAR-RATSbvarshrinkage-p12';
% models(m).prettylabel  = 'standard linear VAR';
% models(m).shortlabel   = 'Standard-p12';
% models(m).fcstType     = 'fcstY';

% % SHADOW-RATE
% m = m + 1;
% models(m).resultlabel  = 'ELBsampling-RATSbvarshrinkage-p12';
% models(m).prettylabel  = 'simple shadow-rate VAR';
% models(m).shortlabel   = 'ShadowRateVAR-p12';
% models(m).fcstType     = 'fcstY';

% BLOCK-HYBRID RATE
m = m + 1;
models(m).resultlabel  = 'ELBblockhybrid-RATSbvarshrinkage-p12';
models(m).prettylabel  = 'hybrid shadow-rate VAR';
models(m).shortlabel   = 'BlockHybridVAR-p12';
models(m).fcstType     = 'fcstY';


%% loop over altlabels
datalabel0   = strcat(datalabel0, '-2022-09');

for dd = 1 : length(ALTDATA)

    datalabel1   = strcat(ALTDATA(dd).dataset, '-2022-09');
    prettylabel1 = ALTDATA(dd).prettylabel;

    %% one big wrapper

    titlename = sprintf('oosCompareDatasetsEvaluationTables-%s-vs-%s', datalabel0, datalabel1);
    initwrap

    %% loop over model sets

    for m = 1 : length(models)


        %% eval window

        for sam = 1 : length(samples)

            evalStart  = samples(sam).evalStart;
            evalStop   = samples(sam).evalStop;

            evaltxt = sprintf('evalStart%sevalEnd%s', datestr(evalStart, 'yyyymm'), datestr(evalStop, 'yyyymm'));

            %% load data
            clear oos0 oos1 ydates Tjumpoffs

            oos0 = matfile(fullfile(resultsdir, sprintf('%s-%s.mat', datalabel0, 'standardVAR-RATSbvarshrinkage-p12')));
            oos1 = matfile(fullfile(resultsdir, sprintf('%s-%s.mat', datalabel1, models(m).resultlabel)));

            if oos0.MCMCdraws ~= oos1.MCMCdraws
                warning('unequal numbers of MCMCdraws, 0 has %d, 1 has %d', oos0.MCMCdraws, oos1.MCMCdraws)
            end


            %% check for identical samples
            if oos0.ydates ~= oos1.ydates
                error('oos estimates based on different samples')
            end
            ydates    = oos0.ydates;

            if ~isequal(oos0.Tjumpoffs, oos1.Tjumpoffs)
                error('oos jumpoffs differ')
            end

            Tjumpoffs = oos0.Tjumpoffs;

            %% cut eval sample if desired
            ndxJumpoff = ismember(Tjumpoffs, find((ydates >= evalStart) & (ydates <= evalStop)));
            Tjumpoffs  = Tjumpoffs(ndxJumpoff);


            dates    = ydates(Tjumpoffs);

            comparisonNote = sprintf('Baseline estimated from linear VAR, alternative forecasts generated by %s. Evaluation window with forecast origins from %s through %s (and outcome data as far as available).', ...
                models(m).prettylabel, datestr(dates(1), 'yyyy:mm'), datestr(dates(end), 'yyyy:mm'));
            shortcomparisonNote = sprintf('from %s through %s.', ...
                datestr(dates(1), 'yyyy:mm'), datestr(dates(end), 'yyyy:mm'));


            %% some parameters
            Nhorizons  = min(oos0.fcstNhorizons,oos1.fcstNhorizons);

            %% find common set of variables
            ncode0 = oos0.ncode;
            ncode1 = oos1.ncode;

            [ncode, ndxY0, ndxY1] = intersect(ncode0, ncode1, 'stable');
            if doCharts
                Ylabels = fredMDshortlabel(ncode);
            else
                Ylabels = fredMDprettylabel(ncode);
            end
            N       = length(Ylabels);

            Ylabels = strrep(Ylabels, '_', '');

            %% setup monthly tables
            if doCharts
                theseHorizons = [3 12 24];
            else
                theseHorizons = [3 6 12 24];
            end

            %% RMSE
            losstype0 = sprintf('%shaterror', models(m).fcstType);
            losstype1 = sprintf('%shaterror', models(m).fcstType);
            mseloss0 = oos0.(losstype0).^2;
            mseloss1 = oos1.(losstype1).^2;

            % match samples
            mseloss0 = mseloss0(ndxY0,:,ndxJumpoff);
            mseloss1 = mseloss1(ndxY1,:,ndxJumpoff);

            RMSE0 = sqrt(nanmean(mseloss0,3));
            RMSE1 = sqrt(nanmean(mseloss1,3));
            relativeRMSE01 =  RMSE1(:,1:Nhorizons) ./ RMSE0(:,1:Nhorizons); % here: RMSE


            %% MAE
            losstype0 = sprintf('%smederror', models(m).fcstType);
            losstype1 = sprintf('%smederror', models(m).fcstType);

            maeloss0     = abs(oos0.(losstype0));
            maeloss1     = abs(oos1.(losstype1));

            % match samples
            maeloss0 = maeloss0(ndxY0,:,ndxJumpoff);
            maeloss1 = maeloss1(ndxY1,:,ndxJumpoff);

            mae0      = nanmean(maeloss0,3);
            mae1      = nanmean(maeloss1,3);
            relativeMAD01 = mae1(:,1:Nhorizons) ./ mae0(:,1:Nhorizons); % here: RMAE


            %% CRPS
            losstype0     = sprintf('%scrps', models(m).fcstType);
            losstype1     = sprintf('%scrps', models(m).fcstType);

            crpsloss0     = oos0.(losstype0);
            crpsloss1     = oos1.(losstype1);
            % match samples
            crpsloss0 = crpsloss0(ndxY0,:,ndxJumpoff);
            crpsloss1 = crpsloss1(ndxY1,:,ndxJumpoff);


            crps1          = nanmean(crpsloss1,3);
            crps0          = nanmean(crpsloss0,3);
            relativeCRPS01 = crps1(:,1:Nhorizons) ./ crps0(:,1:Nhorizons); % here: Relative CRPS

            %% compare all
            statlabels = {'RMSE', 'MAE', 'CRPS'};
            if doCharts
                tabname = sprintf('comparedatasetForecastsChart-%s-%s-vs-%s-%s.tex', models(m).shortlabel, ...
                    datalabel0, datalabel1, evaltxt);
            else
                tabname = sprintf('comparedatasetForecasts-%s-%s-vs-%s-%s.tex', models(m).shortlabel, ...
                    datalabel0, datalabel1, evaltxt);
            end
            tabcaption =  sprintf('%s %s', models(m).shortlabel, shortcomparisonNote);

            compareAllinone(tabname, wrap, doCharts, doBold, ...
                mseloss0, mseloss1, relativeRMSE01, ...
                maeloss0, maeloss1, relativeMAD01, ...
                crpsloss0, crpsloss1, relativeCRPS01, ...
                prettylabel0, prettylabel1, ...
                Ylabels, theseHorizons, tabcaption, statlabels, comparisonNote)


        end
    end

    %% finish script
    finishwrap
end % altlabel loop

finishscript

%% helper function to create tex table
function compareAllinone(tabname, wrap, doCharts, doBold, ...
    mseloss0, mseloss1, relativeRMSE01, ...
    maeloss0, maeloss1, relativeMAE01, ...
    crpsloss0, crpsloss1, relativeCRPS01, ...
    prettylabel0, prettylabel1, ...
    Ylabels, theseHorizons, tabcaption, statlabels, comparisonNote)


%% parse inputs
N = length(Ylabels);
Nhorizons = length(theseHorizons);

%% DM tests

[relativeRMSE01, dmMSEtstat] = dodm(mseloss0, mseloss1, relativeRMSE01, theseHorizons);

[relativeMAE01, dmMADtstat] = dodm(maeloss0, maeloss1, relativeMAE01, theseHorizons);

[relativeCRPS01, dmCRPStstat] = dodm(crpsloss0, crpsloss1, relativeCRPS01, theseHorizons);


%% set up tab
if ~isempty(wrap)
    tabdir = wrap.dir;
    latexwrapper(wrap, 'add', 'sidetab', tabname, tabcaption)
end

%% tabulate
fid = fopen(fullfile(tabdir, tabname), 'wt');
if doCharts
    fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('.3', 1, 3 * Nhorizons));
    fprintf(fid, ' & \\multicolumn{%d}{c}{\\bf %s}  & \\multicolumn{%d}{c}{\\bf %s}  & \\multicolumn{%d}{c}{\\bf %s} \\\\ \\cmidrule(lr){%d-%d}\\cmidrule(lr){%d-%d}\\cmidrule(lr){%d-%d} \n', ...
        Nhorizons, statlabels{1}, Nhorizons, statlabels{2}, Nhorizons, statlabels{3}, ...
        1+1,1+Nhorizons,1+Nhorizons+1,1+2*Nhorizons,1+2*Nhorizons+1, 1+3*Nhorizons);
    fprintf(fid, 'Var. / Hor. ');

else
    fprintf(fid, '\\begin{center}\n');
    fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('.4', 1, 3 * Nhorizons));
    fprintf(fid, ' & \\multicolumn{%d}{c}{%s}  & \\multicolumn{%d}{c}{%s}  & \\multicolumn{%d}{c}{%s} \\\\ \\cmidrule(lr){%d-%d}\\cmidrule(lr){%d-%d}\\cmidrule(lr){%d-%d} \n', ...
        Nhorizons, statlabels{1}, Nhorizons, statlabels{2}, Nhorizons, statlabels{3}, ...
        1+1,1+Nhorizons,1+Nhorizons+1,1+2*Nhorizons,1+2*Nhorizons+1, 1+3*Nhorizons);
    fprintf(fid, 'Variable / Horizon ');
end
for h = 1 : Nhorizons
    fprintf(fid, '& \\multicolumn{1}{c}{$%d$} ', theseHorizons(h));
end
for h = 1 : Nhorizons
    fprintf(fid, '& \\multicolumn{1}{c}{$%d$} ', theseHorizons(h));
end
for h = 1 : Nhorizons
    fprintf(fid, '& \\multicolumn{1}{c}{$%d$} ', theseHorizons(h));
end
fprintf(fid, '\\\\\n');
fprintf(fid, '\\midrule\n');
for n = 1 : N
    fprintf(fid, '%s ', Ylabels{n});
    for h = 1 : Nhorizons
        if isfinite(relativeRMSE01(n,h))
            if doCharts
                switch doColorCode(relativeRMSE01(n,h))
                    case 1
                        fprintf(fid, '& %s ', dcolred(sprintf('%6.2f%s ', relativeRMSE01(n,h), Zstar(dmMSEtstat(n,h)))));
                    case -1
                        fprintf(fid, '& %s ', dcolgreen(sprintf('%6.2f%s ', relativeRMSE01(n,h), Zstar(dmMSEtstat(n,h)))));
                    otherwise
                        fprintf(fid, '& %6.2f%s ', relativeRMSE01(n,h), Zstar(dmMSEtstat(n,h)));
                end
            else
                if doBold && doColorCode(relativeRMSE01(n,h))
                    fprintf(fid, '& %s%s ', dcolbf(relativeRMSE01(n,h), '%6.2f'), Zstar(dmMSEtstat(n,h)));
                else
                    fprintf(fid, '& %6.2f%s ', relativeRMSE01(n,h), Zstar(dmMSEtstat(n,h)));
                end
            end
        else
            fprintf(fid, '& -.- ');
        end
    end
    for h = 1 : Nhorizons
        if isfinite(relativeMAE01(n,h))
            if doCharts
                switch doColorCode(relativeMAE01(n,h))
                    case 1
                        fprintf(fid, '& %s ', dcolred(sprintf('%6.2f%s ', relativeMAE01(n,h), Zstar(dmMADtstat(n,h)))));
                    case -1
                        fprintf(fid, '& %s ', dcolgreen(sprintf('%6.2f%s ', relativeMAE01(n,h), Zstar(dmMADtstat(n,h)))));
                    otherwise
                        fprintf(fid, '& %6.2f%s ', relativeMAE01(n,h), Zstar(dmMADtstat(n,h)));
                end
            else
                if doBold && doColorCode(relativeMAE01(n,h))
                    fprintf(fid, '& %s%s ', dcolbf(relativeMAE01(n,h), '%6.2f'), Zstar(dmMADtstat(n,h)));
                else
                    fprintf(fid, '& %6.2f%s ', relativeMAE01(n,h), Zstar(dmMADtstat(n,h)));
                end
            end
        else
            fprintf(fid, '& -.- ');
        end
    end
    for h = 1 : Nhorizons
        if isfinite(relativeCRPS01(n,h))
            if doCharts
                switch doColorCode(relativeCRPS01(n,h))
                    case 1
                        fprintf(fid, '& %s ', dcolred(sprintf('%6.2f%s ', relativeCRPS01(n,h), Zstar(dmCRPStstat(n,h)))));
                    case -1
                        fprintf(fid, '& %s ', dcolgreen(sprintf('%6.2f%s ', relativeCRPS01(n,h), Zstar(dmCRPStstat(n,h)))));
                    otherwise
                        fprintf(fid, '& %6.2f%s ', relativeCRPS01(n,h), Zstar(dmCRPStstat(n,h)));
                end
            else
                if doBold && doColorCode(relativeCRPS01(n,h))
                    fprintf(fid, '& %s%s ', dcolbf(relativeCRPS01(n,h), '%6.2f'), Zstar(dmCRPStstat(n,h)));
                else
                    fprintf(fid, '& %6.2f%s ', relativeCRPS01(n,h), Zstar(dmCRPStstat(n,h)));
                end
            end
        else
            fprintf(fid, '& -.- ');
        end
    end
    fprintf(fid, '\\\\\n');
end
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
if ~doCharts
    fprintf(fid, '\\end{center}\n');
end
fprintf(fid, '\n');

if ~doCharts
    sigone = cat(1, abs(dmMSEtstat) > norminv(0.95, 0, 1) & (round(relativeRMSE01,2) == 1), ...
        abs(dmMADtstat) > norminv(0.95, 0, 1) & (round(relativeMAE01,2) == 1), ...
        abs(dmCRPStstat) > norminv(0.95, 0, 1) & (round(relativeCRPS01,2) == 1));

    fprintf(fid, 'Note: Comparison of ``%s'''' (baseline, in denominator) against ``%s'''' for horizons', ...
        prettylabel0, prettylabel1);
    fprintf(fid, ' %d, ', theseHorizons(1:end-1));
    fprintf(fid, 'and %d.\n', theseHorizons(end));
    fprintf(fid, 'Values below 1 indicate improvement over baseline. \n');
    fprintf(fid, '%s \n', comparisonNote);

    fprintf(fid, 'Significance assessed by Diebold-Mariano-West test using Newey-West standard errors with $h + 1$ lags.\n');

    if any(sigone, 'all')
        if sum(sigone(:)) > 1
            fprintf(fid, 'Due to the close behavior of some of the models compared, and rounding of the reported values, a few comparisons show significant ratios  of 1.00.\n');
            fprintf(fid, 'These cases arise from persistent differences in performance that are, however, too small to be relevant after rounding.\n');
        else
            fprintf(fid, 'Due to the close behavior of some of the models compared, and rounding of the reported values, one of the comparisons shows a significant ratio of 1.00.\n');
            fprintf(fid, 'This case arises from persistent differences in performance that are, however, too small to be relevant after rounding.\n');
        end
    end

    if doBold
        fprintf(fid, 'Performance differences of 5 percent and more (relative to baseline) are indicated by bold face numbers.\n');
    end

    if ~all(isfinite(relativeMAE01(:)))
        fprintf(fid, 'In some cases, due to strong performance of the baseline model, relative MAD may involve divisions by zero. These cases are reported as blank entries.');
    end
end
fclose(fid);

type(fullfile(tabdir, tabname))

end

function flag = doColorCode(x)

if round(x,2) >= 1.05
    flag = 1;
elseif round(x,2) <= .95
    flag = -1;
else
    flag = 0;
end

end % function

function [deltaLoss, tstat] = dodm(loss0, loss1, deltaLoss, theseHorizons)

loss0 = loss0(:,theseHorizons,:);
loss1 = loss1(:,theseHorizons,:);
deltaLoss = deltaLoss(:,theseHorizons);

[N, Nhorizons,~] = size(loss0);

tstat = NaN(N,Nhorizons);

for h = 1 : Nhorizons
    nwLag = theseHorizons(h) + 1;
    for n = 1 : N
        thisloss0 = squeeze(loss0(n,h,:));
        thisloss1 = squeeze(loss1(n,h,:));

        if isequaln(thisloss0, thisloss1) || any(isinf(thisloss0)) || any(isinf(thisloss1))
            % do noting
        else
            [~,tstat(n,h)] = dmtest(thisloss0,thisloss1, nwLag);
        end
    end
end

end