%% compare OOS results from two estimates
% load quantico*.mat files and assess OOS performance

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

resultsdir = '~/jam/lager/quantico2023logscoresXL/';

doCharts = false;
doBold   = true; % matter only if doCharts = false;

datalabel = 'fredMD20-2022-09';

for p = 12 % [12 6 3]

    %% one wrapper per lag choice

    titlename = sprintf('oosPairwiseEvaluationTables-%s-p%d', datalabel, p);
    initwrap



    %% BASELINE p=12

    m0 = 0;

    % STANDARD
    m = m0 +  1;
    models(m).datalabel    = datalabel; %#ok<*SAGROW>
    models(m).resultlabel  = sprintf('standardVAR-RATSbvarshrinkage-p%d', p);
    if p == 12
        models(m).prettylabel  = sprintf('standard linear VAR');
    else
        models(m).prettylabel  = sprintf('standard linear VAR (p=%d)', p);
    end
    models(m).shortlabel   = sprintf('Standard-p%d', p);
    models(m).fcstType     = 'fcstY';

    % CENSORED
    m = m0 +  2;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = sprintf('standardVAR-RATSbvarshrinkage-p%d', p);
    if p == 12
        models(m).prettylabel  = sprintf('Censored');
    else
        models(m).prettylabel  = sprintf('Censored (p=%d)', p);
    end
    models(m).shortlabel   = sprintf('CensoredVAR-p%d', p);
    models(m).fcstType     = 'fcstYcensor';


    % QUASI-SHADOW
    m = m0 +  3;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = sprintf('standardVAR-RATSbvarshrinkage-p%d', p);
    if p == 12
        models(m).prettylabel  = sprintf('Quasi-shadow-rate');
    else
        models(m).prettylabel  = sprintf('Quasi-shadow-rate (p=%d)', p);
    end
    models(m).shortlabel   = sprintf('QuasiShadowVAR-p%d', p);
    models(m).fcstType     = 'fcstYshadow';

    % SHADOW-RATE
    m = m0 +  4;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = sprintf('ELBsampling-RATSbvarshrinkage-p%d', p);
    if p == 12
        models(m).prettylabel  = sprintf('simple shadow-rate VAR');
    else
        models(m).prettylabel  = sprintf('simple shadow-rate VAR (p=%d)', p);
    end
    models(m).shortlabel   = sprintf('ShadowRateVAR-p%d', p);
    models(m).fcstType     = 'fcstY';

    % SHADOW-RATE (Censored Yields)
    m = m0 +  5;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = sprintf('ELBsampling-RATSbvarshrinkage-p%d', p);
    if p == 12
        models(m).prettylabel  = sprintf('Shadow-rate (w/censored yields');
    else
        models(m).prettylabel  = sprintf('Shadow-rate (w/censored yields, p=%d)', p);
    end
    models(m).shortlabel   = sprintf('ShadowRateCensoredVAR-p%d', p);
    models(m).fcstType     = 'fcstYcensor';

    % BLOCK-HYBRID RATE
    m = m0 +  6;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = sprintf('ELBblockhybrid-RATSbvarshrinkage-p%d', p);
    if p == 12
        models(m).prettylabel  = sprintf('hybrid shadow-rate VAR');
    else
        models(m).prettylabel  = sprintf('hybrid shadow-rate VAR (p=%d)', p);
    end
    models(m).shortlabel   = sprintf('BlockHybridVAR-p%d', p);
    models(m).fcstType     = 'fcstY';


    if strcmpi(datalabel, 'fredMD20-2022-09')

        % BLOCK HYBRID w/Aelb
        m = m0 +  7;
        models(m).datalabel    = datalabel;
        models(m).resultlabel  = sprintf('ELBblockhybridAelbPrec10000-RATSbvarshrinkage-p%d', p);
        if p == 12
            models(m).prettylabel  = sprintf('hybrid shadow-rate VAR (A switching)');
        else
            models(m).prettylabel  = sprintf('hybrid shadow-rate VAR (A switching, p=%d)', p);
        end
        models(m).shortlabel   = sprintf('BlockHybridVARswitchingA10000-p%d', p);
        models(m).fcstType     = 'fcstY';

        % KRIPPNER
        m = m0 +  8;
        models(m).datalabel    = 'fredMD20krippner-2022-09';
        models(m).resultlabel  = sprintf('standardVAR-RATSbvarshrinkage-p%d', p);
        if p == 12
            models(m).prettylabel  = sprintf('Plug-in VAR (Krippner)');
        else
            models(m).prettylabel  = sprintf('Plug-in VAR (Krippner, p=%d)', p);
        end
        models(m).shortlabel   = sprintf('krippnerVAR-p%d', p);
        models(m).fcstType     = 'fcstYshadow';

        % WUXIA
        m = m0 +  9;
        models(m).datalabel    = 'fredMD20wuxia-2022-09';
        models(m).resultlabel  = sprintf('standardVAR-RATSbvarshrinkage-p%d', p);
        if p == 12
            models(m).prettylabel  = sprintf('Plug-in VAR (WuXia)');
        else
            models(m).prettylabel  = sprintf('Plug-in VAR (WuXia, p=%d)', p);
        end
        models(m).shortlabel   = sprintf('wuxiaVAR-p%d', p);
        models(m).fcstType     = 'fcstYshadow';
    end


    %% select models
    if strcmpi(datalabel, 'fredMD20-2022-09')
        MODELS = {[1 4] [1 6] [4 6] [1 7] [6 7]...
            [8 1] [9 1] [8 4] [9 4] [8 6] [9 6]
            };

        MLABELS = {'Linear vs Shadow-rate VAR', 'Linear vs Block-Hybrid VAR', 'Shadow-rate vs Block-Hybrid VAR', 'Linear vs Block-Hybrid A-Switching VAR', 'Block-Hybrid  vs Block-Hybrid A-Switching VAR', ...
            'Kripper-plug-in vs Linear VAR', 'Wu-Xia-plug-in vs Linear VAR', ...
            'Kripper-plug-in vs Shadow-rate VAR', 'Wu-Xia-plug-in vs Shadow-rate VAR', ...
            'Kripper-plug-in vs Block-Hybrid VAR', 'Wu-Xia-plug-in vs Block-Hybrid VAR', ...
            };
    else
        MODELS = {[1 4] [1 6] [4 6] ...
            };

        MLABELS = {'Linear vs Shadow-rate VAR', 'Linear vs Block-Hybrid VAR', 'Shadow-rate vs Block-Hybrid VAR', ...
            };
    end

    %% loop over model sets

    for m = 1 : length(MODELS)



        m0 = MODELS{m}(1);
        m1 = MODELS{m}(2);


        %% eval window

        for sam = 1 : 4

            switch sam
                case 1
                    % baseline
                    evalStart  = datenum(2009,1,1);
                    evalStop   = datenum(2022,8,1);
                case 2
                    % ex COVID
                    evalStart = datenum(2009,1,1);
                    evalStop  = datenum(2017,12,1);
                case 3
                    % COVID
                    evalStart = datenum(2018,1,1);
                    evalStop  = datenum(2022,8,1);
                case 4
                    % first ELB
                    evalStart  = datenum(2009,1,1);
                    evalStop   = datenum(2015,12,1);
                otherwise
                    error('sam %d not defined', sam)
            end


            evaltxt = sprintf('evalStart%sevalEnd%s', datestr(evalStart, 'yyyymm'), datestr(evalStop, 'yyyymm'));




            %#ok<*UNRCH>
            %% load data
            clear oos0 oos1 oos2 ydates Tjumpoffs

            oos0 = matfile(fullfile(resultsdir, sprintf('%s-%s.mat', models(m0).datalabel, models(m0).resultlabel)));
            oos1 = matfile(fullfile(resultsdir, sprintf('%s-%s.mat', models(m1).datalabel, models(m1).resultlabel)));

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

            comparisonNote = sprintf('Evaluation window with forecast origins from %s through %s (and outcome data as far as available).', ...
                datestr(dates(1), 'yyyy:mm'), datestr(dates(end), 'yyyy:mm'));
            shortcomparisonNote = sprintf('from %s through %s.', ...
                datestr(dates(1), 'yyyy:mm'), datestr(dates(end), 'yyyy:mm'));


            %% some parameters
            Nhorizons  = oos0.fcstNhorizons;

            if doCharts
                Ylabels = fredMDshortlabel(oos0.ncode);
            else
                Ylabels = fredMDprettylabel(oos0.ncode);
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
            losstype0 = sprintf('%shaterror', models(m0).fcstType);
            losstype1 = sprintf('%shaterror', models(m1).fcstType);
            mseloss0 = oos0.(losstype0).^2;
            mseloss1 = oos1.(losstype1).^2;

            % match samples
            mseloss0 = mseloss0(:,:,ndxJumpoff);
            mseloss1 = mseloss1(:,:,ndxJumpoff);

            RMSE0 = sqrt(nanmean(mseloss0,3));
            RMSE1 = sqrt(nanmean(mseloss1,3));
            relativeRMSE01 =  RMSE1 ./ RMSE0; % here: RMSE


            %% MAE
            losstype0 = sprintf('%smederror', models(m0).fcstType);
            losstype1 = sprintf('%smederror', models(m1).fcstType);

            maeloss0     = abs(oos0.(losstype0));
            maeloss1     = abs(oos1.(losstype1));

            % match samples
            maeloss0 = maeloss0(:,:,ndxJumpoff);
            maeloss1 = maeloss1(:,:,ndxJumpoff);

            mae0      = nanmean(maeloss0,3);
            mae1      = nanmean(maeloss1,3);
            relativeMAD01 = mae1 ./ mae0; % here: RMAE


            %% CRPS
            losstype0     = sprintf('%scrps', models(m0).fcstType);
            losstype1     = sprintf('%scrps', models(m1).fcstType);

            crpsloss0     = oos0.(losstype0);
            crpsloss1     = oos1.(losstype1);
            % match samples
            crpsloss0 = crpsloss0(:,:,ndxJumpoff);
            crpsloss1 = crpsloss1(:,:,ndxJumpoff);


            relativeCRPS01 = nanmean(crpsloss1,3) ./ nanmean(crpsloss0,3); % here: Relative CRPS

            %% compare all
            statlabels = {'RMSE', 'MAE', 'CRPS'};
            if doCharts
                tabname = sprintf('allinoneChart-%s-%s-vs-%s-%s.tex', models(m0).datalabel, ...
                    models(m0).shortlabel, models(m1).shortlabel, evaltxt);
            else
                tabname = sprintf('allinone-%s-%s-vs-%s-%s.tex', models(m0).datalabel, ...
                    models(m0).shortlabel, models(m1).shortlabel, evaltxt);
            end
            tabcaption =  sprintf('%s %s', MLABELS{m}, shortcomparisonNote);

            compareAllinone(tabname, wrap, doCharts, doBold, ...
                mseloss0, mseloss1, relativeRMSE01, ...
                maeloss0, maeloss1, relativeMAD01, ...
                crpsloss0, crpsloss1, relativeCRPS01, ...
                models(m0).prettylabel, models(m1).prettylabel, ...
                Ylabels, theseHorizons, tabcaption, statlabels, comparisonNote)


        end
    end

    %% finish wrap
    finishwrap

end % p

%% finish script
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

    fprintf(fid, 'Note: Comparison of ``%s'''' (baseline, in denominator) against ``%s.'''' Values below 1 indicate improvement over baseline. \n', ...
        prettylabel0, prettylabel1);
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