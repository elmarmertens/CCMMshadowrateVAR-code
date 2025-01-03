%% compare OOS results from two pairs of estimates
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


resultsdir = pwd;


for DATALABEL = {'fredsxMD20exYield-2022-09', 'fredsxMD20-2022-09'}

    datalabel = DATALABEL{:};

    %% one wrapper per data choice

    titlename = sprintf('oosTripleEvaluationTablesQE-%s', datalabel);
    initwrap



    %% BASELINE p=12

    doBold      = true;
    doFlipOrder = true;
    p           = 12;

    m0 = 0;

    % STANDARD
    m = m0 +  1;
    models(m).datalabel    = datalabel; %#ok<*SAGROW>
    models(m).resultlabel  = sprintf('standardVARAR1SV-RATSbvarshrinkage-p%d', p);
    if p == 12
        models(m).prettylabel  = sprintf('standard linear VAR');
        models(m).prettyshortlabel = 'Standard';
    else
        models(m).prettylabel  = sprintf('standard linear VAR (p=%d)', p);
        models(m).prettyshortlabel = sprintf('Standard (p=%d)', p);
    end
    models(m).shortlabel   = sprintf('Standard-p%d', p);
    models(m).fcstType     = 'fcstY';

    % CENSORED
    m = m0 +  2;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = sprintf('standardVARAR1SV-RATSbvarshrinkage-p%d', p);
    if p == 12
        models(m).prettylabel  = sprintf('Censored');
        models(m).prettyshortlabel = 'Censored';
    else
        models(m).prettylabel  = sprintf('Censored (p=%d)', p);
        models(m).prettyshortlabel = sprintf('Censored (p=%d)', p);
    end
    models(m).shortlabel   = sprintf('CensoredVAR-p%d', p);
    models(m).fcstType     = 'fcstYcensor';


    % QUASI-SHADOW
    m = m0 +  3;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = sprintf('standardVARAR1SV-RATSbvarshrinkage-p%d', p);
    if p == 12
        models(m).prettylabel  = sprintf('Quasi-shadow-rate');
        models(m).prettyshortlabel = 'Quasi-shadow-rate';
    else
        models(m).prettylabel  = sprintf('Quasi-shadow-rate (p=%d)', p);
        models(m).prettyshortlabel = sprintf('Quasi-shadow-rate (p=%d)', p);
    end
    models(m).shortlabel   = sprintf('QuasiShadowVAR-p%d', p);
    models(m).fcstType     = 'fcstYshadow';

    % SHADOW-RATE
    m = m0 +  4;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = sprintf('ELBsamplingAR1SV-RATSbvarshrinkage-p%d', p);
    if p == 12
        models(m).prettylabel  = sprintf('simple shadow-rate VAR');
        models(m).prettyshortlabel = 'simple';
    else
        models(m).prettylabel  = sprintf('simple shadow-rate VAR (p=%d)', p);
        models(m).prettyshortlabel = sprintf('simple (p=%d)', p);
    end
    models(m).shortlabel   = sprintf('ShadowRateVAR-p%d', p);
    models(m).fcstType     = 'fcstY';

    % SHADOW-RATE (Censored Yields)
    m = m0 +  5;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = sprintf('ELBsamplingAR1SV-RATSbvarshrinkage-p%d', p);
    if p == 12
        models(m).prettylabel  = sprintf('Shadow-rate (w/censored yields');
        models(m).prettyshortlabel = 'Censored yields';
    else
        models(m).prettylabel  = sprintf('Shadow-rate (w/censored yields, p=%d)', p);
        models(m).prettyshortlabel = sprintf('Censored yields (p=%d)', p);
    end
    models(m).shortlabel   = sprintf('ShadowRateCensoredVAR-p%d', p);
    models(m).fcstType     = 'fcstYcensor';

    % BLOCK-HYBRID RATE
    m = m0 +  6;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = sprintf('ELBblocknonstructuralAR1SV-RATSbvarshrinkage-p%d', p);
    if p == 12
        models(m).prettylabel       = sprintf('Restricted non-structural shadow-rate VAR');
        models(m).prettyshortlabel  = sprintf('Restricted');
    else
        models(m).prettylabel       = sprintf('Restricted non-structural shadow-rate VAR (p=%d)', p);
        models(m).prettyshortlabel  = sprintf('Restricted (p=%d)', p);
    end
    models(m).shortlabel   = sprintf('blocknonstructuralVAR-p%d', p);
    models(m).fcstType     = 'fcstY';


    % HYBRID
    m = m0 +  7;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = sprintf('ELBnonstructuralAR1SV-RATSbvarshrinkage-p%d', p);
    if p == 12
        models(m).prettylabel  = sprintf('Non-structural shadow-rate VAR');
        models(m).prettyshortlabel  = sprintf('Non-structural');
    else
        models(m).prettylabel  = sprintf('Non-structural shadow-rate VAR (p=%d)', p);
        models(m).prettyshortlabel  = sprintf('Non-structural (p=%d)', p);
    end
    models(m).shortlabel   = sprintf('nonstructuralVAR-p%d', p);
    models(m).fcstType     = 'fcstY';

    % General HYBRID
    m = m0 +  8;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = sprintf('ELBshadowrateGeneralAR1SV-RATSbvarshrinkage-p%d', p);
    if p == 12
        models(m).prettylabel = sprintf('General shadow-rate VAR');
        models(m).prettyshortlabel  = sprintf('General');
    else
        models(m).prettylabel = sprintf('General shadow-rate VAR (p=%d)', p);
        models(m).prettyshortlabel  = sprintf('General (p=%d)', p);
    end
    models(m).shortlabel   = sprintf('shadowrateGeneralVAR-p%d', p);
    models(m).fcstType     = 'fcstY';

    %% compare models

    if strcmpi(datalabel, 'fredsxMD20exYield-2022-09')
        TRIPLETS = {[1 8 6] [1 8 7] [1 7 6]};
    else
        TRIPLETS = {[1 7 6]};
    end

    for tt = 1 : length(TRIPLETS)

        m0 = TRIPLETS{tt}(1);
        m1 = TRIPLETS{tt}(2);
        m2 = TRIPLETS{tt}(3);

        thisTRIPLETlabel = sprintf('%s and %s', models(m1).prettylabel, models(m2).prettylabel);

        %% eval window

        for sam = [3 4]

            switch sam
                case 1
                    % baseline
                    evalStart  = datenum(2009,1,1);
                    evalStop   = datenum(2022,8,1);
                case 2
                    % ex COVID
                    evalStart = datenum(2009,1,1);
                    evalStop  = datenum(2017,12,1);
                case 3 % baseline since 2010
                    evalStart  = datenum(2010,1,1);
                    evalStop   = datenum(2022,8,1);
                case 4
                    % ex COVID since 2010
                    evalStart = datenum(2010,1,1);
                    evalStop  = datenum(2017,12,1);
                otherwise
                    error('sam %d not defined', sam)
            end


            evaltxt = sprintf('evalStart%sevalEnd%s', datestr(evalStart, 'yyyymm'), datestr(evalStop, 'yyyymm'));




            %#ok<*UNRCH>
            %% load data
            clear oos0 oos1 oos2 ydates Tjumpoffs

            oos0 = matfile(fullfile(resultsdir, sprintf('%s-%s.mat', models(m0).datalabel, models(m0).resultlabel)));
            oos1 = matfile(fullfile(resultsdir, sprintf('%s-%s.mat', models(m1).datalabel, models(m1).resultlabel)));
            oos2 = matfile(fullfile(resultsdir, sprintf('%s-%s.mat', models(m2).datalabel, models(m2).resultlabel)));

            if oos0.MCMCdraws ~= oos1.MCMCdraws
                warning('unequal numbers of MCMCdraws, 0 has %d, 1 has %d', oos0.MCMCdraws, oos1.MCMCdraws)
            end
            if oos0.MCMCdraws ~= oos2.MCMCdraws
                warning('unequal numbers of MCMCdraws, 0 has %d, 2 has %d', oos0.MCMCdraws, oos2.MCMCdraws)
            end

            %% check for identical samples
            if oos0.ydates ~= oos1.ydates
                error('oos estimates based on different samples')
            end
            if oos0.ydates ~= oos2.ydates
                error('oos estimates based on different samples')
            end
            ydates    = oos0.ydates;

            if ~isequal(oos0.Tjumpoffs, oos1.Tjumpoffs)
                error('oos jumpoffs differ')
            end
            if ~isequal(oos0.Tjumpoffs, oos2.Tjumpoffs)
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
            Nhorizons  = min(oos0.fcstNhorizons,oos1.fcstNhorizons);
            Nhorizons  = min(Nhorizons, oos2.fcstNhorizons);

            ncode = oos0.ncode;

            Ylabels = fredMDveryshortlabel(oos0.ncode);
            N       = length(Ylabels);

            Ylabels = strrep(Ylabels, '_', '');

            %% setup monthly tables
            theseHorizons = [6 12 24]; % QE

            %% MAE
            losstype0 = sprintf('%smederror', models(m0).fcstType);
            losstype1 = sprintf('%smederror', models(m1).fcstType);
            losstype2 = sprintf('%smederror', models(m2).fcstType);

            maeloss0     = abs(oos0.(losstype0));
            maeloss1     = abs(oos1.(losstype1));
            maeloss2     = abs(oos2.(losstype2));

            % match samples
            maeloss0 = maeloss0(:,:,ndxJumpoff);
            maeloss1 = maeloss1(:,:,ndxJumpoff);
            maeloss2 = maeloss2(:,:,ndxJumpoff);

            mae0      = nanmean(maeloss0,3);
            mae1      = nanmean(maeloss1,3);
            relativeMAD01 = mae1(:,1:Nhorizons) ./ mae0(:,1:Nhorizons); % here: RMAE
            mae2      = nanmean(maeloss2,3);
            relativeMAD02 = mae2(:,1:Nhorizons) ./ mae0(:,1:Nhorizons); % here: RMAE

            %% CRPS
            losstype0     = sprintf('%scrps', models(m0).fcstType);
            losstype1     = sprintf('%scrps', models(m1).fcstType);
            losstype2     = sprintf('%scrps', models(m2).fcstType);

            crpsloss0     = oos0.(losstype0);
            crpsloss1     = oos1.(losstype1);
            crpsloss2     = oos2.(losstype2);

            % match samples
            crpsloss0 = crpsloss0(:,:,ndxJumpoff);
            crpsloss1 = crpsloss1(:,:,ndxJumpoff);
            crpsloss2 = crpsloss2(:,:,ndxJumpoff);

            crps1          = nanmean(crpsloss1,3);
            crps0          = nanmean(crpsloss0,3);
            relativeCRPS01 = crps1(:,1:Nhorizons) ./ crps0(:,1:Nhorizons); % here: Relative CRPS
            crps2          = nanmean(crpsloss2,3);
            relativeCRPS02 = crps2(:,1:Nhorizons) ./ crps0(:,1:Nhorizons); % here: Relative CRPS

            %% flip order of variables
            if doFlipOrder
                relativeMAD01  = flip(relativeMAD01,1);
                relativeMAD02  = flip(relativeMAD02,1);
                maeloss0       = flip(maeloss0,1);
                maeloss1       = flip(maeloss1,1);
                maeloss2       = flip(maeloss2,1);
                relativeCRPS01 = flip(relativeCRPS01,1);
                relativeCRPS02 = flip(relativeCRPS02,1);
                crpsloss0      = flip(crpsloss0,1);
                crpsloss1      = flip(crpsloss1,1);
                crpsloss2      = flip(crpsloss2,1);
                Ylabels        = flip(Ylabels);
            end
            %% compare all
            statlabels = {'MAE', 'CRPS'};
            tabname = sprintf('tripleComparisonQE-%s-%s-vs-%s-vs-%s-%s.tex', models(m0).datalabel, ...
                models(m0).shortlabel, models(m1).shortlabel, models(m2).shortlabel, evaltxt);
            tabcaption =  sprintf('%s %s', thisTRIPLETlabel, shortcomparisonNote);

            compareTriple(tabname, wrap, doBold, ...
                maeloss0, maeloss1, relativeMAD01, maeloss2, relativeMAD02, ...
                crpsloss0, crpsloss1, relativeCRPS01, crpsloss2, relativeCRPS02, ...
                models(m0).prettylabel, models(m1).prettylabel, models(m2).prettylabel, ...
                models(m0).prettyshortlabel, models(m1).prettyshortlabel, models(m2).prettyshortlabel, ...
                Ylabels, theseHorizons, tabcaption, statlabels, comparisonNote)


        end

    end % triplets
    %% finish wrap
    finishwrap

end % p

%% finish script
finishscript

%% helper function to create tex table
function compareTriple(tabname, wrap, doBold, ...
    maeloss0, maeloss1, relativeMAE01, maeloss2, relativeMAE02, ...
    crpsloss0, crpsloss1, relativeCRPS01, crpsloss2, relativeCRPS02, ...
    prettylabel0, prettylabel1, prettylabel2, ...
    ~, prettyshortlabel1, prettyshortlabel2, ...
    Ylabels, theseHorizons, tabcaption, statlabels, comparisonNote)


%% parse inputs
N         = length(Ylabels);
Nhorizons = length(theseHorizons);

%% DM tests

[relativeMAE01, dmMADtstat1] = dodm(maeloss0, maeloss1, relativeMAE01, theseHorizons);
[relativeCRPS01, dmCRPStstat1] = dodm(crpsloss0, crpsloss1, relativeCRPS01, theseHorizons);

[relativeMAE02, dmMADtstat2] = dodm(maeloss0, maeloss2, relativeMAE02, theseHorizons);
[relativeCRPS02, dmCRPStstat2] = dodm(crpsloss0, crpsloss2, relativeCRPS02, theseHorizons);

%% set up tab
if ~isempty(wrap)
    tabdir = wrap.dir;
    latexwrapper(wrap, 'add', 'sidetab', tabname, tabcaption)
end

%% tabulate
fid = fopen(fullfile(tabdir, tabname), 'wt');
fprintf(fid, '\\begin{center}\n');
fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('.2', 1, 4 * Nhorizons));
fprintf(fid, '\\toprule\n');
fprintf(fid, ' & \\multicolumn{%d}{c}{%s}  & \\multicolumn{%d}{c}{%s}   \\\\ \\cmidrule(lr){%d-%d}\\cmidrule(lr){%d-%d} \n', ...
    2 * Nhorizons, statlabels{1}, 2 * Nhorizons, statlabels{2}, ...
    1+1, 1+2*Nhorizons, 1+2*Nhorizons+1, 1+4*Nhorizons);
fprintf(fid, ' & \\multicolumn{%d}{c}{%s}  & \\multicolumn{%d}{c}{%s}  & \\multicolumn{%d}{c}{%s}  & \\multicolumn{%d}{c}{%s}  \\\\ \\cmidrule(lr){%d-%d}\\cmidrule(lr){%d-%d} \\cmidrule(lr){%d-%d}\\cmidrule(lr){%d-%d} \n', ...
    Nhorizons, prettyshortlabel1, Nhorizons, prettyshortlabel2, ...
    Nhorizons, prettyshortlabel1, Nhorizons, prettyshortlabel2, ...
    1+1,1+Nhorizons,1+Nhorizons+1,1+2*Nhorizons,1+2*Nhorizons+1,1+3*Nhorizons,1+3*Nhorizons+1,1+4*Nhorizons);


for h = 1 : Nhorizons
    fprintf(fid, '& \\multicolumn{1}{c}{$%d$} ', theseHorizons(h));
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
        if isfinite(relativeMAE01(n,h))
            if doBold && doColorCode(relativeMAE01(n,h))
                fprintf(fid, '& %s%s ', dcolbf(relativeMAE01(n,h), '%6.2f'), Zstar1(dmMADtstat1(n,h)));
            else
                fprintf(fid, '& %6.2f%s ', relativeMAE01(n,h), Zstar1(dmMADtstat1(n,h)));
            end
        else
            fprintf(fid, '& -.- ');
        end
    end
    for h = 1 : Nhorizons
        if isfinite(relativeMAE02(n,h))
            if doBold && doColorCode(relativeMAE02(n,h))
                fprintf(fid, '& %s%s ', dcolbf(relativeMAE02(n,h), '%6.2f'), Zstar1(dmMADtstat2(n,h)));
            else
                fprintf(fid, '& %6.2f%s ', relativeMAE02(n,h), Zstar1(dmMADtstat2(n,h)));
            end
        else
            fprintf(fid, '& -.- ');
        end
    end
    for h = 1 : Nhorizons
        if isfinite(relativeCRPS01(n,h))
            if doBold && doColorCode(relativeCRPS01(n,h))
                fprintf(fid, '& %s%s ', dcolbf(relativeCRPS01(n,h), '%6.2f'), Zstar1(dmCRPStstat1(n,h)));
            else
                fprintf(fid, '& %6.2f%s ', relativeCRPS01(n,h), Zstar1(dmCRPStstat1(n,h)));
            end
        else
            fprintf(fid, '& -.- ');
        end
    end
    for h = 1 : Nhorizons
        if isfinite(relativeCRPS02(n,h))
            if doBold && doColorCode(relativeCRPS02(n,h))
                fprintf(fid, '& %s%s ', dcolbf(relativeCRPS02(n,h), '%6.2f'), Zstar1(dmCRPStstat2(n,h)));
            else
                fprintf(fid, '& %6.2f%s ', relativeCRPS02(n,h), Zstar1(dmCRPStstat2(n,h)));
            end
        else
            fprintf(fid, '& -.- ');
        end
    end
    fprintf(fid, '\\\\\n');

end
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{center}\n');
fprintf(fid, '\n');

sigone = cat(1, abs(dmMADtstat1) > norminv(0.95, 0, 1) & (round(relativeMAE01,2) == 1), ...
    abs(dmCRPStstat1) > norminv(0.95, 0, 1) & (round(relativeCRPS01,2) == 1));

fprintf(fid, '\\legend{\n');
fprintf(fid, 'Comparison of %s and %s against %s (baseline, in denominator) for horizons', ...
    prettylabel1, prettylabel2, prettylabel0);
fprintf(fid, ' %d, ', theseHorizons(1:end-1));
fprintf(fid, 'and %d.\n', theseHorizons(end));
fprintf(fid, 'Values below 1 indicate improvement over baseline. \n');
fprintf(fid, '%s \n', comparisonNote);

fprintf(fid, 'Significance assessed by Diebold-Mariano-West test using Newey-West standard errors with $h + 1$ lags, and stars indicating $p$ values of 10\\%% and below.\n');

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
    fprintf(fid, 'Relative differences of 5 percent and more (compared to baseline) are indicated by bold face numbers.\n');
end

if ~all(isfinite(relativeMAE01(:)))
    fprintf(fid, 'In some cases, due to strong performance of the baseline model, relative MAE may involve divisions by zero. These cases are reported as blank entries.\n');
end
fprintf(fid, '}\n'); % close legend

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