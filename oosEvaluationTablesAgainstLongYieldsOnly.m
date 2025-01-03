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

doBold   = true;

p = 12;

model0.datalabel        = 'fredsxMD14longyields-2022-09';
model0.resultlabel      = 'standardVARAR1SV-RATSbvarshrinkage-p12';
model0.prettylabel      = 'Linear VAR without short-term yields';
model0.prettyshortlabel = 'Linear';
model0.shortlabel       = 'standardVAR-p12';
model0.fcstType         = 'fcstY';

model1.datalabel    = 'fredsxMD20exYield-2022-09';
model1.resultlabel  = sprintf('ELBnonstructuralAR1SV-RATSbvarshrinkage-p%d', p);
model1.prettylabel  = sprintf('Non-structural shadow-rate VAR (w/o yields)');
model1.prettyshortlabel  = sprintf('w/o yields');
model1.shortlabel   = sprintf('nonstructuralVAR-p%d', p);
model1.fcstType     = 'fcstY';

model2.datalabel    = 'fredsxMD20-2022-09';
model2.resultlabel  = sprintf('ELBblocknonstructuralAR1SV-RATSbvarshrinkage-p%d', p);
model2.prettylabel  = sprintf('Restricted non-structural shadow-rate VAR (w/yields)');
model2.prettyshortlabel  = sprintf('w/yields');
model2.shortlabel   = sprintf('blocknonstructuralVAR-p%d', p);
model2.fcstType     = 'fcstY';

doFlipOrder = true;

% NOTE: model2 is assumed to comprise all variables


%% one wrapper per lag choice

titlename = 'oosEvaluationTablesLongyieldsOnlyQE';
initwrap



thisTRIPLETlabel = sprintf('%s and %s', model1.prettylabel, model2.prettylabel);

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

    oos0 = matfile(fullfile(resultsdir, sprintf('%s-%s.mat', model0.datalabel, model0.resultlabel)));
    oos1 = matfile(fullfile(resultsdir, sprintf('%s-%s.mat', model1.datalabel, model1.resultlabel)));
    oos2 = matfile(fullfile(resultsdir, sprintf('%s-%s.mat', model2.datalabel, model2.resultlabel)));

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

    %% find common set of variables
    ncode0   = oos0.ncode;
    ncode1   = oos1.ncode;
    ncode2   = oos2.ncode;   

    ncode     = ncode2; % ncode2 is assumed to be the all-encompassing list of variables
    
    Ylabels   = fredMDveryshortlabel(ncode);
    Ylabels   = strrep(Ylabels, '_', '');
    N         = length(Ylabels);
    ndxY0     = ismember(ncode, ncode0);
    ndxY1     = ismember(ncode, ncode1);
    ndxY2     = ismember(ncode, ncode2);
    
    if ~isequal(ncode(ndxY0), ncode0)
        error houston0
    end
    if ~isequal(ncode(ndxY1), ncode1)
        error houston1
    end
    if ~isequal(ncode(ndxY2), ncode2)
        error houston2
    end
    
    %% setup monthly tables
    theseHorizons = [6 12 24];

    %% MAE
    losstype0 = sprintf('%smederror', model0.fcstType);
    losstype1 = sprintf('%smederror', model1.fcstType);
    losstype2 = sprintf('%smederror', model2.fcstType);


    [maeloss0, maeloss1, maeloss2] = deal(NaN(N, Nhorizons, length(ndxJumpoff)));
    maeloss0(ndxY0,:,:) = abs(oos0.(losstype0));
    maeloss1(ndxY1,:,:) = abs(oos1.(losstype1));
    maeloss2(ndxY2,:,:) = abs(oos2.(losstype2));
    % cut sample
    maeloss0 = maeloss0(:,:,ndxJumpoff);
    maeloss1 = maeloss1(:,:,ndxJumpoff);
    maeloss2 = maeloss2(:,:,ndxJumpoff);
    % averages
    mae0          = nanmean(maeloss0,3);
    mae1          = nanmean(maeloss1,3);
    mae2          = nanmean(maeloss2,3);
    % relatives
    relativeMAD01 = mae1(:,1:Nhorizons) ./ mae0(:,1:Nhorizons);
    relativeMAD02 = mae2(:,1:Nhorizons) ./ mae0(:,1:Nhorizons);
    

    %% CRPS
    losstype0     = sprintf('%scrps', model0.fcstType);
    losstype1     = sprintf('%scrps', model1.fcstType);
    losstype2     = sprintf('%scrps', model2.fcstType);

    [crpsloss0, crpsloss1, crpsloss2] = deal(NaN(N, Nhorizons, length(ndxJumpoff)));
    crpsloss0(ndxY0,:,:) = oos0.(losstype0);
    crpsloss1(ndxY1,:,:) = oos1.(losstype1);
    crpsloss2(ndxY2,:,:) = oos2.(losstype2);
    % cut sample
    crpsloss0 = crpsloss0(:,:,ndxJumpoff);
    crpsloss1 = crpsloss1(:,:,ndxJumpoff);
    crpsloss2 = crpsloss2(:,:,ndxJumpoff);
    % averages
    crps0          = nanmean(crpsloss0,3);
    crps1          = nanmean(crpsloss1,3);
    crps2          = nanmean(crpsloss2,3);
    % relatives
    relativeCRPS01 = crps1(:,1:Nhorizons) ./ crps0(:,1:Nhorizons);
    relativeCRPS02 = crps2(:,1:Nhorizons) ./ crps0(:,1:Nhorizons);
    
    %% prune variables with all NaNs
    nanny = all(isnan(relativeMAD01),2) & all(isnan(relativeMAD02),2) & all(isnan(relativeCRPS01),2) & all(isnan(relativeCRPS02),2);
    if any(nanny)
        fprintf('Dropping %d variables with NaNs\n', sum(nanny))
        Ylabels = Ylabels(~nanny);
        relativeMAD01 = relativeMAD01(~nanny,:);
        relativeMAD02 = relativeMAD02(~nanny,:);
        maeloss0 = maeloss0(~nanny,:,:);
        maeloss1 = maeloss1(~nanny,:,:);
        maeloss2 = maeloss2(~nanny,:,:);
        relativeCRPS01 = relativeCRPS01(~nanny,:);
        relativeCRPS02 = relativeCRPS02(~nanny,:);
        crpsloss0 = crpsloss0(~nanny,:,:);
        crpsloss1 = crpsloss1(~nanny,:,:);
        crpsloss2 = crpsloss2(~nanny,:,:);
    end
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
    tabname = sprintf('tripleComparisonQE-%s-%s-vs-%s-vs-%s-%s.tex', model0.datalabel, ...
        model0.shortlabel, model1.shortlabel, model2.shortlabel, evaltxt);
    tabcaption =  sprintf('%s %s', thisTRIPLETlabel, shortcomparisonNote);

    compareTriple(tabname, wrap, doBold, ...
        maeloss0, maeloss1, relativeMAD01, maeloss2, relativeMAD02, ...
        crpsloss0, crpsloss1, relativeCRPS01, crpsloss2, relativeCRPS02, ...
        model0.prettylabel, model1.prettylabel, model2.prettylabel, ...
        model0.prettyshortlabel, model1.prettyshortlabel, model2.prettyshortlabel, ...
        Ylabels, theseHorizons, tabcaption, statlabels, comparisonNote)


end

%% finish wrap
finishwrap


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
fprintf(fid, 'Comparison of ``%s'''' and ``%s'''' against ``%s'' (baseline, in denominator) for horizons', ...
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