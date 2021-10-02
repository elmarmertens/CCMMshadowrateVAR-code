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

%#ok<*NANMEAN>

%% setup

resultsdir = '~/jam/lager/quanticoELBmatfiles2021cum';

doCharts = false;

datalabel = 'fredMD20baa-2021-07';
% datalabel = 'fredMD20baaExYieldBAA-2021-07';


%% one big wrapper

titlename = sprintf('oosPairwiseEvaluationTables-%s', datalabel);
initwrap

%% eval window

for sam = [2 1]
    
    switch sam
        case 1
            % baseline
            evalStart  = datenum(2009,1,1);
            evalStop   = datenum(2021,6,1);
        case 2
            % ex COVID
            evalStart = datenum(2009,1,1);
            evalStop  = datenum(2017,12,1);
        case 3
            % COVID
            evalStart = datenum(2018,1,1);
            evalStop  = datenum(2021,6,1);
        case 4
            % middle
            evalStart  = datenum(2011,1,1);
            evalStop   = datenum(2017,12,1);
        otherwise
            error('sam %d not defined', sam)
    end
    
    
    evaltxt = sprintf('evalStart%sevalEnd%s', datestr(evalStart, 'yyyymm'), datestr(evalStop, 'yyyymm'));

    %% def models
    
    % STANDARD
    m = 1;
    models(m).datalabel    = datalabel; %#ok<*SAGROW>
    models(m).resultlabel  = 'standardVAR-tightBVARshrinkage-p12';
    models(m).prettylabel  = 'standard linear VAR';
    models(m).shortlabel   = 'Standard';
    models(m).fcstType     = 'fcstYcum';
    
    % SIMPLE SHADOW-RATE
    m = 2;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = 'simpleshadowrateVAR-tightBVARshrinkage-p12';
    models(m).prettylabel  = 'simple shadow-rate VAR';
    models(m).shortlabel   = 'ShadowRateVAR';
    models(m).fcstType     = 'fcstYcum';
    
    % BLOCK-HYBRID SHADOW RATE
    m = 3;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = 'hybridshadowrateVAR-tightBVARshrinkage-p12';
    models(m).prettylabel  = 'hybrid shadow-rate VAR';
    models(m).shortlabel   = 'BlockHybridVAR';
    models(m).fcstType     = 'fcstYcum';

    % KRIPPNER
    m = 4;
    models(m).datalabel    = 'fredMD20krippner-2021-07';
    models(m).resultlabel  = 'standardVAR-tightBVARshrinkage-p12';
    models(m).prettylabel  = 'Plug-in VAR (Krippner)';
    models(m).shortlabel   = 'krippnerVAR';
    models(m).fcstType     = 'fcstYcumshadow';
    
    % WUXIA
    m = 5;
    models(m).datalabel    = 'fredMD20wuxia-2021-07';
    models(m).resultlabel  = 'standardVAR-tightBVARshrinkage-p12';
    models(m).prettylabel  = 'plug-in VAR (WuXia)';
    models(m).shortlabel   = 'wuxiaVAR';
    models(m).fcstType     = 'fcstYcumshadow';
    
    
    %% define model sets
    
    if strcmpi(datalabel, 'fredMD20baa-2021-07')
        MODELS  = {[1 2] ...
            [1 3] [2 3] ...
            [4 2] [5 2] ...
            [4 3] [5 3] ...
            };       
        MLABELS = {'Linear vs Simple shadow-rate VAR', ...
            'Linear vs Hybrid shadow-rate VAR', ...
            'Simple vs Hybrid shadow-rate VAR', ...
            'Plug-in Krippner vs Simple shadow-rate VAR', 'Plug-in Wu-Xia vs Simple shadow-rate VAR', ...
            'Plug-in Krippner vs Hybrid shadow-rate VAR', 'Plug-in Wu-Xia vs Hybird shadow-rate VAR', ...
            }; 
    else
        MODELS  = {[1 2] ...
            [1 3] [2 3] ...
            };       
        MLABELS = {'Linear vs Simple shadow-rate VAR', ...
            'Linear vs Hybrid shadow-rate VAR', ...
            'Simple vs Hybrid shadow-rate VAR', ...
            }; 
    end
    
    %% loop over model sets
    
    for m = 1 : length(MODELS) % mset % [2 8 9] % 2 : 4 % 10 % length(MODELS)
        
        
        
        m0 = MODELS{m}(1);
        m1 = MODELS{m}(2);
        
        
        doQuarterly = false;
        
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
        
        comparisonNote = sprintf('Evaluation window from %s through %s.', ...
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
        compareAllinone(tabname, wrap, doCharts, ...
            mseloss0, mseloss1, relativeRMSE01, ...
            maeloss0, maeloss1, relativeMAD01, ...
            crpsloss0, crpsloss1, relativeCRPS01, ...
            models(m0).prettylabel, models(m1).prettylabel, ...
            Ylabels, theseHorizons, tabcaption, statlabels, comparisonNote)
        
        
    end
end

%% finish script
finishwrap
finishscript


function compareAllinone(tabname, wrap, doCharts, ...
    mseloss0, mseloss1, relativeRMSE01, ...
    maeloss0, maeloss1, relativeMAE01, ...
    crpsloss0, crpsloss1, relativeCRPS01, ...
    prettylabel0, prettylabel1, ...
    Ylabels, theseHorizons, tabcaption, statlabels, comparisonNote)


%% parse inputs
N = length(Ylabels);
Nhorizons = length(theseHorizons);

%% DM tests

[relativeRMSE01, dmRMSEtstat] = dodm(mseloss0, mseloss1, relativeRMSE01, theseHorizons);

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
            switch doBold(relativeRMSE01(n,h))
                case 1
                    fprintf(fid, '& %s ', dcolred(sprintf('%6.2f%s ', relativeRMSE01(n,h), Zstar(dmRMSEtstat(n,h)))));
                case -1
                    fprintf(fid, '& %s ', dcolgreen(sprintf('%6.2f%s ', relativeRMSE01(n,h), Zstar(dmRMSEtstat(n,h)))));
                otherwise
                    fprintf(fid, '& %6.2f%s ', relativeRMSE01(n,h), Zstar(dmRMSEtstat(n,h)));
            end
        else
            fprintf(fid, '& -.- ');
        end
    end
    for h = 1 : Nhorizons
        if isfinite(relativeMAE01(n,h))
            switch doBold(relativeMAE01(n,h))
                case 1
                    fprintf(fid, '& %s ', dcolred(sprintf('%6.2f%s ', relativeMAE01(n,h), Zstar(dmMADtstat(n,h)))));
                case -1
                    fprintf(fid, '& %s ', dcolgreen(sprintf('%6.2f%s ', relativeMAE01(n,h), Zstar(dmMADtstat(n,h)))));
                otherwise
                    fprintf(fid, '& %6.2f%s ', relativeMAE01(n,h), Zstar(dmMADtstat(n,h)));
            end
        else
            fprintf(fid, '& -.- ');
        end
    end
    for h = 1 : Nhorizons
        if isfinite(relativeCRPS01(n,h))
            switch doBold(relativeCRPS01(n,h))
                case 1
                    fprintf(fid, '& %s ', dcolred(sprintf('%6.2f%s ', relativeCRPS01(n,h), Zstar(dmCRPStstat(n,h)))));
                case -1
                    fprintf(fid, '& %s ', dcolgreen(sprintf('%6.2f%s ', relativeCRPS01(n,h), Zstar(dmCRPStstat(n,h)))));
                otherwise
                    fprintf(fid, '& %6.2f%s ', relativeCRPS01(n,h), Zstar(dmCRPStstat(n,h)));
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
    sigone = cat(1, abs(dmRMSEtstat) > norminv(0.95, 0, 1) & (round(relativeRMSE01,2) == 1), ...
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
    
    
    if ~all(isfinite(relativeMAE01(:)))
        fprintf(fid, 'In some cases, due to strong performance of the baseline model, relative MAD may involve divisions by zero. These cases are reported as blank entries.');
    end
end
fclose(fid);

type(fullfile(tabdir, tabname))

end

function flag = doBold(x)

if x > 1.05
    flag = 1;
elseif x < .95
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