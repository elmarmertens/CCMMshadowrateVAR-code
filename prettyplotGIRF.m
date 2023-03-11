%#ok<*NOSEL>
%#ok<*DISPLAYPROG>
%#ok<*UNRCH>
%#ok<*ASGLU>
%#ok<*DATNM>
%#ok<*DATST>

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
DATALABEL           = {'fredMD20VXO-2022-09', 'fredMD20VXOexYield-2022-09','fredMD20EBP-2022-09', 'fredMD20EBPexYield-2022-09'};
modeltype           = 'BlockHybrid';
modelpretty         = 'Hybrid Shadow-Rate VAR';
resultsdir          = 'irfVXOEBP';

resultsdirALT       = 'irfVXOEBP';
modeltypeALT        = 'Linear';
modelprettyALT      = 'Linear VAR';

p                   = 12;                    % Number of lags on dependent variables
irfHorizon          = 25;

samStart            = [];

irfDate0            = datenum(2007,1,1);
irfDATES            = [datenum(2009,1,1) datenum([2010 2012 2014],12,1)]; %#ok<NASGU> 
irfSCALES           = [1 10]; %#ok<NASGU> 

% SED-PARAMETERS-HERE
irfDATES            = datenum(2012,12,1);
irfSCALES           = 1;

quicky = false;

np = 12;

fontsize = 24;

%% GIRF

for dd = 1  :  length(DATALABEL)

    datalabel = DATALABEL{dd};
    if contains(datalabel, 'PreCovid')
        jumpDate            = datenum(2019,12,01);
    else
        jumpDate            = datenum(2022,08,01);
    end

    for IRF1scale = irfSCALES

        display(IRF1scale)

        for irfDate = irfDATES

            display(datestr(irfDate, 'yyyymmm'))


            % prepare wrap

            titlename=sprintf('prettyGIRF-%s-%s-p%d-IRF1scale%d-jumpoff%s-irfDate%s-vs-%s', modeltype, datalabel, p, IRF1scale, ...
                datestr(jumpDate, 'yyyymmm'), datestr(irfDate0, 'yyyymmm'), datestr(irfDate, 'yyyymmm'));


            if ~isempty(samStart)
                titlename = strcat(titlename,'-', datestr(samStart, 'yyyymmm'));
            end

            if quicky
                wrap = [];
            else
                initwrap
            end

            %% load IRF
            % load 0
            filename=sprintf('irf%s-%s-p%d-IRF1scale%d-jumpoff%s-irfDate%s', modeltype, ...
                datalabel, p, IRF1scale, ...
                datestr(jumpDate, 'yyyymmm'), datestr(irfDate0, 'yyyymmm'));
            irf0 = matfile(fullfile(resultsdir, filename));

            % load 1
            filename=sprintf('irf%s-%s-p%d-IRF1scale%d-jumpoff%s-irfDate%s', modeltype, ...
                datalabel, p, IRF1scale, ...
                datestr(jumpDate, 'yyyymmm'), datestr(irfDate, 'yyyymmm'));
            irf1 = matfile(fullfile(resultsdir, filename));

           % load ALT
            filename=sprintf('irf%s-%s-p%d-IRF1scale%d-jumpoff%s-irfDate%s', modeltypeALT, ...
                datalabel, p, IRF1scale, ...
                datestr(jumpDate, 'yyyymmm'), datestr(irfDate, 'yyyymmm'));
            irfALT = matfile(fullfile(resultsdirALT, filename));

            ncode = irf0.ncode;
            if ~isequal(ncode, irf1.ncode)
                error('data set mismatch')
            end
            if ~isequal(ncode, irfALT.ncode)
                error('data set mismatch')
            end
            Ylabels = irf0.Ylabels;
            N = length(ncode);

            %% PLOT RESULTS
            color0     = Colors4Plots(1);
            color1     = Colors4Plots(2);
            color2     = Colors4Plots(1);
            color3     = Colors4Plots(3);
            colorALT   = Colors4Plots(8);

            %% plot ELB IRF
            for n = 1 : N

                % FIGURE TYPE 1
                % hybrid at different dates
                thisfig1 = figure;
                hold on
                ax1 = gca;
                set(ax1, 'FontSize', fontsize)
                h1 = plot(1:irfHorizon-1, irf1.IRF1plus(n,1:end-1), '-', 'color', color1, 'linewidth', 4);
                plot(1:irfHorizon-1, squeeze(irf1.IRF1plusTails(n,1:end-1,:,:)), '-', 'color', color1, 'linewidth', 2);
                h0  = plot(1:irfHorizon-1, irf0.IRF1plus(n,1:end-1), '-.', 'color', color0, 'linewidth', 3);             
                plot(1:irfHorizon-1, squeeze(irf0.IRF1plusTails(n,1:end-1,:,:)), '-.', 'color', color0, 'linewidth', 1);
                xlim([1 irfHorizon-1])
                xticks([1, 6:6:irfHorizon-1])
                grid on
                YLIM1 = ylim;
               
                
                % FIGURE TYPE 2
                % baseline vs ALT
                thisfig2 = figure;
                hold on
                ax2 = gca;
                set(ax2, 'FontSize', fontsize)
                h2BASE = plot(1:irfHorizon-1, irf1.IRF1plus(n,1:end-1), '-', 'color', color1, 'linewidth', 4);
                plot(1:irfHorizon-1, squeeze(irf1.IRF1plusTails(n,1:end-1,:,:)), '-', 'color', color1, 'linewidth', 2);
                h2ALT  = plot(1:irfHorizon-1, irfALT.IRF1plus(n,1:end-1), '-.', 'color', colorALT, 'linewidth', 3);             
                plot(1:irfHorizon-1, squeeze(irfALT.IRF1plusTails(n,1:end-1,:,:)), '-.', 'color', colorALT, 'linewidth', 1);
                xlim([1 irfHorizon-1])
                xticks([1, 6:6:irfHorizon-1])
                % yline(0, 'k:')
                grid on
                YLIM2 = ylim;

                % align axis limits
                YLIM = [min([YLIM1(1), YLIM2(1)]), max([YLIM1(2), YLIM2(2)])];
                ylim(ax1, YLIM);
                ylim(ax2, YLIM);
                
                % enforce common yticks
                yticks(ax2, yticks(ax1))
                
                % STORE FIGURES
                wrapthisfigure(thisfig1, sprintf('IRF1scale%d-%s-%s-%s-jumpoff%s-irfDate%s-vs-%s', IRF1scale, modeltype, datalabel, ncode{n}, ...
                    datestr(jumpDate, 'yyyymmm'), datestr(irfDate0, 'yyyymmm'), datestr(irfDate, 'yyyymmm')), wrap, [], [], [], [], true);
                legend([h0, h1], datestr(irfDate0, 'yyyy mmm'), datestr(irfDate, 'yyyy mmm'), 'location', 'best')
                wrapthisfigure(thisfig1, sprintf('IRF1scale%d-%s-%s-%s-jumpoff%s-irfDate%s-vs-%s-WITHLEGEND', IRF1scale, modeltype, datalabel, ncode{n}, ...
                    datestr(jumpDate, 'yyyymmm'), datestr(irfDate0, 'yyyymmm'), datestr(irfDate, 'yyyymmm')), wrap)
     
               
                wrapthisfigure(thisfig2, sprintf('IRF1scale%d-%s-vs-%s-%s-%s-jumpoff%s-irfDate%s', IRF1scale, modeltype, modeltypeALT, datalabel, ncode{n}, ...
                    datestr(jumpDate, 'yyyymmm'), datestr(irfDate, 'yyyymmm')), wrap, [], [], [], [], true)
                legend([h2BASE, h2ALT], modelpretty, modelprettyALT, 'location', 'best')
                wrapthisfigure(thisfig2, sprintf('IRF1scale%d-%s-vs-%s-%s-%s-jumpoff%s-irfDate%s-WITHLEGEND', IRF1scale, modeltype, modeltypeALT, datalabel, ncode{n}, ...
                    datestr(jumpDate, 'yyyymmm'), datestr(irfDate, 'yyyymmm')), wrap)
     
                if ~quicky
                    close(thisfig1)
                    close(thisfig2)
                end

            end

            %% wrap up
            if quicky
                dockAllFigures
                return
            else
                close all
                finishwrap
            end

        end % irfDate
    end % irfscale
end % datalabel

dockAllFigures
finishscript
