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
resultsdir          = 'irfVXOEBP';

p                   = 12;                    % Number of lags on dependent variables
irfHorizon          = 25;

samStart            = [];

irfDate0            = datenum(2007,1,1);
irfDATES            = [datenum(2009,1,1) datenum([2010 2012 2014],12,1)];
irfSCALES           = [1 10];

% SED-PARAMETERS-HERE


np = 12;

fontsize = 12;

%% GIRF

for dd = 3 :  length(DATALABEL)

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

            titlename=sprintf('%s-%s-p%d-IRF1scale%d-jumpoff%s-irfDate%s-vs-%s', modeltype, datalabel, p, IRF1scale, ...
                datestr(jumpDate, 'yyyymmm'), datestr(irfDate0, 'yyyymmm'), datestr(irfDate, 'yyyymmm'));


            if ~isempty(samStart)
                titlename = strcat(titlename,'-', datestr(samStart, 'yyyymmm'));
            end

            wrap = [];
            initwrap

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

            ncode = irf0.ncode;
            if ~isequal(ncode, irf1.ncode)
                error('data set mismatch')
            end
            Ylabels = irf0.Ylabels;
            N = length(ncode);

            %% PLOT RESULTS
            color0     = Colors4Plots(1);
            color1     = Colors4Plots(2);
            %         colorBase     = Colors4Plots(8);

            %% plot ELB IRF
            for n = 1 : N

                thisfig = figure;
                subplot(2,1,1)
                hold on
                set(gca, 'FontSize', fontsize)
                h0  = plot(0:irfHorizon-1, irf0.IRF1plus(n,:), '-', 'color', color0, 'linewidth', 3);             plot(0:irfHorizon-1, squeeze(irf0.IRF1plusTails(n,:,:,:)), '-', 'color', color0, 'linewidth', 1);
                h1 = plot(0:irfHorizon-1, irf1.IRF1plus(n,:), '-.', 'color', color1, 'linewidth', 3);
                plot(0:irfHorizon-1, squeeze(irf1.IRF1plusTails(n,:,:,:)), '-.', 'color', color1, 'linewidth', 1);
                xlim([0 irfHorizon-1])
                xticks(0:6:irfHorizon-1)
                yline(0, 'k:')
                legend([h0, h1], datestr(irfDate0, 'yyyy mmm'), datestr(irfDate, 'yyyy mmm'), 'location', 'best')
                title('positive shock', 'FontWeight', 'normal')

                subplot(2,1,2)
                hold on
                set(gca, 'FontSize', fontsize)
                h0  = plot(0:irfHorizon-1, irf0.IRF1minus(n,:), '-', 'color', color0, 'linewidth', 3);
                plot(0:irfHorizon-1, squeeze(irf0.IRF1minusTails(n,:,:,:)), '-', 'color', color0, 'linewidth', 1);
                h1 = plot(0:irfHorizon-1, irf1.IRF1minus(n,:), '-.', 'color', color1, 'linewidth', 3);
                plot(0:irfHorizon-1, squeeze(irf1.IRF1minusTails(n,:,:,:)), '-.', 'color', color1, 'linewidth', 1);
                xlim([0 irfHorizon-1])
                xticks(0:6:irfHorizon-1)
                yline(0, 'k:')
                %             legend([h0, h1], datestr(irfDate0, 'yyyy mmm'), datestr(irfDate, 'yyyy mmm'), 'location', 'best')
                title('negative shock', 'FontWeight', 'normal')

                sgtitle(sprintf('%s\n(shocks of size %d)', Ylabels{n}, IRF1scale), 'FontSize', 18', 'FontWeight', 'bold')

                wrapthisfigure(thisfig, sprintf('IRF1plusminus-%s-%s-IRF1scale%d-%s-jumpoff%s-irfDate%s-vs-%s', modeltype, datalabel, IRF1scale, ncode{n}, ...
                    datestr(jumpDate, 'yyyymmm'), datestr(irfDate0, 'yyyymmm'), datestr(irfDate, 'yyyymmm')), wrap)
            end


            %% wrap up
            close all
            finishwrap
        end % irfDate
    end % irfscale
end % datalabel

dockAllFigures
finishscript
