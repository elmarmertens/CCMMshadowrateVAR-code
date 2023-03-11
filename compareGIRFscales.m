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

irfDATES            = [datenum(2009,1,1) datenum([2010 2012 2014],12,1)];

irfScale0           = 1;
irfScale1           = 5;
irfScale2           = 10;


% SED-PARAMETERS-HERE
irfDATES            = [datenum(2007,1,1) datenum([2010 2014],12,1)];


np = 12;

fontsize = 18;

%% GIRF

for dd = 1 %  1 :  length(DATALABEL)

    datalabel = DATALABEL{dd};
    jumpDate  = datenum(2022,08,01);



    for irfDate = irfDATES

        display(datestr(irfDate, 'yyyymmm'))


        % prepare wrap

        titlename=sprintf('CompareIRFscales-%s-%s-p%d-jumpoff%s-irfDate%s', modeltype, datalabel, p, ...
            datestr(jumpDate, 'yyyymmm'), datestr(irfDate, 'yyyymmm'));


        if ~isempty(samStart)
            titlename = strcat(titlename,'-', datestr(samStart, 'yyyymmm'));
        end

        wrap = [];
        initwrap

        %% load IRF
        % load 0
        filename=sprintf('irf%s-%s-p%d-IRF1scale%d-jumpoff%s-irfDate%s', modeltype, ...
            datalabel, p, irfScale0, ...
            datestr(jumpDate, 'yyyymmm'), datestr(irfDate, 'yyyymmm'));
        irf0 = matfile(fullfile(resultsdir, filename));

        % load 1
        filename=sprintf('irf%s-%s-p%d-IRF1scale%d-jumpoff%s-irfDate%s', modeltype, ...
            datalabel, p, irfScale1, ...
            datestr(jumpDate, 'yyyymmm'), datestr(irfDate, 'yyyymmm'));
        irf1 = matfile(fullfile(resultsdir, filename));

        % load 2
        filename=sprintf('irf%s-%s-p%d-IRF1scale%d-jumpoff%s-irfDate%s', modeltype, ...
            datalabel, p, irfScale2, ...
            datestr(jumpDate, 'yyyymmm'), datestr(irfDate, 'yyyymmm'));
        irf2 = matfile(fullfile(resultsdir, filename));

        ncode = irf0.ncode;
        if ~isequal(ncode, irf1.ncode)
            error('data set mismatch')
        end
        if ~isequal(ncode, irf2.ncode)
            error('data set mismatch')
        end
        Ylabels = irf0.Ylabels;
        N = length(ncode);

        %% PLOT RESULTS
        color0     = Colors4Plots(1);
        color1     = Colors4Plots(2);
        color2     = Colors4Plots(8);
        %         colorBase     = Colors4Plots(8);

        %% plot ELB IRF
        for n = 1 : N

            thisfig = figure;
            subplot(2,1,1)
            hold on
            set(gca, 'FontSize', fontsize)
            h0  = plot(1:irfHorizon-1, irf0.IRF1plus(n,1:end-1), '-', 'color', color0, 'linewidth', 4);             %#ok<NASGU> 
            h1  = plot(1:irfHorizon-1, irf1.IRF1plus(n,1:end-1) / irfScale1 * irfScale0, '--', 'color', color1, 'linewidth', 3);             %#ok<NASGU> 
            h2  = plot(1:irfHorizon-1, irf2.IRF1plus(n,1:end-1) / irfScale2 * irfScale0, 'd', 'color', color2, 'linewidth', 3);             %#ok<NASGU> 
            xlim([0 irfHorizon-1])
            xticks(0:6:irfHorizon-1)
            yline(0, 'k:')
            title('positive shock', 'FontWeight', 'normal')

            subplot(2,1,2)
            hold on
            set(gca, 'FontSize', fontsize)
            h0  = plot(1:irfHorizon-1, irf0.IRF1minus(n,1:end-1), '-', 'color', color0, 'linewidth', 4);             
            h1  = plot(1:irfHorizon-1, irf1.IRF1minus(n,1:end-1) / irfScale1 * irfScale0, '--', 'color', color1, 'linewidth', 3);             
            h2  = plot(1:irfHorizon-1, irf2.IRF1minus(n,1:end-1) / irfScale2 * irfScale0, 'd', 'color', color2, 'linewidth', 3);             
            xlim([0 irfHorizon-1])
            xticks(0:6:irfHorizon-1)
            yline(0, 'k:')
            title('negative shock', 'FontWeight', 'normal')

            wrapthisfigure(thisfig, sprintf('CompareScales-%s-%s-%s-jumpoff%s-irfDate%s', modeltype, datalabel, ncode{n}, ...
                datestr(jumpDate, 'yyyymmm'), datestr(irfDate, 'yyyymmm')), wrap, [], [], [], [], true);
            hl = legend([h0, h1, h2], sprintf('scale %d', irfScale0), sprintf('scale %d (and rescaled)', irfScale1), ...
                sprintf('scale %d (and rescaled)', irfScale2), 'location', 'best', 'box', 'off');
            hl.NumColumns = 1;
            wrapthisfigure(thisfig, sprintf('CompareScales-%s-%s-%s-jumpoff%s-irfDate%s-WITHLEGEND', modeltype, datalabel, ncode{n}, ...
                datestr(jumpDate, 'yyyymmm'), datestr(irfDate, 'yyyymmm')), wrap, [], [], [], [], true);
            sgtitle(sprintf('%s\n', Ylabels{n}), 'FontSize', 18', 'FontWeight', 'bold')
            wrapthisfigure(thisfig, sprintf('CompareScales-%s-%s-%s-jumpoff%s-irfDate%s-WITHTITLE', modeltype, datalabel, ncode{n}, ...
                datestr(jumpDate, 'yyyymmm'), datestr(irfDate, 'yyyymmm')), wrap);

        end


        %% wrap up
        close all
        finishwrap

    end % irfDate
end % datalabel

dockAllFigures
finishscript
