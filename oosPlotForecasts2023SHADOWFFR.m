%% compare OOS results from two estimates
% load quantico*.mat files and assess OOS performance

clear
close all
fclose all;

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

%% setup

resultsdir = '../matfilesShadowrateVAR/lagerFREDblock';

doDateTicks        = true;
doLinearShadowOnly = true;
doLinearActualOnly = true;
doShadowActualOnly = true;
doShadowOnly       = true;
doLinearOnly       = true;

doActualMean       = true;

Nhorizons = 24;

for DATALABEL = {'fredblockMD20-2022-09'} % {'fredMD20-2022-09', 'fredMD20exYield-2022-09'}

    datalabel = DATALABEL{:};

    %% list of models
    % STANDARD
    m = 1;
    models(m).datalabel    = datalabel; %#ok<*SAGROW>
    models(m).resultlabel  = 'standardVAR-RATSbvarshrinkage-p12';
    models(m).prettylabel  = 'Standard';
    models(m).shortlabel   = 'Standard';
    models(m).fcstType     = 'fcstY';

    m = 2;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = 'ELBblockhybrid-RATSbvarshrinkage-p12';
    models(m).prettylabel  = 'Hybrid sadow-rate';
    models(m).shortlabel   = 'HybridShadowRateVAR';
    models(m).fcstType     = 'fcstY';

    % SIMPLE SHADOWRATE (SY)
    m = 3;
    models(m).datalabel    = datalabel;
    models(m).resultlabel  = 'ELBsampling-RATSbvarshrinkage-p12';
    models(m).prettylabel  = 'Simple shadow-rate';
    models(m).shortlabel   = 'SimpleShadowRateVAR';
    models(m).fcstType     = 'fcstY';


    %% loop over hybrid and simple shadow-rateVAR
    for doHybrid = true % [false true]

        %% parameters

        m1  = 1;

        if doHybrid
            m0  = 2;
            shadowshortlabel = 'hybridshadowrate';
        else
            m0  = 3;
            shadowshortlabel = 'simpleshadowrate';
        end

        fontsize = 16;

        if isdesktop
            wrap = [];
        else
            titlename = sprintf('oosPlotForecasts2023chartsSHADOWFFR-%s-%s', datalabel, shadowshortlabel);
            initwrap
        end

        %#ok<*UNRCH>
        %% load data

        mat0 = matfile(fullfile(resultsdir, sprintf('%s-%s.mat', models(m0).datalabel, models(m0).resultlabel)));
        mat1 = matfile(fullfile(resultsdir, sprintf('%s-%s.mat', models(m1).datalabel, models(m1).resultlabel)));

        ydates    = mat0.ydates;

        Tjumpoffs = mat0.Tjumpoffs;

        %% cut eval sample if desired
        ndxJumpoff = true(size(Tjumpoffs)); % legacy variable

        Tjumpoffs   = Tjumpoffs(ndxJumpoff);


        dates    = ydates(Tjumpoffs);

        %% some parameters
        % Nhorizons  = mat0.fcstNhorizons;

        ncode = mat0.ncode;

        Ylabels = fredMDprettylabel(ncode);
        N       = length(Ylabels);

        Ylabels = strrep(Ylabels, '_', '');


        %% choose dates to plot
        % if isdesktop
        %     datechoice = datenum(2020,3,1);
        % else
        %     % datechoice = sort([datenum(2021,1:6,1) datenum(2020,1:12,1) datenum(2012:2019,10,1), datenum(2009:2019,12,1), datenum(2009:2019,1,1), datenum(2009:2019,6,1)]);
        %     datechoice = sort([datenum(2022,1:8,1) datenum(2021,1:12,1) datenum(2020,1:12,1) datenum(2012:2019,10,1), datenum(2009:2019,12,1), datenum(2009:2019,1,1), datenum(2009:2019,6,1)]);
        % end

        datechoice = [datenum(2021,6,1), datenum(2020,3,1), datenum(2020,9,1), datenum(2022,3,1), datenum(2022,8,1)];
        datechoice = datenum(2021,8,1);

        %% loop over forecast origins
        for thisdate = datechoice

            thisT = find(dates == thisdate);
            thisY = find(ismember(ncode, 'FEDFUNDS'));

            if doDateTicks
                [thisYear, thisMonth] = datevec(dates(thisT));
                forecastDates = datenum(thisYear, thisMonth + (1 : Nhorizons), 1);
            else
                forecastDates = 1:Nhorizons;
            end

            %% realized
            realized = squeeze(mat0.fcstYrealized(thisY,:,thisT)); % actuals from shadow-rate model are censored at ELB

            %% load draws

            ndxCI = find(ismember(mat0.setQuantiles, normcdf([-1 1]) * 100));


            fcstYmedian0 = mat0.(sprintf('%smedian', models(m0).fcstType));
            med0 = fcstYmedian0(thisY,1:Nhorizons,thisT);
            fcstYhat0 = mat0.(sprintf('%shat', models(m0).fcstType));
            mid0 = fcstYhat0(thisY,1:Nhorizons,thisT);

            fcstYquantiles0 = mat0.(sprintf('%squantiles', models(m0).fcstType));
            tail0 = squeeze(fcstYquantiles0(thisY,1:Nhorizons,ndxCI,thisT));

            n = 1;

            shadowmid   = mat0.fcstShadowYmedian(n,1:Nhorizons,thisT);
            these       = mat0.fcstShadowYquantiles;
            shadowtails = squeeze(these(n,1:Nhorizons,ndxCI,thisT));

            linearmid   = mat1.fcstYmedian(thisY,1:Nhorizons,thisT);
            these       = mat1.fcstYquantiles;
            lineartails = squeeze(these(thisY,1:Nhorizons,ndxCI,thisT));

            %% plot ShadowActualOnly
            if doShadowActualOnly
                thisfig = figure;
                ax = gca;
                set(ax, 'fontsize', fontsize)

                hold on

                plotCIaltcolor(shadowmid, shadowtails, forecastDates, [], 'w--', 'linewidth', 3);

                hs = plot(forecastDates, shadowmid, '-',  'color', [0.5843    0.8157    0.9882], 'linewidth', 4);
                plot(forecastDates, shadowmid, 'w--',  'linewidth', 4);


                %     hrealized = plot(forecastDates, realized, 'd', 'color', [0 .5 0], 'linewidth', 3);

                h1 = plot(forecastDates, med0, 'b-', 'linewidth', 4);
                h1tail = plot(forecastDates, tail0, 'b-', 'linewidth', 2); %#ok<NASGU>
                if doActualMean
                    h1hat = plot(forecastDates, mid0, 'bd', 'linewidth', 4);
                end

                if doDateTicks
                    xticks(forecastDates([1, 6 : 6 : end]))
                    xlim(forecastDates([1 end]));
                    datetick('x', 'yyyy:mm', 'keeplimits', 'keepticks')
                else
                    xticks(forecastDates([1, 3 : 3 : end]))
                    xlim(forecastDates([1 end]));
                end


                if doActualMean

                    wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity-%s-MEAN', ...
                        ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                        wrap, [], [], [], [], true)
                    ht = title(datestr(dates(thisT), 'yyyy:mm'));
                    wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity-%s-MEAN-WITHDATE', ...
                        ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                        wrap, [], [], [], [], true)
                    box off
                    legend([hs h1 h1hat], 'shadow rate density', 'actual rate (median)',   'actual rate (mean)', 'location', 'northwest')
                    wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity-%s-MEAN-WITHDATELEGEND', ...
                        ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                        wrap, [], [], [], [], true)
                    delete(ht)
                    wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity-%s-MEAN-WITHLEGEND', ...
                        ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                        wrap, [], [], [], [], true)
                else

                    wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity-%s', ...
                        ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                        wrap, [], [], [], [], true)
                    ht = title(datestr(dates(thisT), 'yyyy:mm'));
                    wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity-%s-WITHDATE', ...
                        ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                        wrap, [], [], [], [], true)
                    box off
                    legend([hs h1 ], 'shadow rate density', 'actual rate (median)', 'location', 'northwest')
                    wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity-%s-WITHDATELEGEND', ...
                        ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                        wrap, [], [], [], [], true)
                    delete(ht)
                    wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity-%s-WITHLEGEND', ...
                        ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                        wrap, [], [], [], [], true)
                end
            end

            %% plot ShadowOnly
            if doShadowOnly
                thisfig = figure;
                ax = gca;
                set(ax, 'fontsize', fontsize)

                hold on

                plotCIaltcolor(shadowmid, shadowtails, forecastDates, [], 'w--', 'linewidth', 3);

                hs = plot(forecastDates, shadowmid, '-',  'color', [0.5843    0.8157    0.9882], 'linewidth', 4); %#ok<NASGU>
                plot(forecastDates, shadowmid, 'w--',  'linewidth', 4);


                %     hrealized = plot(forecastDates, realized, 'd', 'color', [0 .5 0], 'linewidth', 3);

                % h1 = plot(forecastDates, med0, 'b-', 'linewidth', 4);
                % h1tail = plot(forecastDates, tail0, 'b-', 'linewidth', 2);
                % if doShadowMean
                %     h1hat = plot(forecastDates, mid0, 'bd', 'linewidth', 4);
                % end

                if doDateTicks
                    xticks(forecastDates([1, 6 : 6 : end]))
                    xlim(forecastDates([1 end]));
                    datetick('x', 'yyyy:mm', 'keeplimits', 'keepticks')
                else
                    xticks(forecastDates([1, 3 : 3 : end]))
                    xlim(forecastDates([1 end]));
                end


                wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensityShadowOnly-%s', ...
                    ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                    wrap, [], [], [], [], true)
                ht = title(datestr(dates(thisT), 'yyyy:mm')); %#ok<NASGU>
                wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensityShadowOnly-%s-WITHDATE', ...
                    ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                    wrap, [], [], [], [], true)
                box off
                % wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensityShadowOnly-%s-WITHDATELEGEND', ...
                %     ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                %     wrap, [], [], [], [], true)
                % delete(ht)
                % wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensityShadowOnly-%s-WITHLEGEND', ...
                %     ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                %     wrap, [], [], [], [], true)
            end

            %% plot w/linear
            thisfig = figure;
            ax = gca;
            set(ax, 'fontsize', fontsize)

            hold on

            plotCIaltcolor(shadowmid, shadowtails, forecastDates, [], 'w--', 'linewidth', 3);
            hs = plot(forecastDates, shadowmid, '-',  'color', [0.5843    0.8157    0.9882], 'linewidth', 4);
            plot(forecastDates, shadowmid, 'w--',  'linewidth', 4);

            % linear
            hlin = plot(forecastDates, linearmid, 'k-.',  'linewidth', 4);
            plot(forecastDates, lineartails, 'k-.',  'linewidth', 2);

            %     hrealized = plot(forecastDates, realized, 'd', 'color', [0 .5 0], 'linewidth', 3);

            h1 = plot(forecastDates, med0, 'b-', 'linewidth', 4);
            h1tail = plot(forecastDates, tail0, 'b-', 'linewidth', 2);
            % if doActualMean
            %     h1hat = plot(forecastDates, mid0, 'bd', 'linewidth', 4);
            % end

            if doDateTicks
                xticks(forecastDates([1, 6 : 6 : end]))
                xlim(forecastDates([1 end]));
                datetick('x', 'yyyy:mm', 'keeplimits', 'keepticks')
            else
                xticks(forecastDates([1, 3 : 3 : end]))
                xlim(forecastDates([1 end]));
            end


            wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity1-%s', ...
                ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                wrap, [], [], [], [], true)
            ht = title(datestr(dates(thisT), 'yyyy:mm'));
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity1-%s-WITHDATE', ...
                ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                wrap, [], [], [], [], true)
            box off
            % if doActualMean
            %     hl = legend([hs h1 h1hat hlin], 'shadow rate density', 'actual rate (median)', 'actual rate (mean)', 'linear VAR', 'location', 'northwest');
            % else
            hl = legend([hs h1  hlin], 'shadow rate density', 'actual rate (median)', 'linear VAR', 'location', 'northwest');
            % end
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity1-%s-WITHDATELEGEND', ...
                ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                wrap)
            delete(ht)
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity1-%s-WITHLEGEND', ...
                ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                wrap, [], [], [], [], true)

            if doLinearActualOnly
                %% plot w/linear
                thisfig = figure;
                ax = gca;
                set(ax, 'fontsize', fontsize)

                hold on

                % plotCIaltcolor(shadowmid, shadowtails, forecastDates, [], 'w--', 'linewidth', 3);
                % hs = plot(forecastDates, shadowmid, '-',  'color', [0.5843    0.8157    0.9882], 'linewidth', 4);
                % plot(forecastDates, shadowmid, 'w--',  'linewidth', 4);

                % linear
                hlin = plot(forecastDates, linearmid, 'k-.',  'linewidth', 4);
                plot(forecastDates, lineartails, 'k-.',  'linewidth', 2);

                %     hrealized = plot(forecastDates, realized, 'd', 'color', [0 .5 0], 'linewidth', 3);

                h1 = plot(forecastDates, med0, 'b-', 'linewidth', 4);
                h1tail = plot(forecastDates, tail0, 'b-', 'linewidth', 2);

                if doDateTicks
                    xticks(forecastDates([1, 6 : 6 : end]))
                    xlim(forecastDates([1 end]));
                    datetick('x', 'yyyy:mm', 'keeplimits', 'keepticks')
                else
                    xticks(forecastDates([1, 3 : 3 : end]))
                    xlim(forecastDates([1 end]));
                end


                wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity2-%s', ...
                    ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                    wrap, [], [], [], [], true)
                ht = title(datestr(dates(thisT), 'yyyy:mm'));
                wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity2-%s-WITHDATE', ...
                    ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                    wrap, [], [], [], [], true)
                box off
                hl = legend([h1  hlin], 'actual rate (median)', 'linear VAR', 'location', 'northwest');
                wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity2-%s-WITHDATELEGEND', ...
                    ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                    wrap)
                delete(ht)
                wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity2-%s-WITHLEGEND', ...
                    ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                    wrap, [], [], [], [], true)

            end

            if doLinearOnly
                % plot lin only
                YLIM = ylim;
                thisfig = figure;
                ax = gca;
                set(ax, 'fontsize', fontsize)

                hold on


                % linear
                hlin = plot(forecastDates, linearmid, 'k-.',  'linewidth', 4);
                plot(forecastDates, lineartails, 'k-.',  'linewidth', 2);

                ylim(YLIM)
                if doDateTicks
                    xticks(forecastDates([1, 6 : 6 : end]))
                    xlim(forecastDates([1 end]));
                    datetick('x', 'yyyy:mm', 'keeplimits', 'keepticks')
                else
                    xticks(forecastDates([1, 3 : 3 : end]))
                    xlim(forecastDates([1 end]));
                end
                wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity1lin-%s', ...
                    ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                    wrap, [], [], [], [], true)
                title(datestr(dates(thisT), 'yyyy:mm'))
                wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity1lin-%s-WITHDATE', ...
                    ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                    wrap, [], [], [], [], true)
            end

            if doLinearShadowOnly
                % plot lin only
                YLIM = ylim;
                thisfig = figure;
                ax = gca;
                set(ax, 'fontsize', fontsize)

                hold on

                % shadow
                plotCIaltcolor(shadowmid, shadowtails, forecastDates, [], 'w--', 'linewidth', 3);
                hs = plot(forecastDates, shadowmid, '-',  'color', [0.5843    0.8157    0.9882], 'linewidth', 4);
                plot(forecastDates, shadowmid, 'w--',  'linewidth', 4);

                % linear
                hlin = plot(forecastDates, linearmid, 'k-.',  'linewidth', 4);
                plot(forecastDates, lineartails, 'k-.',  'linewidth', 2);

                ylim(YLIM)
                if doDateTicks
                    xticks(forecastDates([1, 6 : 6 : end]))
                    xlim(forecastDates([1 end]));
                    datetick('x', 'yyyy:mm', 'keeplimits', 'keepticks')
                else
                    xticks(forecastDates([1, 3 : 3 : end]))
                    xlim(forecastDates([1 end]));
                end
                wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity1shadow-%s', ...
                    ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                    wrap, [], [], [], [], true)
                hl = legend([hs hlin], 'shadow rate density', 'linear VAR', 'location', 'northwest');
                wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity1shadow-%s-WITHLEGEND', ...
                    ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                    wrap, [], [], [], [], true)

                title(datestr(dates(thisT), 'yyyy:mm'))
                wrapthisfigure(thisfig, sprintf('%s-%s-%s-predictivedensity1shadow-%s-WITHDATELEGEND', ...
                    ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm')), ...
                    wrap, [], [], [], [], true)

            end


            %% tabulate
            tabname = sprintf('%s-%s-%s-predictivedensity1lin-%s', ncode{thisY}, datalabel, shadowshortlabel, datestr(dates(thisT), 'yyyy-mm'));
            %             call tabulateFigure(wrap, thisPlotLabel, n, ndxT, thisdateT, theseFcstdates, prettylabel1, prettylabel2, prettylabel3, fcstMid1, fcstMid2, fcstMid3, fcstTails1, fcstTails2, fcstTails3, ncode, ndxCI, setQuantiles)

            labels4table = {sprintf('%s median', 'SR'), sprintf('%s %4.2f\\%%', 'SR', mat0.setQuantiles(1,ndxCI(1))), sprintf('%s %4.2f\\%%', 'SR', mat0.setQuantiles(1,ndxCI(2))), ...
                sprintf('%s median', 'LIN'), sprintf('%s %4.2f\\%%', 'LIN', mat0.setQuantiles(1,ndxCI(1))), sprintf('%s %4.2f\\%%', 'LIN', mat0.setQuantiles(1,ndxCI(2)))};

            data4table = [shadowmid(:) shadowtails linearmid(:) lineartails];
            writedatatable(wrap, tabname, forecastDates, data4table, labels4table, 'yyyy:mm');
            writedatatable2tex(wrap, tabname, forecastDates, data4table, labels4table, 'yyyy:mm');

        end


        %% finish script
        if ~isdesktop
            close all
        end
        finishwrap

    end % doHybrid
end % datalabel

dockAllFigures
finishscript
