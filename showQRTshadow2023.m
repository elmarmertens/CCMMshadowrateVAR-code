%% load em toolboxes
warning('off','MATLAB:handle_graphics:exceptions:SceneNode')

path(pathdef)

addpath matlabtoolbox/emtools/
addpath matlabtoolbox/emtexbox/
addpath matlabtoolbox/emgibbsbox/
addpath matlabtoolbox/emeconometrics/
addpath matlabtoolbox/emstatespace/

%#ok<*DATNM>
%#ok<*DATST>

%% prep
clear; close all; clc;

fontsize = 14;


ELBbound = .125;
ELBcolor = Colors4Plots(8);

%% get LSAP dates
[lsapDates, lsapLabels] = getLSAPdates();
lsapDates = datenum(lsapDates);
ndxDropTaper = ~strcmpi(lsapLabels, 'taper begins');
lsapLabels   = lsapLabels(ndxDropTaper);
lsapDates    = lsapDates(ndxDropTaper);


tailNDX = [1 4]; % [1 4] for 90% or [2 3] for IQR

%% load stuff
%#ok<*UNRCH>


for MODELTYPE = {'ELBblockhybrid'}

    modeltype = MODELTYPE{:};


    for DATALABEL = {'fredblockMD20-2022-09', 'fredblockMD20exYield-2022-09'}

        datalabel = DATALABEL{:};

        switch ELBbound
            case .25
                resultsdir = '../matfilesShadowrateVAR/lagerFREDblock';
                ELBtag    = '';
                ELBlegend = '25 bp';
                switch datalabel
                    case 'fredblockMD20-2022-09'
                        YLIM      = [-8 4];
                    case {'fredMD20exYield-2022-09', 'fredblockMD20exYield-2022-09'}
                        if strcmpi(modeltype, 'ELBhybrid')
                            YLIM      = [-8 4];
                        else
                            YLIM      = [-30 15];
                        end
                    otherwise
                        YLIM      = [-4 4];
                end
            case .125
                resultsdir = '../matfilesShadowrateVAR/lagerFREDblock';
                ELBtag    = '-ELB125';
                ELBlegend = '12.5 bp';
                switch datalabel
                    case 'fredMD20-2022-09'
                        YLIM      = [-8 4];
                    case 'fredMD20exYield-2022-09'
                        YLIM      = [-12 4];
                    otherwise
                        YLIM      = [-12 4];
                end
            otherwise
                error('ELBbound value of %5.2f not recognized', ELBbound)
        end


        titlename = sprintf('showQRTshadow-%s-%s%s', datalabel, modeltype, ELBtag);
        initwrap


        mat = matfile(fullfile(resultsdir, sprintf('%s-%s%s-RATSbvarshrinkage-p12.mat', datalabel, modeltype, ELBtag)));

        shadowshortlabel = 'shadowrate';

        Nvin = size(mat.shadowrateQRTmid,3);

        firstQRTobs = mat.Tjumpoffs(1,1);
        ydates      = mat.ydates;
        T           = length(ydates);
        p           = mat.p;

        shadowrateQRTmid        = mat.shadowrateQRTmid;
        shadowrateQRTtails      = mat.shadowrateQRTtails;
        shadowrateVintagesMid   = mat.shadowrateVintagesMid;
        shadowrateVintagesTails = mat.shadowrateVintagesTails;

        ndxSHADOWRATE = mat.ndxSHADOWRATE;

        %% patch in actual rate data
        FREDmd = importdata(sprintf('%s.csv', datalabel),',');
        checkdiff(ydates, FREDmd.data(3:end,1));

        data           = FREDmd.data(3:end,2:end);
        ActualRateData = data(:,mat.ndxSHADOWRATE);

        % patch QRT
        ndxActual      = isnan(shadowrateQRTmid);
        shadowrateQRTmid(ndxActual) = ActualRateData(ndxActual);
        for nn = 1 : size(shadowrateQRTtails,3)
            this                       = shadowrateQRTtails(:,:,nn);
            this(ndxActual)            = ActualRateData(ndxActual);
            shadowrateQRTtails(:,:,nn) = this;
        end

        % patch FINAL
        ndxActual      = isnan(shadowrateVintagesMid(:,:,end));
        this = shadowrateVintagesMid(:,:,end);
        this(ndxActual) = ActualRateData(ndxActual);
        shadowrateVintagesMid(:,:,end) = this;
        for nn = 1 : size(shadowrateVintagesTails,3)
            this                       = shadowrateVintagesTails(:,:,nn,end);
            this(ndxActual)            = ActualRateData(ndxActual);
            shadowrateVintagesTails(:,:,nn,end) = this;
        end

        %% qrt vs final
        for n = 1 : length(ndxSHADOWRATE)

            thisfig = figure;
            set(gca, 'fontsize', fontsize)

            hold on
            hfinal = plotCI(shadowrateVintagesMid(:,n,end), squeeze(shadowrateVintagesTails(:,n,tailNDX,end)), ydates, [], 'k-', 'linewidth', 3);

            hqrt      = plot(ydates, shadowrateQRTmid(:,n), 'r-', 'linewidth', 3);
            hqrttails = plot(ydates, squeeze(shadowrateQRTtails(:,n,tailNDX)), 'r-.', 'linewidth', 2);


            % xtickdates(ydates([firstQRTobs end]))
            xtickdates([datenum(2008,8,1) ydates(end)])
            ylim(YLIM)
            hELB = yline(ELBbound, ':', 'color', ELBcolor);
            hl = legend([hfinal hqrt hELB], 'full sample', 'quasi-real time', ELBlegend, ...
                'location', 'northwest', 'box', 'off', 'AutoUpdate', 'off');
            wrapthisfigure(thisfig, sprintf('%s%d-QRTp%d-%s-%s%s', shadowshortlabel, n, p, datalabel, modeltype, ELBtag), wrap, [], [], [], [], true)
            hLSAP = xline(lsapDates, '--', lsapLabels, 'fontsize', fontsize, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center'); 
            set(hl, 'location', 'southeast', 'box','on');
            wrapthisfigure(thisfig, sprintf('%s%d-QRTp%d-%s-%s%s-LSAP', shadowshortlabel, n, p, datalabel, modeltype, ELBtag), wrap)

            delete(hLSAP); delete(hqrt); delete(hqrttails); delete(hl);
            ylim(YLIM);
            wrapthisfigure(thisfig, sprintf('%s%d-FINALp%d-%s-%s%s', shadowshortlabel, n, p, datalabel, modeltype, ELBtag), wrap, [], [], [], [], true)
            hLSAP = xline(lsapDates, '--', lsapLabels, 'fontsize', fontsize, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center'); %#ok<NASGU>
            wrapthisfigure(thisfig, sprintf('%s%d-FINALp%d-%s-%s%s-LSAP', shadowshortlabel, n, p, datalabel, modeltype, ELBtag), wrap, [], [], [], [], true)

        end

        %% patch in wuxia
        wuxia = importdata('WUXIASHADOWRATE.csv');
        wuxiaDates = wuxia.data(:,1);
        wuxiaRate  = wuxia.data(:,2);

        %% patch in krippner
        krippner = importdata('KRIPPNERSHADOWRATE.csv');
        krippnerDates = krippner.data(:,1);
        krippnerRate  = krippner.data(:,2);


        %% plot QRT vs WUXIA and Krippner
        n = 1;
        thisfig = figure;
        set(gca, 'fontsize', fontsize)

        hold on
        hfinal = plotCIredshades(shadowrateQRTmid(:,n,end), squeeze(shadowrateQRTtails(:,n,tailNDX,end)), ydates, [], 'r-', 'linewidth', 3);

        % xtickdates(ydates([firstQRTobs end]))
        xtickdates([datenum(2008,8,1) ydates(end)])
        hwx = plot(wuxiaDates, wuxiaRate, '-.', 'color', Colors4Plots('lightblue'), 'linewidth', 3);
        hkrip = plot(krippnerDates, krippnerRate, ':', 'color', Colors4Plots('darkblue'), 'linewidth', 3);
        ylim(YLIM)
        hELB = yline(ELBbound, ':', 'color', ELBcolor);
        hl = legend([hfinal hwx hkrip hELB], 'Shadow-rate VAR', 'Wu-Xia', 'Krippner', ELBlegend, ...
            'location', 'southeast', 'box', 'off', 'AutoUpdate', 'off');
        wrapthisfigure(thisfig, sprintf('%s%d-QRTp%d-%s-%s%s-wuxiakrippner', shadowshortlabel, n, p, datalabel, modeltype, ELBtag), wrap, [], [], [], [], true)

        hLSAP = xline(lsapDates, '--', lsapLabels, 'fontsize', fontsize, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center'); %#ok<NASGU>
        set(hl, 'location', 'southeast', 'box','on');
        wrapthisfigure(thisfig, sprintf('%s%d-QRTp%d-%s-%s%s-wuxiakrippner-LSAP', shadowshortlabel, n, p, datalabel, modeltype, ELBtag), wrap)

        %% plot FINAL vs WUXI and Krippner
        n = 1;
        thisfig = figure;
        set(gca, 'fontsize', fontsize)

        hold on
        hfinal = plotCI(shadowrateVintagesMid(:,n,end), squeeze(shadowrateVintagesTails(:,n,tailNDX,end)), ydates, [], 'k-', 'linewidth', 3);

        % xtickdates(ydates([firstQRTobs end]))
        xtickdates([datenum(2008,8,1) ydates(end)])
        hwx = plot(wuxiaDates, wuxiaRate, '-.', 'color', Colors4Plots('lightblue'), 'linewidth', 3);
        hkrip = plot(krippnerDates, krippnerRate, ':', 'color', Colors4Plots('darkblue'), 'linewidth', 3);
        ylim(YLIM)
        hELB = yline(ELBbound, ':', 'color', ELBcolor);
        hl = legend([hfinal hwx hkrip hELB], 'Shadow-rate VAR', 'Wu-Xia', 'Krippner', ELBlegend, ...
            'location', 'southeast', 'box', 'off', 'AutoUpdate', 'off');
        wrapthisfigure(thisfig, sprintf('%s%d-FINALp%d-%s-%s%s-wuxiakrippner', shadowshortlabel, n, p, datalabel, modeltype, ELBtag), wrap, [], [], [], [], true)

        hLSAP = xline(lsapDates, '--', lsapLabels, 'fontsize', fontsize, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center'); 
        set(hl, 'location', 'southeast', 'box','on');
        wrapthisfigure(thisfig, sprintf('%s%d-FINALp%d-%s-%s%s-wuxiakrippner-LSAP', shadowshortlabel, n, p, datalabel, modeltype, ELBtag), wrap)


        %% tabulate values
        tabname      = sprintf('%s%d-p%d-%s-%s', shadowshortlabel, n, p, datalabel, modeltype);
        data4table   = [shadowrateQRTmid(:,n) shadowrateVintagesMid(:,n,end)];
        labels4table = {'Shadow rate (QRT)', 'Shadow rate (final)'};
        sam          = firstQRTobs : T;
        writedatatable(wrap, tabname, ydates(sam), data4table(sam,:), labels4table, 'yyyy-mm');
        %         writedatatable2tex(wrap, tabname, ydates(sam), data4table(sam,:), labels4table, 'yyyy-mm');


        %% wrap up
        dockAllFigures
        finishwrap

        % wrap = latexwrapper(wrap, 'close');

    end
end
