%% load em toolboxes
warning('off','MATLAB:handle_graphics:exceptions:SceneNode')

path(pathdef)

addpath matlabtoolbox/emtools/
addpath matlabtoolbox/emtexbox/
addpath matlabtoolbox/emgibbsbox/
addpath matlabtoolbox/emeconometrics/
addpath matlabtoolbox/emstatespace/

%% prep
clear; close all; clc;

fontsize = 16;

%% load stuff
%#ok<*UNRCH>


for MODELTYPE = {'simpleshadowrateVAR', 'hybridshadowrateVAR'} 
    
    modeltype = MODELTYPE{:};
    
    
    for DATALABEL =  {'fredMD20baa-2021-07'} %{'fredMD20baaExYield-2021-07''}
        
        datalabel = DATALABEL{:};
        
        
        titlename = sprintf('showQRTshadow-%s-%s', datalabel, modeltype);
        initwrap
        
        
        mat = matfile(sprintf('~/jam/lager/quanticoELBmatfiles2021cum/%s-%s-tightBVARshrinkage-p12.mat', datalabel, modeltype));
        
        shadowshortlabel = 'shadowrate';
        
        Nvin = size(mat.shadowrateQRTmid,3);
        
        firstQRTobs = mat.Tjumpoffs(1,1);
        ydates      = mat.ydates;
        p           = mat.p;
        
        shadowrateQRTmid = mat.shadowrateQRTmid;
        shadowrateQRTtails = mat.shadowrateQRTtails;
        shadowrateVintagesMid = mat.shadowrateVintagesMid;
        shadowrateVintagesTails = mat.shadowrateVintagesTails;
        
        ndxSHADOWRATE = mat.ndxSHADOWRATE;
        
        %% qrt vs final
        for n = 1 : length(ndxSHADOWRATE)
            
            thisfig = figure;
            set(gca, 'fontsize', fontsize)
            
            hold on
            hfinal = plotCI(shadowrateVintagesMid(:,n,end), squeeze(shadowrateVintagesTails(:,n,[1 4],end)), ydates, [], 'k-', 'linewidth', 3);
            
            hqrt      = plot(ydates, shadowrateQRTmid(:,n), 'r-', 'linewidth', 3);
            hqrttails = plot(ydates, squeeze(shadowrateQRTtails(:,n,[1 4])), 'r-.', 'linewidth', 2);
            
            hl = legend([hfinal hqrt], 'full sample', 'quasi-real time', 'location', 'northwest', 'box', 'off');
            
            xtickdates(ydates([firstQRTobs end]))
            YLIM = ylim;
            wrapthisfigure(thisfig, sprintf('%s%d-QRTp%d-%s-%s', shadowshortlabel, n, p, datalabel, modeltype), wrap)
            
            delete(hqrt); delete(hqrttails); delete(hl);
            ylim(YLIM);
            wrapthisfigure(thisfig, sprintf('%s%d-FINALp%d-%s-%s', shadowshortlabel, n, p, datalabel, modeltype), wrap)
            
            %             ylim([-11 3])
            %             yticks(-12:2:2)
            %             wrapthisfigure(thisfig, sprintf('%s%d-FINALp%d-%s-commonlim', shadowshortlabel, n, p, datalabel), wrap)
            %             hqrt      = plot(ydates, shadowrateQRTmid(:,n), 'r-', 'linewidth', 3);
            %             hqrttails = plot(ydates, squeeze(shadowrateQRTtails(:,n,[1 4])), 'r-.', 'linewidth', 2);
            %             ylim([-16 3])
            %             yticks(-16:2:2)
            %             wrapthisfigure(thisfig, sprintf('%s%d-QRTp%d-%s-commonlim', shadowshortlabel, n, p, datalabel), wrap)
            %             hl = legend([hfinal hqrt], 'full sample', 'quasi-real time', 'location', 'northwest', 'box', 'off');
            %             wrapthisfigure(thisfig, sprintf('%s%d-QRTp%d-%s-commonlim-WITHLEGEND', shadowshortlabel, n, p, datalabel), wrap)
            
        end
        
        %% tabulate some numbes
        n = 1;
        SRfinal = shadowrateVintagesMid(:,n,end);
        T = length(ydates);
        tELB = find(ydates == datenum(2009,1,1));
        for t = tELB : T
            fprintf('%s \t %6.2f \n', datestr(ydates(t), 'yyyy:mm'), SRfinal(t));
        end
        
        %% patch in wuxia
        wuxia = importdata('data/shadowrate_US.xls');
        wuxiaY = floor(wuxia(:,1) / 100);
        wuxiaM = wuxia(:,1) - wuxiaY * 100;
        wuxiaDates = datenum(wuxiaY, wuxiaM, 1);
        wuxiaRate = wuxia(:,2);
        
        %% patch in krippner
        krippner         = importdata('data/SSR_Estimates_20210726.xlsx');
        krippnerDatesEOM = xls2mdate(krippner.data.D0x2EMonthlyAverageSSRSeries(:,1));
        ndxSSR = 2;
        if ~strcmpi('US SSR', krippner.textdata.D0x2EMonthlyAverageSSRSeries(7,ndxSSR))
            error('check out Krippner XLS sheet')
        end
        krippnerRate = krippner.data.D0x2EMonthlyAverageSSRSeries(:,ndxSSR);
        
        [y, m] = datevec(krippnerDatesEOM);
        krippnerDates = datenum(y,m,1);
        
        %% plot FINAL vs WUXI and Krippner
        thisfig = figure;
        set(gca, 'fontsize', fontsize)
        
        hold on
        hfinal = plotCI(shadowrateVintagesMid(:,n,end), squeeze(shadowrateVintagesTails(:,n,[1 4],end)), ydates, [], 'k-', 'linewidth', 3);
        
       
        xtickdates(ydates([firstQRTobs end]))
        YLIM = ylim;
        hwx = plot(wuxiaDates, wuxiaRate, 'm-.', 'linewidth', 2);
        hkrip = plot(krippnerDates, krippnerRate, 'b-.', 'linewidth', 2);
        ylim([-6 3])
        legend([hfinal hwx hkrip], 'Shadow-rate VAR', 'Wu-Xia', 'Krippner', 'location', 'northwest', 'box', 'off')
        wrapthisfigure(thisfig, sprintf('%s%d-FINALp%d-%s-%s-wuxiakrippner', shadowshortlabel, n, p, datalabel, modeltype), wrap)
        
        
        
        %% wrap up
        dockAllFigures
        finishwrap
        % wrap = latexwrapper(wrap, 'close');
        
    end
end
