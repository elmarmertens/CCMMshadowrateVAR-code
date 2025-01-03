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

%% Initial operations
clear; close all; clc;


resultsdir = pwd;

DATALABEL  = {'fredsxMD20-2022-09', 'fredsxMD20exYield-2022-09'};
MODELTYPES = {'standardVARAR1SV', 'ELBblocknonstructuralAR1SV'};


resultsdir0 = resultsdir;
modeltype0  = 'standardVARAR1SV200811';

doMedian = false;
doTitle  = true;

doLags        = false;
doEquations   = false;

doWrap = true;
%#ok<*NASGU>

maxIntercept  = [];
maxLag1       = [];
maxLagOTHER   = [];
maxALL        = [];

maxEq         = cell(20,1);


%% loop over DATALABEL
for thisd = 1 : length(DATALABEL)
    
    datalabel = DATALABEL{thisd};
    %% load model0
    mat0  = matfile(fullfile(resultsdir0, sprintf('%s-%s-RATSbvarshrinkage-p12.mat', datalabel, modeltype0)));
    
    %% loop over models
    for thism =  1 : length(MODELTYPES)
        
        modeltype = MODELTYPES{thism};
        
        mat   = matfile(fullfile(resultsdir, sprintf('%s-%s-RATSbvarshrinkage-p12.mat', datalabel, modeltype)));
        
        wrap = [];
        titlename = sprintf('PAIsinceGFC0-%s-%s', datalabel, modeltype);
        if doMedian
            titlename = strcat(titlename, '-median');
        end
        
        if doWrap
            initwrap
        end
        
        modeltype = strcat(modeltype, '-', datalabel);
        %#ok<*UNRCH>
        
        %% pull some objects
        Tjumpoffs = mat.Tjumpoffs;
        ydates    = mat.ydates;
        N         = mat.N;
        p         = mat.p;
        K         = N*p + 1;
        ncode     = mat.ncode;
        Ylabels   = fredMDprettylabel(ncode);
        
        % if ~isequal(Tjumpoffs, mat0.Tjumpoffs)
        %     error houston
        % end
        % if ~isequal(ydates, mat0.ydates)
        %     error houston
        % end
        % if ~isequal(ncode, mat0.ncode)
        %     error houston
        % end
        
        %% settings
        jumpoff0 = datenum(2009,1,1);
        T0       = find(ydates == jumpoff0);
        ndxT0    = find(Tjumpoffs == T0); % find to make matfile calls work
        
        %% collect PAI
        
        if ~doMedian
            PAI0   = mat0.PAImean(:,:,1);
            PAI0se = mat0.PAIstdev(:,:,1);
            PAIdevs = (mat.PAImean(:,:,ndxT0:end) - PAI0) ./ PAI0se;
        else
            PAI0   = mat0.PAImedian(:,:,1);
            PAI0se = mat0.PAIstdev(:,:,1);
            PAIdevs = (mat.PAImedian(:,:,ndxT0:end) - PAI0) ./ PAI0se;
        end
        
        maxRED = 2.5;
        
        %% plot PAI
        
        %     maxZ   = ceil(max(abs(PAIdevs),[], 'all') * 10) / 10;
        PAIdevs2D = reshape(PAIdevs(1:K,:,:), N * K, size(PAIdevs,3));
        thisfig = figure;
        
        hb = surf(ydates(T0:end),1:N*K,abs(PAIdevs2D));
        
        if isempty(maxALL)
            maxALL = max(zlim);
        else
            zlim([0 maxALL])
        end
        ylim([1 N * K])
        clim([0 maxRED])
        
        datetick('x')
        xlabel('end of estimation window')
        ylabel('parameters')
        
        % color by height
        shading interp
        colorbar
        colormap turbo
        set(gca, 'fontsize', 12)
        wrapthisfigure(thisfig, sprintf('PAIall-%s', modeltype), wrap);
        
        %% plot PAI per equation
        if doEquations
            for n = 1 : N
                
                % maxZ   = ceil(max(abs(PAIdevs),[], 'all') * 10) / 10;
                
                thisfig = figure;
                
                hb = surf(ydates(T0:end),1:K,squeeze(abs(PAIdevs(1:K,n,:))));
                
                maxthis = maxEq{n};
                if isempty(maxthis)
                    maxEq{n} = max(zlim);
                else
                    zlim([0 maxthis])
                end
                ylim([1 K])
                clim([0 maxRED])
                datetick('x')
                xlabel('end of estimation window')
                ylabel('parameters')
                
                
                
                % color by height
                shading interp
                colorbar
                colormap turbo
                
                wrapthisfigure(thisfig, sprintf('PAI-%s-%s', ncode{n}, modeltype), wrap, [], [], [], [], true);
                if doTitle
                    title(sprintf('%s equation', Ylabels{n}))
                    wrapthisfigure(thisfig, sprintf('PAI-%s-%s-WITHTITLE', ncode{n}, modeltype), wrap);
                end
            end
        end
        %% plot intercepts per equation
        thisfig = figure;
        
        thesedevs = squeeze(abs(PAIdevs(1,:,:)));
        
        % find variables with largest devs
        maxdev = max(thesedevs, [], 2);
        [maxdev,ndx] = sort(maxdev, 'desc');
        thesendx = sort(ndx(maxdev > 1.5));
        
        hb = surf(ydates(T0:end),1:N,thesedevs);
        
        if isempty(maxIntercept)
            maxIntercept = max(zlim);
        else
            zlim([0 maxIntercept])
        end
        ylim([1 N])
        clim([0 maxRED])
        
        datetick('x')
        yticks(thesendx)
        yticklabels(Ylabels(thesendx));
        
        xlabel('end of estimation window')
        
        
        % color by height
        shading interp
        colorbar
        colormap turbo
        wrapthisfigure(thisfig, sprintf('PAI-intercept-%s', modeltype), wrap, [], [], [], [], true);
        if doTitle
            title(sprintf('intercepts'))
            wrapthisfigure(thisfig, sprintf('PAI-intercept-%s-WITHTITLE', modeltype), wrap);
        end
        
        %% plot lag1 per equation
        thisfig = figure;
        
        thesedevs = reshape(abs(PAIdevs(1+(1:N),:,:)), N * N, []);
        
        % find variables with largest devs
        maxdev = max(thesedevs, [], 2);
        [maxdev,ndx] = sort(maxdev, 'desc');
        thesendx = sort(ndx(maxdev > 1.5));
        
        hb = surf(ydates(T0:end),1:N*N,thesedevs);
        
        if isempty(maxLag1)
            maxLag1 = max(zlim);
        else
            zlim([0 maxLag1])
        end
        
        clim([0 maxRED])
        
        datetick('x')
        yticks(1:N:N*N)
        ylim([1 N*N])
        yticklabels(Ylabels);
        
        
        xlabel('end of estimation window')
        
        
        % color by height
        shading interp
        colorbar
        colormap turbo
        wrapthisfigure(thisfig, sprintf('PAI-lag1-%s', modeltype), wrap, [], [], [], [], true);
        if doTitle
            title(sprintf('lag 1'))
            wrapthisfigure(thisfig, sprintf('PAI-lag1-%s-WITHTITLE', modeltype), wrap);
        end
        
        %% plot other lags per equation
        thisfig = figure;
        
        thesedevs = reshape(abs(PAIdevs(1+N+1:end,:,:)), N * N * (p - 1), []);
        
        % find variables with largest devs
        %     maxdev = max(thesedevs, [], 2);
        %     [maxdev,ndx] = sort(maxdev, 'desc');
        %     thesendx = sort(ndx(maxdev > 1.5));
        
        hb = surf(ydates(T0:end),1:N*N * (p - 1),thesedevs);
        
        if isempty(maxLagOTHER)
            maxLagOTHER = max(zlim);
        else
            zlim([0 maxLagOTHER])
        end
        clim([0 maxRED])
        
        datetick('x')
        yticks(1:N:N*N)
        ylim([1 N*N])
        yticklabels(Ylabels);
        
        xlabel('end of estimation window')
        
        
        % color by height
        shading interp
        colorbar
        colormap turbo
        wrapthisfigure(thisfig, sprintf('PAI-lagOTHER-%s', modeltype), wrap, [], [], [], [], true);
        if doTitle
            title(sprintf('other lags'))
            wrapthisfigure(thisfig, sprintf('PAI-lagOTHER-%s-WITHTITLE', modeltype), wrap);
        end
        
        %% plot PAI per block of each equation
        if doLags
            for n = 1 : N
                
                maxZ   = ceil(max(abs(PAIdevs(:,n,:)),[], 'all') * 10) / 10;
                
                for lag = 1 : p
                    
                    ndxlag = 1 + N * (lag-1) + (1 : N);
                    thisfig = figure;
                    
                    hb = surf(ydates(T0:end),1:N,squeeze(abs(PAIdevs(ndxlag,n,:))));
                    
                    zlim([0 maxZ])
                    clim([0 maxRED])
                    datetick('x')
                    xlabel('end of estimation window')
                    %             ylabel('parameters')
                    yticks(1:N)
                    ylim([1 N])
                    yticklabels(Ylabels)
                    
                    
                    
                    % color by height
                    shading interp
                    colorbar
                    colormap turbo
                    
                    wrapthisfigure(thisfig, sprintf('PAI-%s-lag%d-%s', ncode{n}, lag, modeltype), wrap, [], [], [], [], true);
                    if doTitle
                        title(sprintf('%s equation, lag %d', Ylabels{n}, lag))
                        wrapthisfigure(thisfig, sprintf('PAI-%s-lag%d-%s-WITHTITLE', ncode{n}, lag, modeltype), wrap);
                    end
                end
                if doWrap
                    close all
                end
            end
        end
        
        %% finish loop
        if doWrap
            close all % dockAllFigures
        end
        finishwrap
    end % MODELTYPES
end % DATALABEL

%% finish script
finishscript