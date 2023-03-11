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

rng(01012023)

%% Initial operations
clear; close all; clc;


%% set parameters for VAR and MCMC
datalabel           = 'fredMD20VXO-2022-09';
jumpDate            = datenum(2022,08,01);
check_stationarity  = 0;                  % Truncate nonstationary draws? (1=yes)

% resultsdir           = 'mcmcMatfiles';
resultsdir          = 'mcmcVXOEBP';

irfSCALES           = [1 5 10];
irfDATES            = [datenum(2007,1,1) datenum(2009,1,1) datenum([2010 2012 2014],12,1)];

p                   = 12;                    % Number of lags on dependent variables
irfNdraws           = 1e3;                   % per MCMC node
irfHorizon          = 25;

% SED-PARAMETERS-HERE

doRATSprior         = true;

samStart            = [];                 % truncate start of sample if desired (leave empty if otherwise)

doELBsampling       = true;
ELBbound            = 0.25;

np = 12;

%% load data
% load CSV file
dum=importdata(sprintf('%s.csv', datalabel),',');


ydates=dum.data(3:end,1);
% Variable names
ncode=dum.textdata(1,2:end);
% Transformation codes (data are already transformed)
tcode  =dum.data(1,2:end);
cumcode=logical(dum.data(2,2:end));
cumcode(tcode == 5) = 1;
% Data
data=dum.data(3:end,2:end);

setShadowYields

ndxYIELDS = union(ndxSHADOWRATE, ndxOTHERYIELDS);
Nyields   = length(ndxYIELDS);

Nshadowrates = length(ndxSHADOWRATE);
Tdata = length(ydates);

Ylabels = fredMDprettylabel(ncode);

%% process settings
N     = size(data,2);
K     = N * p + 1; % number of regressors per equation
Nstates = K + p * Nyields;

% ndxSHADOWRATELAGS = cat(2, false, repmat(ismember(1:N, ndxSHADOWRATE), 1, p)); % prepend by false for CONST
ndxYIELDLAGS = cat(2, false, repmat(ismember(1:N, ndxYIELDS), 1, p)); % prepend by false for CONST

% truncate start of sample (if desired)
if ~isempty(samStart)
    ndx    = ydates >= samStart;
    data   = data(ndx,:);
    ydates = ydates(ndx);
    Tdata  = length(ydates);
end

% define oos jump offs

ELBdummy = data(:,ndxSHADOWRATE) <= ELBbound;

startELB       = find(any(ELBdummy,2), 1);
elbT0          = startELB - 1 - p;
% elbT0: first obs prior to missing obs, this is the jump off for the state space
% note: startELB is counted against the available obs in sample, which include
% p additional obs compared to the VAR



%% some parameters

fontsize = 12;

TID   = parid;

thisT = find(ydates == jumpDate);
T     = thisT - p;

setQuantiles    = [.5, 2.5, 5, normcdf(-1) * 100, 25 , 75,  (1 - normcdf(-1)) * 100, 95, 97.5, 99.5];
Nquantiles      = length(setQuantiles);
ndxCI68         = ismember(setQuantiles, [normcdf(-1) * 100, 100 - normcdf(-1) * 100]);
ndxCI90         = ismember(setQuantiles, [5 95]);
ndxCI           = ndxCI68 | ndxCI90;

rndStream   = getDefaultStream;

actualrateBlock = ~ismember(1:N, ndxYIELDS);

%% collect MCMC results
titlename=sprintf('%s-p%d-jumpoff%s', datalabel, p, datestr(jumpDate, 'yyyymmm'));

mcmc = matfile(fullfile(resultsdir, sprintf('mcmcBlockHybrid-%s', titlename)));

MCMCdraws = mcmc.MCMCdraws;

% reload and clear the following with each simulation:
% PAIdraws  = permute(mcmc.PAI_all, [3 2 1]);
% invAdraws = permute(mcmc.invA_all, [2 3 1]);
% PHIdraws  = permute(mcmc.PHI_all, [2 1]);
%
% shadowrate_all = mcmc.shadowrate_all;
% sqrtht_all     = mcmc.sqrtht_all;

%% shock draws
zdraws   = randn(rndStream, N, irfHorizon, irfNdraws);
zSVdraws = randn(rndStream, N, irfHorizon * irfNdraws);

%% shock size
switch datalabel
    case {'fredMD20VXO-2022-09', 'fredMD20VXOexYield-2022-09'}
        shocksize = 1.27;
    case {'fredMD20EBP-2022-09', 'fredMD20EBPexYield-2022-09'}
        shocksize = 0.11;
    otherwise
        shocksize = 1;
end

%% GIRF

for IRF1scale = irfSCALES

    display(IRF1scale)

    for irfDate = irfDATES

        display(datestr(irfDate, 'yyyymmm'))


        % prepare wrap
        titlename=sprintf('%s-p%d-IRF1scale%d-jumpoff%s-irfDate%s', datalabel, p, IRF1scale, ...
            datestr(jumpDate, 'yyyymmm'), datestr(irfDate, 'yyyymmm'));

        if ~isempty(samStart)
            titlename = strcat(titlename,'-', datestr(samStart, 'yyyymmm'));
        end

        wrap = [];
        initwrap

        %% load mcmc
        PAIdraws  = permute(mcmc.PAI_all, [3 2 1]);
        invAdraws = permute(mcmc.invA_all, [2 3 1]);
        PHIdraws  = permute(mcmc.PHI_all, [2 1]);
        [~, ndxvech]    = ivech(PHIdraws(:,1));

        ndxIRFT0               = find(ydates == irfDate);
        % to accomodate pre 2009 jumpoffs, need to patch in actual rate history
        if ndxIRFT0 - elbT0 - p > 0
            shadowrateJumpoffDraws = permute(mcmc.shadowrate_all(:,:,1:ndxIRFT0 - elbT0 - p), [3 2 1]);
        else
            shadowrateJumpoffDraws = [];
        end

        SVjumpoffDraws         = permute(mcmc.sqrtht_all(:,ndxIRFT0 - p,:), [3 1 2]);

        %% allocate memory
        [fcstYdraws, fcstYdraws1plus, fcstYdraws1minus] = deal(NaN(N, irfHorizon, irfNdraws, MCMCdraws));

        % prepare state space for forecast simulation
        ndxfcstY          = 1+(1:N);
        fcstB             = zeros(Nstates,N);
        fcstB(ndxfcstY,:) = eye(N);

        fcstA0                  = zeros(Nstates,Nstates);
        fcstA0(1,1)             = 1; % unit root for the constant
        fcstA0(1+N+1:K,2:K)     = [eye(N*(p-1)),zeros(N*(p-1),N)]; % fill in lower part of companion form
        fcstA0(K+Nyields+(1:Nyields*(p-1)),K+(1:Nyields*(p-1))) = eye(Nyields*(p-1)); % companion for actual rates
        ndxfcstActual          = K+(1:Nyields);
        ndxfcstShadow          = 1+ndxYIELDS;


        % construct forecast jumpoff (with placeholders for shadow rates)
        jumpoffDate       = ydates(thisT);
        ndx               = ydates <= jumpoffDate;
        jumpoffData       = data(ndx,:);

        prc70 = normcdf([-1 1]) * 100;

        %% loop over MCMC draws
        %         progressbar(0)
        parfor mm = 1 : MCMCdraws

            TID       = parid;

            % parfor preps (better to do inside parfor loop)
            thisData               = jumpoffData;
            fcstA                  = fcstA0;

            % construct jump off vector
            if ~isempty(shadowrateJumpoffDraws) % identical to p+elbT0+1 < 0
                thisData(p+elbT0+1:ndxIRFT0,ndxSHADOWRATE) = shadowrateJumpoffDraws(:,:,mm);
            end
            Xjumpoff    = zeros(Nstates,1);
            Xjumpoff(1) = 1;
            for l=1:p
                Xjumpoff(1+(l-1)*N+(1:N)) = thisData(ndxIRFT0-(l-1),1:N);
            end
            % add lagged actual rates to jumpoff
            for l=1:p
                Xjumpoff(K+(l-1)*Nyields+(1:Nyields)) = thisData(ndxIRFT0-(l-1),ndxYIELDS);
            end
            fcstX0                        = Xjumpoff;

            % map into state space
            PAIactual                     = PAIdraws(:,ndxYIELDLAGS,mm);
            PAIactual(~actualrateBlock,:) = 0;

            PAIshadow = PAIdraws(:,:,mm); %#ok<PFBNS> 
            PAIshadow(actualrateBlock,ndxYIELDLAGS) = 0;

            fcstA(ndxfcstY, 1:K)         = PAIshadow;
            fcstA(ndxfcstY, K+1:Nstates) = PAIactual;

            PHI     = ivech(PHIdraws(:,mm), ndxvech);

            % generate SV paths
            sqrtPHI     = chol(PHI, 'lower');
            logSV       = sqrtPHI * zSVdraws;
            logSV       = reshape(logSV, N, irfHorizon, irfNdraws);
            logSV       = cumsum(logSV,2);
            fcstSVdraws = exp(logSV * 0.5) .* SVjumpoffDraws(:,mm);
            nushocks    = fcstSVdraws .* zdraws;

            % baseline
            fcstYdraws(:,:,:,mm)       = simVARshadowrateBlockHybrid(N, fcstX0, ndxfcstY, ndxfcstActual, ndxfcstShadow, cumcode, np, fcstA, fcstB, invAdraws(:,:,mm), ...
                nushocks, ELBbound, ndxYIELDS, irfHorizon, irfNdraws);
            % positive shock
            nushocks(1,1,:)            = IRF1scale * shocksize;
            fcstYdraws1plus(:,:,:,mm)  = simVARshadowrateBlockHybrid(N, fcstX0, ndxfcstY, ndxfcstActual, ndxfcstShadow, cumcode, np, fcstA, fcstB, invAdraws(:,:,mm), ...
                nushocks, ELBbound, ndxYIELDS, irfHorizon, irfNdraws);
            % negative shock
            nushocks(1,1,:)            = -1 * IRF1scale * shocksize;
            fcstYdraws1minus(:,:,:,mm) = simVARshadowrateBlockHybrid(N, fcstX0, ndxfcstY, ndxfcstActual, ndxfcstShadow, cumcode, np, fcstA, fcstB, invAdraws(:,:,mm), ...
                nushocks, ELBbound, ndxYIELDS, irfHorizon, irfNdraws);

            %             progressbar(mm / MCMCdraws)

        end % parfor

        clear logSV fcstSVdraws nushocks % obsolete in case of parfor
        clear PAIdraws invAdraws PHIdraws shadowrateJumpoffDraws SVjumpoffDraws
        clear fcstA* fcstB*

        %% forecast moments
        % first, integrate only per MCMC node
        fcstYhat       = mean(fcstYdraws, 3);
        fcstYhat1plus  = mean(fcstYdraws1plus, 3);
        fcstYhat1minus = mean(fcstYdraws1minus, 3);

        clear fcstYdraws fcstYdraws1plus fcstYdraws1minus

        %% IRF
        IRFdraws1plus   = fcstYhat1plus  - fcstYhat;
        IRFdraws1minus  = fcstYhat1minus - fcstYhat;
        
        % integrate over MCMC nodes
        fcstYhat       = mean(fcstYhat, 4);
        fcstYhat1plus  = mean(fcstYhat1plus, 4);
        fcstYhat1minus = mean(fcstYhat1minus, 4);       

        IRF1plus        = mean(IRFdraws1plus, 4);
        IRF1plusTails   = prctile(IRFdraws1plus, prc70, 4);

        IRF1minus       = mean(IRFdraws1minus, 4);
        IRF1minusTails  = prctile(IRFdraws1minus, prc70, 4);

        deltaIRFdraws   = IRFdraws1plus + IRFdraws1minus;
        deltaIRF        = median(deltaIRFdraws,4);
        deltaIRFtails   = prctile(deltaIRFdraws, prc70, 4);

        clear deltaIRFdraws IRFdraws1plus IRFdraws1minus

        %% PLOT RESULTS
        colorPlus     = Colors4Plots(1);
        colorMinus    = Colors4Plots(2);
        colorBase     = Colors4Plots(8);

        %% plot ELB IRF
        for n = 1 : N

            thisfig = figure;
            subplot(2,1,1)
            hold on
            set(gca, 'FontSize', fontsize)
            hplus  = plot(0:irfHorizon-1, IRF1plus(n,:), '-', 'color', colorPlus, 'linewidth', 3);
            plot(0:irfHorizon-1, squeeze(IRF1plusTails(n,:,:,:)), '-', 'color', colorPlus, 'linewidth', 1);
            hminus = plot(0:irfHorizon-1, -1 * IRF1minus(n,:), '-.', 'color', colorMinus, 'linewidth', 3);
            plot(0:irfHorizon-1, -1 * squeeze(IRF1minusTails(n,:,:,:)), '-.', 'color', colorMinus, 'linewidth', 1);
            xlim([0 irfHorizon-1])
            yline(0, 'k:')
            legend([hplus, hminus], 'response', 'inverted response to negative shock', 'location', 'southoutside')
            title('positive shock', 'FontWeight', 'normal')

            subplot(2,1,2)
            hold on
            set(gca, 'FontSize', fontsize)
            hminus = plot(0:irfHorizon-1, IRF1minus(n,:), '-', 'color', colorMinus, 'linewidth', 3);
            plot(0:irfHorizon-1, squeeze(IRF1minusTails(n,:,:,:)), '-', 'color', colorMinus, 'linewidth', 1);
            hplus  = plot(0:irfHorizon-1, -1 * IRF1plus(n,:), '-.', 'color', colorPlus, 'linewidth', 3);
            plot(0:irfHorizon-1, -1 * squeeze(IRF1plusTails(n,:,:,:)), '-.', 'color', colorPlus, 'linewidth', 1);
            xlim([0 irfHorizon-1])
            yline(0, 'k:')
            legend([hminus, hplus], 'response', 'inverted response to positive shock', 'location', 'southoutside')
            title('negative shock', 'FontWeight', 'normal')

            sgtitle(sprintf('%s per %s', Ylabels{n}, datestr(irfDate, 'yyyymmm')), 'FontSize', 18', 'FontWeight', 'bold')

            wrapthisfigure(thisfig, sprintf('IRF1plusminus-%s-IRF1scale%d-%s-jumpoff%s-irfDate%s', datalabel, IRF1scale, ncode{n}, ...
                datestr(jumpDate, 'yyyymmm'), datestr(irfDate, 'yyyymmm')), wrap)
        end
        close all


        %% plot delta ELB IRF
        for n = 1 : N

            thisfig = figure;
            hold on
            set(gca, 'FontSize', fontsize)

            plot(0:irfHorizon-1, deltaIRF(n,:), '-', 'color', colorPlus, 'linewidth', 2);
            plot(0:irfHorizon-1, squeeze(deltaIRFtails(n,:,:,:)), '-', 'color', colorPlus, 'linewidth', 1);

            xlim([0 irfHorizon-1])
            yline(0, 'k:')

            title(sprintf('%s', Ylabels{n}))

            wrapthisfigure(thisfig, sprintf('deltaIRFplusminus-%s-IRF1scale%d-%s-jumpoff%s-irfDate%s', datalabel, IRF1scale, ncode{n}, ...
                datestr(jumpDate, 'yyyymmm'), datestr(irfDate, 'yyyymmm')), wrap)

        end
        close all


        %% Response paths
        for n = 1 : N

            thisfig = figure;
            hold on
            set(gca, 'FontSize', fontsize)

            hbase  = plot(0:irfHorizon-1, fcstYhat(n,:), '-', 'color', colorBase, 'linewidth', 2);
            hplus  = plot(0:irfHorizon-1, fcstYhat1plus(n,:), '-', 'color', colorPlus, 'linewidth', 2);
            hminus = plot(0:irfHorizon-1, fcstYhat1minus(n,:), '-', 'color', colorMinus, 'linewidth', 2);

            xlim([0 irfHorizon-1])

            if any(ndxYIELDS == n) && ~(all(ylim > ELBbound) || all(ylim < ELBbound))
                hELB = yline(ELBbound, ':', 'color', Colors4Plots(5), 'LineWidth',2);
                legend([hbase hplus hminus hELB], 'baseline', 'positive shock', 'negative shock', 'ELB', ...
                    'location', 'best')
            else
                legend([hbase hplus hminus], 'baseline', 'positive shock', 'negative shock', 'location', 'best')
            end
            title(sprintf('%s per %s', Ylabels{n}, datestr(irfDate, 'yyyymmm')), 'FontWeight', 'normal')
            wrapthisfigure(thisfig, sprintf('pathResponses1plusminus-%s-IRF1scale%d-%s-jumpoff%s-irfDate%s', datalabel, IRF1scale, ncode{n}, ...
                datestr(jumpDate, 'yyyymmm'), datestr(irfDate, 'yyyymmm')), wrap)
        end
        close all

        %% wrap up
        allw = whos;
        ndx = contains({allw.class}, 'Figure');
        if any(ndx)
            clear(allw(ndx).name)
        end
        save(sprintf('irfBlockHybrid-%s.mat', titlename), '-v7.3')

        close all
        finishwrap

    end % irfDate
end % irfscale

finishscript


%% define forecast simulation as function
function ydraws = simVARshadowrateBlockHybrid(N, fcstX0, ndxfcstY, ndxfcstActual, ndxfcstShadow,  cumcode, np, fcstA, fcstB, invA, nushocks, ELBbound, ndxYIELDS, irfHorizon, irfNdraws)

ydraws      = NaN(N,irfHorizon,irfNdraws);

for nn = 1 : irfNdraws

    xlag        = fcstX0;
    theseShocks = invA * nushocks(:,:,nn);

    for hh = 1 : irfHorizon
        xdraw    = fcstA * xlag + fcstB * theseShocks(:,hh);
        ydraws(:,hh,nn)        = xdraw(ndxfcstY);
        % prepare next iteration
        xlag                   = xdraw;
        xlag(ndxfcstActual)    = max(xlag(ndxfcstShadow), ELBbound);
    end
end

yieldDraws                   = ydraws(ndxYIELDS,:,:);
ndx                          = yieldDraws < ELBbound;
yieldDraws(ndx)              = ELBbound;
ydraws(ndxYIELDS,:,:)        = yieldDraws;

ydraws(cumcode, :,:)         = cumsum(ydraws(cumcode,:,:), 2) / np;

end % function simVARshadowrate