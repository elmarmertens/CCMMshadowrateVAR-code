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
addpath matlabtoolbox/empbsbox/

rng(01012023)

%% Initial operations
clear; close all; clc;

%% set parameters for VAR and MCMC
datalabel           = 'fredblockMD20EBP-2022-09';
jumpDate            = datenum(2022,08,01);
irfDATES            = [datenum(2006,12,1) datenum(2012,12,1)];

irfSCALES           = 1;

irfNdraws           = 1e3;
irfHorizon          = 120;

resultsdir          = pwd;

doPlots = false;

ELBbound            = 0.25;
p                   = 12; 
np                  = 12;

samStart            = [];                 % truncate start of sample if desired (leave empty if otherwise)

% SED-PARAMETERS-HERE

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


%% shock size
switch datalabel
    case {'fredMD20VXO-2022-09', 'fredMD20VXOexYield-2022-09'}
        shocksize = 1.27;
    case {'fredblockMD20EBP-2022-09', ...
            'fredMD20EBP-2022-09', 'fredMD20EBPexYield-2022-09', 'fredMD3EBP-2022-09'}
        shocksize = 0.11;
    otherwise
        shocksize = 1;
end

%% GIRF

for IRF1scale = irfSCALES

    display(IRF1scale)

    for irfDate = irfDATES

        close all
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
        [fcstYHATdraws, fcstYHATdraws1plus, fcstYHATdraws1minus] = deal(NaN(N, irfHorizon, MCMCdraws));

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

        %% GIRF simulation across MCMC nodes
        parfor mm = 1 : MCMCdraws

            TID       = parid;

            % parfor preps (better to do inside parfor loop)
            thisData               = jumpoffData;

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

            % draw shocks
            zdraws      = randn(rndStream, N, irfHorizon, irfNdraws);
            fcstSVdraws = randn(rndStream, N, irfHorizon * irfNdraws);

            % map into state space
            PAIactual                     = PAIdraws(:,ndxYIELDLAGS,mm);
            PAIactual(~actualrateBlock,:) = 0;

            PAIshadow = PAIdraws(:,:,mm); %#ok<PFBNS>
            PAIshadow(actualrateBlock,ndxYIELDLAGS) = 0;

            fcstA               = fcstA0;
            fcstA(ndxfcstY, 1:K)         = PAIshadow;
            fcstA(ndxfcstY, K+1:Nstates) = PAIactual;

            PHI     = ivech(PHIdraws(:,mm), ndxvech);

            SV0         = SVjumpoffDraws(:,mm);
            invA        = invAdraws(:,:,mm);

            % generate SV paths
            sqrtPHI     = chol(PHI, 'lower');
            fcstSVdraws = sqrtPHI * fcstSVdraws;
            fcstSVdraws = reshape(fcstSVdraws, N, irfHorizon, irfNdraws);
            fcstSVdraws = cumsum(fcstSVdraws,2);
            fcstSVdraws = exp(fcstSVdraws * 0.5);

            shock11 = IRF1scale * shocksize;

            simfun = @(nushocks) simVARshadowrateBlockHybrid(nushocks, N, fcstX0, ndxfcstY, ndxfcstActual, ndxfcstShadow, cumcode, np, fcstA, fcstB, invA, ELBbound, ndxYIELDS, irfHorizon, irfNdraws);

            % baseline
            fcstYHATdraws(:,:,mm)       = antitheticSim(0, zdraws, fcstSVdraws, SV0, N, irfHorizon, irfNdraws, simfun);
            % positive shock
            fcstYHATdraws1plus(:,:,mm)  = antitheticSim(shock11, zdraws, fcstSVdraws, SV0, N, irfHorizon, irfNdraws, simfun);
            % negative shock
            fcstYHATdraws1minus(:,:,mm) = antitheticSim(-1 * shock11, zdraws, fcstSVdraws, SV0, N, irfHorizon, irfNdraws, simfun);
            
        end % parfor

        %% clean up workspace
        clear PAIdraws invAdraws PHIdraws shadowrateJumpoffDraws SVjumpoffDraws
        clear fcstA fcstB


        %% IRF
        IRFdraws1plus   = fcstYHATdraws1plus   - fcstYHATdraws;
        IRFdraws1minus  = fcstYHATdraws1minus  - fcstYHATdraws;

        % integrate over MCMC nodes
        fcstYhat        = median(fcstYHATdraws, 3);
        fcstYhat1plus   = median(fcstYHATdraws1plus, 3);
        fcstYhat1minus  = median(fcstYHATdraws1minus, 3);

        IRF1plus        = median(IRFdraws1plus, 3);
        IRF1plusTails   = prctile(IRFdraws1plus, prc70, 3);

        IRF1minus       = median(IRFdraws1minus, 3);
        IRF1minusTails  = prctile(IRFdraws1minus, prc70, 3);

        clear IRFdraws1plus IRFdraws1minus

        %% PLOT RESULTS
        if doPlots
            colorPlus     = Colors4Plots(1);
            colorMinus    = Colors4Plots(2);
            colorBase     = Colors4Plots(8);
            colorRB       = Colors4Plots('green');

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

        end

        %% wrap up
        allw = whos;
        ndx = contains({allw.class}, 'Figure');
        if any(ndx)
            clear(allw(ndx).name)
        end
        save(sprintf('irf2BlockHybrid-%s.mat', titlename), '-v7.3')

        close all
        finishwrap

    end % irfDate
end % irfscale

finishscript


%% define forecast simulation as function
function ydraws = simVARshadowrateBlockHybrid(nushocks, N, fcstX0, ndxfcstY, ndxfcstActual, ndxfcstShadow, cumcode, np, fcstA, fcstB, invA, ELBbound, ndxYIELDS, irfHorizon, irfNdraws)

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

%% wrapper function for antithetic simulation
function yhatdraws = antitheticSim(shock11, zdraws, SVdraws, SV0, N, irfHorizon, irfNdraws, simfun)

ydraws  = NaN(N,irfHorizon,irfNdraws * 4);

mcndxpp = 1: irfNdraws;
mcndxpm = irfNdraws + (1: irfNdraws);
mcndxmp = 2 * irfNdraws + (1: irfNdraws);
mcndxmm = 3 * irfNdraws + (1: irfNdraws);

%% pp
nushocks            = zdraws .* SVdraws .* SV0;
nushocks(1,1,:)     = shock11 + nushocks(1,1,:);
ydraws(:,:,mcndxpp) = simfun(nushocks);

%% pm
nushocks            = -1 .* zdraws .* SVdraws .* SV0;
nushocks(1,1,:)     = shock11 + nushocks(1,1,:);
ydraws(:,:,mcndxpm) = simfun(nushocks);

%% mp
nushocks            = zdraws ./ SVdraws .* SV0;
nushocks(1,1,:)     = shock11 + nushocks(1,1,:);
ydraws(:,:,mcndxmp) = simfun(nushocks);
%% mm
nushocks            = -1 .* zdraws ./ SVdraws .* SV0;
nushocks(1,1,:)     = shock11 + nushocks(1,1,:);
ydraws(:,:,mcndxmm) = simfun(nushocks);

yhatdraws = mean(ydraws,3);

end % function antitheticSim