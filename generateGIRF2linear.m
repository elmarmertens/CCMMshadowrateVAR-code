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


% truncate start of sample (if desired)
if ~isempty(samStart)
    ndx    = ydates >= samStart;
    data   = data(ndx,:);
    ydates = ydates(ndx);
    Tdata  = length(ydates);
end



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


%% collect MCMC results
titlename=sprintf('%s-p%d-jumpoff%s', datalabel, p, datestr(jumpDate, 'yyyymmm'));

mcmc = matfile(fullfile(resultsdir, sprintf('mcmcLinear-%s', titlename)));

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
        SVjumpoffDraws         = permute(mcmc.sqrtht_all(:,ndxIRFT0 - p,:), [3 1 2]);

        %% allocate memory
        [fcstYHATdraws, fcstYHATdraws1plus, fcstYHATdraws1minus] = deal(NaN(N, irfHorizon, MCMCdraws));

        % prepare state space for forecast simulation
        ndxfcstY          = 1+(1:N);
        fcstB             = zeros(K,N);
        fcstB(ndxfcstY,:) = eye(N);

        fcstA0                  = zeros(K,K);
        fcstA0(1,1)             = 1; % unit root for the constant
        fcstA0(1+N+1:end,2:end) = [eye(N*(p-1)),zeros(N*(p-1),N)]; % fill in lower part of companion form

        % construct forecast jumpoff (with placeholders for shadow rates)
        jumpoffDate       = ydates(thisT);
        ndx               = ydates <= jumpoffDate;
        jumpoffData       = data(ndx,:);

        prc70 = normcdf([-1 1]) * 100;

        %% generate draws RB for linear IRF
        IRF1RBdraws = NaN(N, irfHorizon, MCMCdraws);
        YhatRBdraws = NaN(N, irfHorizon, MCMCdraws);
        parfor mm = 1 : MCMCdraws

            % map into state space
            fcstA               = fcstA0;
            fcstA(ndxfcstY, :)  = PAIdraws(:,:,mm);

            % baseline path
            thisData               = jumpoffData;
            % construct jump off vector
            Xjumpoff    = zeros(K,1);
            Xjumpoff(1) = 1;
            for l=1:p
                Xjumpoff(1+(l-1)*N+(1:N)) = thisData(ndxIRFT0-(l-1),1:N);
            end
            fcstX0                   = Xjumpoff;
            for hh = 1 : irfHorizon
                fcstX0               = fcstA * fcstX0;
                YhatRBdraws(:,hh,mm) = fcstX0(ndxfcstY);
            end

            % IRF
            thisResponse              = zeros(K,1);
            thisResponse(ndxfcstY(1)) = IRF1scale * shocksize;
            thisResponse(ndxfcstY)    = invAdraws(:,:,mm) * thisResponse(ndxfcstY);
            for hh = 1 : irfHorizon
                IRF1RBdraws(:,hh,mm) = thisResponse(ndxfcstY);
                thisResponse         = fcstA * thisResponse;
            end
        end

        % cumulate and collect moments

        YhatRBdraws(cumcode, :,:) = cumsum(YhatRBdraws(cumcode,:,:), 2) / np;
        YhatRB        = median(YhatRBdraws, 3);
        YhatRBtails   = prctile(YhatRBdraws, prc70, 3);

        IRF1RBdraws(cumcode, :,:) = cumsum(IRF1RBdraws(cumcode,:,:), 2) / np;
        IRF1RB        = median(IRF1RBdraws, 3);
        IRF1RBtails   = prctile(IRF1RBdraws, prc70, 3);

        %% GIRF simulation across MCMC nodes
        parfor mm = 1 : MCMCdraws

            TID       = parid;

            % parfor preps (better to do inside parfor loop)
            thisData               = jumpoffData;

            % construct jump off vector
            Xjumpoff    = zeros(K,1);
            Xjumpoff(1) = 1;
            for l=1:p
                Xjumpoff(1+(l-1)*N+(1:N)) = thisData(ndxIRFT0-(l-1),1:N);
            end
            fcstX0                        = Xjumpoff;

            % draw shocks
            zdraws      = randn(rndStream, N, irfHorizon, irfNdraws);
            fcstSVdraws = randn(rndStream, N, irfHorizon * irfNdraws);

            % map into state space
            fcstA               = fcstA0;
            fcstA(ndxfcstY, :)  = PAIdraws(:,:,mm);

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

            simfun = @(nushocks) simVAR(nushocks, N, fcstX0, ndxfcstY, cumcode, np, fcstA, fcstB, invA, irfHorizon, irfNdraws);
            
            % baseline
            fcstYHATdraws(:,:,mm)       = antitheticSim(0, zdraws, fcstSVdraws, SV0, N, irfHorizon, irfNdraws, simfun);
            % positive shock
            fcstYHATdraws1plus(:,:,mm)  = antitheticSim(shock11, zdraws, fcstSVdraws, SV0, N, irfHorizon, irfNdraws, simfun);
            % negative shock
            fcstYHATdraws1minus(:,:,mm) = antitheticSim(-1 * shock11, zdraws, fcstSVdraws, SV0, N, irfHorizon, irfNdraws, simfun);
            
        end % parfor

        %% clean up workspace
        clear PAIdraws invAdraws PHIdraws SVjumpoffDraws
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
                hold on
                set(gca, 'FontSize', fontsize)
                hplus  = plot(0:irfHorizon-1, IRF1plus(n,:), '-', 'color', colorPlus, 'linewidth', 3);
                plot(0:irfHorizon-1, squeeze(IRF1plusTails(n,:,:)), '-', 'color', colorPlus, 'linewidth', 1);
                % hminus = plot(0:irfHorizon-1, -1 * IRF1minus(n,:), '-.', 'color', colorMinus, 'linewidth', 3);
                % plot(0:irfHorizon-1, -1 * squeeze(IRF1minusTails(n,:,:,:)), '-.', 'color', colorMinus, 'linewidth', 1);
                hRB = plot(0:irfHorizon-1, IRF1RB(n,:), '-.', 'color', colorRB, 'linewidth', 3);
                plot(0:irfHorizon-1, squeeze(IRF1RBtails(n,:,:)), '-.', 'color', colorRB, 'linewidth', 1);
                xlim([0 irfHorizon-1])
                yline(0, 'k:')
                legend([hplus, hRB], 'simulated GIRF', 'RB IRF', 'location', 'southoutside')

                sgtitle(sprintf('%s per %s', Ylabels{n}, datestr(irfDate, 'yyyymmm')), 'FontSize', 18', 'FontWeight', 'bold')

                wrapthisfigure(thisfig, sprintf('IRF1rb-%s-IRF1scale%d-%s-jumpoff%s-irfDate%s', datalabel, IRF1scale, ncode{n}, ...
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

                hRB    = plot(0:irfHorizon-1, YhatRB(n,:), '--', 'color', colorRB, 'linewidth', 2);

                xlim([0 irfHorizon-1])

                legend([hbase hplus hminus hRB], 'baseline', 'positive shock', 'negative shock', 'baseline (RB)', 'location', 'best')
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
        save(sprintf('irf2Linear-%s.mat', titlename), '-v7.3')

        close all
        finishwrap

    end % irfDate
end % irfscale

finishscript


%% define forecast simulation as function
function ydraws = simVAR(nushocks, N, fcstX0, ndxfcstY, cumcode, np, fcstA, fcstB, invA, irfHorizon, irfNdraws)

ydraws      = NaN(N,irfHorizon,irfNdraws);
theseShocks = zeros(N, irfHorizon+1); % padded with zeros for use with ltitr

for nn = 1 : irfNdraws
    theseShocks(:,1:irfHorizon)  = invA * nushocks(:,:,nn);
    xdraws         = ltitr(fcstA, fcstB, theseShocks', fcstX0); % faster forecast simulation using ltitr
    ydraws(:,:,nn) = xdraws(2:end,ndxfcstY)';
end

ydraws(cumcode, :,:)         = cumsum(ydraws(cumcode,:,:), 2) / np;

end % function simVARshadowrate

%% wrapper function for antithetic sumulation
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