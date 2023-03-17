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

%% shock draws
zdraws   = randn(rndStream, N, irfHorizon, irfNdraws);
zSVdraws = randn(rndStream, N, irfHorizon * irfNdraws);

%% shock size
switch datalabel
    case {'fredMD20VXO-2022-09', 'fredMD20VXOexYield-2022-09'}
        shocksize = 1.27;
    case {'fredMD20EBP-2022-09', 'fredMD20EBPexYield-2022-09', 'fredMD3EBP-2022-09'}
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
        SVjumpoffDraws         = permute(mcmc.sqrtht_all(:,ndxIRFT0 - p,:), [3 1 2]);

        %% allocate memory
        [fcstYdraws, fcstYdraws1plus, fcstYdraws1minus] = deal(NaN(N, irfHorizon, irfNdraws, MCMCdraws));

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

        %% loop over MCMC draws
        %         progressbar(0)
        parfor mm = 1 : MCMCdraws

            TID       = parid;

            % parfor preps (better to do inside parfor loop)
            thisData               = jumpoffData;
            %             fcstA                  = zeros(K,K);
            %             fcstA(1,1)             = 1; % unit root for the constant
            %             fcstA(1+N+1:end,2:end) = [eye(N*(p-1)),zeros(N*(p-1),N)]; % fill in lower part of companion form

            % construct jump off vector
            Xjumpoff    = zeros(K,1);
            Xjumpoff(1) = 1;
            for l=1:p
                Xjumpoff(1+(l-1)*N+(1:N)) = thisData(ndxIRFT0-(l-1),1:N);
            end
            fcstX0                        = Xjumpoff;

            % map into state space
            fcstA               = fcstA0;
            fcstA(ndxfcstY, :)  = PAIdraws(:,:,mm);

            PHI     = ivech(PHIdraws(:,mm), ndxvech);

            % generate SV paths
            sqrtPHI     = chol(PHI, 'lower');
            logSV       = sqrtPHI * zSVdraws;
            logSV       = reshape(logSV, N, irfHorizon, irfNdraws);
            logSV       = cumsum(logSV,2);
            fcstSVdraws = exp(logSV * 0.5) .* SVjumpoffDraws(:,mm);
            nushocks    = fcstSVdraws .* zdraws;

            % baseline
            fcstYdraws(:,:,:,mm)       = simVAR(N, fcstX0, ndxfcstY, cumcode, np, fcstA, fcstB, invAdraws(:,:,mm), ...
                nushocks, irfHorizon, irfNdraws);
            % positive shock
            nushocks(1,1,:)            = IRF1scale * shocksize;
            fcstYdraws1plus(:,:,:,mm)  = simVAR(N, fcstX0, ndxfcstY, cumcode, np, fcstA, fcstB, invAdraws(:,:,mm), ...
                nushocks, irfHorizon, irfNdraws);
            % negative shock
            nushocks(1,1,:)            = -1 * IRF1scale * shocksize;
            fcstYdraws1minus(:,:,:,mm) = simVAR(N, fcstX0, ndxfcstY, cumcode, np, fcstA, fcstB, invAdraws(:,:,mm), ...
                nushocks, irfHorizon, irfNdraws);

            %             progressbar(mm / MCMCdraws)

        end % parfor

        clear logSV fcstSVdraws nushocks % obsolete in case of parfor
        clear PAIdraws invAdraws PHIdraws SVjumpoffDraws
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

            legend([hbase hplus hminus], 'baseline', 'positive shock', 'negative shock', 'location', 'best')
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
        save(sprintf('irfLinear-%s.mat', titlename), '-v7.3')

        close all
        finishwrap

    end % irfDate
end % irfscale

finishscript


%% define forecast simulation as function
function ydraws = simVAR(N, fcstX0, ndxfcstY, cumcode, np, fcstA, fcstB, invA, nushocks, irfHorizon, irfNdraws)

ydraws      = NaN(N,irfHorizon,irfNdraws);
theseShocks = zeros(N, irfHorizon+1); % padded with zeros for use with ltitr

for nn = 1 : irfNdraws
    theseShocks(:,1:irfHorizon)  = invA * nushocks(:,:,nn);
    xdraws         = ltitr(fcstA, fcstB, theseShocks', fcstX0); % faster forecast simulation using ltitr
    ydraws(:,:,nn) = xdraws(2:end,ndxfcstY)';
end

ydraws(cumcode, :,:)         = cumsum(ydraws(cumcode,:,:), 2) / np;

end % function simVARshadowrate