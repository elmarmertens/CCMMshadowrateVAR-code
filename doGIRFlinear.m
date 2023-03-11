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

doIRF               = true;
irfSCALES           = [1 10];
irfDATES            = [datenum(2007,1,1) datenum(2009,1,1) datenum([2010 2012 2014],12,1)];

p                   = 12;                    % Number of lags on dependent variables
MCMCdraws           = 1e3;                   % Final number of MCMC draws after burn in
irfNdraws           = 1e3;                   % per MCMC node

% SED-PARAMETERS-HERE

irfHorizon          = 25;
fcstNhorizons       = irfHorizon;                 % irrelevant here
fcstNdraws          = MCMCdraws;             % irrelevant here


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
setMinnesotaMean

fontsize = 12;

TID   = parid;

thisT = find(ydates == jumpDate);
T     = thisT - p;

setQuantiles    = [.5, 2.5, 5, normcdf(-1) * 100, 25 , 75,  (1 - normcdf(-1)) * 100, 95, 97.5, 99.5];
Nquantiles      = length(setQuantiles);
ndxCI68         = ismember(setQuantiles, [normcdf(-1) * 100, 100 - normcdf(-1) * 100]);
ndxCI90         = ismember(setQuantiles, [5 95]);
ndxCI           = ndxCI68 | ndxCI90;

%% MCMC sampler

rndStream   = getDefaultStream;

ELBbound = .25; % legacy argument, not relevant for GIRF
[PAI_all, PHI_all, invA_all, sqrtht_all, ...
    ] = mcmcVAR(thisT, MCMCdraws,...
    p, np, data, ydates, ...
    minnesotaPriorMean, doRATSprior, ...
    ndxYIELDS, ELBbound, check_stationarity, ...
    [], ... % yrealized
    fcstNdraws, fcstNhorizons, rndStream, true);

%% store MCMC
titlename=sprintf('%s-p%d-jumpoff%s', datalabel, p, datestr(jumpDate, 'yyyymmm'));

allw = whos;
ndx = contains({allw.class}, 'Figure');
if any(ndx)
    clear(allw(ndx).name)
end
save(sprintf('mcmcLinear-%s.mat', titlename), '-v7.3')


%% collect MCMC draws
PAIdraws  = permute(PAI_all, [3 2 1]);
invAdraws = permute(invA_all, [2 3 1]);
PHIdraws  = permute(PHI_all, [2 1]);
[~, ndxvech]    = ivech(PHIdraws(:,1));
SVdraws   = sqrtht_all;

clear PAI_all invA_all PHI_all sqrtht_all

%% shock draws
zdraws   = randn(rndStream, N, irfHorizon, irfNdraws);
zSVdraws = randn(rndStream, N, irfHorizon * irfNdraws);

%% GIRF
if doIRF
    for IRF1scale = irfSCALES

        display(IRF1scale)

        for irfDate = irfDATES

            display(datestr(irfDate, 'yyyymmm'))


            % prepare wrap
            titlename=sprintf('linear-%s-p%d-IRF1scale%d-jumpoff%s-irfDate%s', datalabel, p, IRF1scale, ...
                datestr(jumpDate, 'yyyymmm'), datestr(irfDate, 'yyyymmm'));

            if ~isempty(samStart)
                titlename = strcat(titlename,'-', datestr(samStart, 'yyyymmm'));
            end

            wrap = [];
            initwrap

            ndxIRFT0               = find(ydates == irfDate);
            SVjumpoffDraws         = permute(SVdraws(:,ndxIRFT0 - p,:), [3 1 2]);
            % allocate memory
            [fcstYdraws, fcstYdraws1plus, fcstYdraws1minus] = deal(NaN(N, irfHorizon, irfNdraws, MCMCdraws));

            % prepare state space for forecast simulation

            ndxfcstY          = 1+(1:N);
            fcstB             = zeros(K,N);
            fcstB(ndxfcstY,:) = eye(N);

            % construct forecast jumpoff (with placeholders for shadow rates)
            jumpoffDate       = ydates(thisT);
            ndx               = ydates <= jumpoffDate;
            jumpoffData       = data(ndx,:);

            prc70 = normcdf([-1 1]) * 100;

            %% loop over MCMC draws
            parfor mm = 1 : MCMCdraws

                TID       = parid;

                % parfor preps (better to do inside parfor loop)
                thisData               = jumpoffData;
                fcstA                  = zeros(K,K);
                fcstA(1,1)             = 1; % unit root for the constant
                fcstA(1+N+1:end,2:end) = [eye(N*(p-1)),zeros(N*(p-1),N)]; % fill in lower part of companion form

                % construct jump off vector
                Xjumpoff    = zeros(K,1);
                Xjumpoff(1) = 1;
                for l=1:p
                    Xjumpoff(1+(l-1)*N+(1:N)) = thisData(ndxIRFT0-(l-1),1:N);
                end
                fcstX0                        = Xjumpoff;

                % map into state space
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
                nushocks(1,1,:)            = IRF1scale;
                fcstYdraws1plus(:,:,:,mm)  = simVAR(N, fcstX0, ndxfcstY, cumcode, np, fcstA, fcstB, invAdraws(:,:,mm), ...
    				nushocks, irfHorizon, irfNdraws);
                % negative shock
                nushocks(1,1,:)            = -1 * IRF1scale;
                fcstYdraws1minus(:,:,:,mm) = simVAR(N, fcstX0, ndxfcstY, cumcode, np, fcstA, fcstB, invAdraws(:,:,mm), ...
    				nushocks, irfHorizon, irfNdraws);

            end % parfor

            % obsolete since these variables are in parfor:
            % clear logSV fcstSVdraws nushocks
            clear SVjumpoffDraws
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
            save(sprintf('irfLinear-%s.mat', titlename), '-v7.3')

            close all
            finishwrap

        end % irfDate
    end % irfscale
end % doIRF

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