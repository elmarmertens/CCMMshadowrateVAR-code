function [PAI_all, hRHO_all, hBAR_all, PHI_all, invA_all, sqrtht_all, shadowrate_all, missingrate_all, ...
    fcstYdraws, fcstYhat, ...
    fcstShadowrateDraws, fcstShadowrateHat, ...
    fcstLogscoreDraws, ...
    fcstLogscoreXdraws, fcstLogscoreIdraws, ...
    IRF1plus, IRF1minus, ...
    IRF1drawsPlus, IRF1drawsMinus, fcstY1hatPlus, fcstY1hatMinus ...
    ] = mcmcVARshadowrateBlockNonstructuralAR1SV(thisT, MCMCdraws, ...
    p, np, data0, ydates0, ...
    actualrateBlock, ...
    minnesotaPriorMean, doRATSprior, ...
    ndxSHADOWRATE, ndxOTHERYIELDS, doELBsampling, doELBsampleAlternate, ELBbound, elbT0, check_stationarity, ...
    IRF1scale, IRFcumcode, ...
    yrealized, fcstNdraws, fcstNhorizons, rndStream, doprogress)

    % stackAccept, ...


if nargin < 23
    doprogress = false;
end

%% get TID
% used to provide context for warning messages
TID   = parid;

doPredictiveDensity = nargout > 8;
doIRF1              = nargout > 15;

%% truncate sample
samEnd = ydates0(thisT);
ndx    = ydates0 <= samEnd;

data   = data0(ndx,:);

%% --------------------------- OPTIONS ----------------------------------
if doRATSprior
    theta=[0.04 0.25 100 2];       % hyperparameters of Minnesota prior:
else
    theta=[0.05 0.5 100 2];       % hyperparameters of Minnesota prior:
end

% [theta1 theta2 int theta3], int is the
% prior on the intercept.
% theta1 is the overall shrinkage, theta2 the
% cross shrinkage and lambda 3 the lag decay
% (quadratic if =2). Note theta2~=1 implies
% the prior becomes asymmetric across eqation,
% so this would not be implementable in the
% standard conjugate setup.
% Minn_pmean = 0;              % Prior mean of the 1-st own lag for each
% equation. For nonstationary variables, this
% is usually set to 1. For transformed
% stationary variables this is set to 0.

MCMCburnin        = 2 * MCMCdraws;
MCMCthinstride    = 1;
MCMCreps          = MCMCthinstride * MCMCdraws + MCMCburnin;        % total MCMC draws
thisMCMCdraw      = 0;
thiniter          = 0;

%% -------------------------Create data matrices-------------------------
% pointers
[Nobs,N]=size(data);
% matrix X
lags=zeros(Nobs,N*p);
for l=1:p
    lags(p+1:Nobs,(N*(l-1)+1):N*l) = data(p+1-l:Nobs-l,1:N);
end
X = [ones(Nobs-p,1) lags(p+1:Nobs,:)];
% trim Y
Y = data(p+1:end,:);
% update pointers
[T,K]=size(X);
Klagreg    = K - 1; % number of lag regressors (without intercept)
ndxKlagreg = 1 + (1 : Klagreg); % location of the Klag regressors in X

% store original data matrices
Y0 = Y;
X0 = X;
% generate state vector for forecast jumpoff
Nshadowrates  = length(ndxSHADOWRATE);
ndxYIELDS     = union(ndxSHADOWRATE, ndxOTHERYIELDS);
ndxYIELDLAGS  = cat(2, false, repmat(ismember(1:N, ndxYIELDS), 1, p)); % prepend by false for CONST
ndxSHADOWRATELAGS = cat(2, false, repmat(ismember(1:N, ndxSHADOWRATE), 1, p)); % used for shadow-rate sampling 
Nyields       = length(ndxYIELDS);
Nstates       = K + Nyields * p; % track lagged actual rates as well
Xjumpoff      = zeros(Nstates,1);
Xjumpoff(1)   = 1;
for l=1:p
    Xjumpoff(1+(l-1)*N+(1:N)) = data(Nobs-(l-1),1:N);
end
% add lagged actual rates to jumpoff
for l=1:p
    Xjumpoff(K+(l-1)*Nyields+(1:Nyields)) = data(Nobs-(l-1),ndxYIELDS);
end


%% prepare some objects for logscore evalation
doLogscores = doPredictiveDensity && ~isempty(yrealized);
if doLogscores
    % macro and yields part of Y vector
    ndxYx      = ~ismember(1:N,ndxYIELDS);
    ndxYi      = ~ndxYx;
    Yxrealized = yrealized(ndxYx,1);
    Yirealized = yrealized(ndxYi,1);

    yNatELB    = sum(Yirealized <= ELBbound);
end

%% create X block for each equation

Xactual = X;
% Xshadow = NaN(T,K);
NactualBlocks = sum(actualrateBlock);
NshadowBlocks = N - NactualBlocks;
XX      = NaN(T,K,N);
XX(:,:,actualrateBlock)  = repmat(Xactual, [1 1 NactualBlocks]);
% XX(:,:,~actualrateBlock) = repmat(Xshadow, [1 1 NshadowBlocks]);


%% allocate memory for out-of-sample forecasts
Ndraws      = fcstNdraws / MCMCdraws;
if mod(fcstNdraws, MCMCdraws) ~= 0
    error('fcstNdraws must be multiple of MCMCdraws')
end
fcstYdraws        = NaN(N,fcstNhorizons, Ndraws, MCMCdraws); % see reshape near end of script
if doLogscores
    fcstLogscoreDraws              = NaN(Ndraws, MCMCdraws); % see reshape near end of script
    [fcstLogscoreXdraws, fcstLogscoreIdraws] = deal(NaN(Ndraws, MCMCdraws)); % see reshape near end of script
else
    [fcstLogscoreDraws, fcstLogscoreXdraws, fcstLogscoreIdraws] = deal([]);
end

if isempty(IRF1scale)
    IRF1scale  = 1;
end

if doIRF1
    [fcstYdraws1plus, fcstYdraws1minus ] = deal(NaN(N,fcstNhorizons, Ndraws, MCMCdraws));
    % do nothing [IRFdraws1plus, IRFdraws1minus ]     = deal(NaN(N,fcstNhorizons, MCMCdraws));
else
    [IRF1plus, IRF1minus]           = deal([]);
    [IRF1drawsPlus, IRF1drawsMinus] = deal([]);
end

%% prepare state space for forecasting
% note: there are some redundacies between the forecast state space and the elb state space
% the redundancies allow, however, for easier comparison and maintenance of codes for different model version
fcstA                  = zeros(Nstates,Nstates);
fcstA(1,1)             = 1; % unit root for the constant
fcstA(1+N+1:K,2:K)     = [eye(N*(p-1)),zeros(N*(p-1),N)]; % fill in lower part of companion form
fcstA(K+Nyields+(1:Nyields*(p-1)),K+(1:Nyields*(p-1))) = eye(Nyields*(p-1)); % companion for actual rates
ndxfcstActual          = K+(1:Nyields);
ndxfcstShadow          = 1+ndxYIELDS;

ndxfcstY          = 1+(1:N);
fcstB             = zeros(Nstates,N);
fcstB(ndxfcstY,:) = eye(N);


%% EM, new: prepare ELB data
Ydata                    = data;
shadowrateData           = data(:,ndxSHADOWRATE);
% actualrateData           = data(:,ndxSHADOWRATE); %only needed for checks
ndxELB                   = shadowrateData <= ELBbound;
shadowrateData(ndxELB)   = NaN;
Ydata(:,ndxSHADOWRATE)   = shadowrateData;

% yNaN
YdataNaN                 = isnan(Ydata);
Ydata(YdataNaN)          = 0;
yNaN                     = YdataNaN(p+1:end,:);


% EM: prepare state space for shadowrate sampling
% state space with K states, N observables and N shocks
% constant transition matrix elb.A, constant shock loadings elb.B (ex SV)
% and time-varying measurement loadings elb.C (to account for missing
% shadow-rate values)

%#ok<*UNRCH>
hasELBdata = T > elbT0;
elbT       = max(0,T - elbT0);

if hasELBdata

    if doELBsampling
        elb.Nproposals     = 1e3;
    end
    elb.gibbsburn      = 1e2;
    elb.ndxS           = ismember(1:N, ndxSHADOWRATE);

    elb.ndxY           = 1+(1:N);
    %     elb.T0             = elbT0;
    %     elb.T              = T - elb.T0;
    elb.yNaN           = yNaN(elbT0+1:end,:)'; % note the transpose, needed for call of sampler later
    elb.sNaN           = elb.yNaN(elb.ndxS,:);

    if any(any(yNaN(1:elbT0,ndxSHADOWRATE)))
        error('something off about elbT0')
    end

    if ~any(any(yNaN(elbT0+1,ndxSHADOWRATE)))
        warning('elbT0 + 1 should be a missing obs, but it is not ...')
    end

    elb.bound = ELBbound;


    % init transition matrix
    elb.A                = zeros(K,K);
    elb.A(1,1)           = 1; % unit root for the constant
    elb.A(1+N+1:end,2:end) = [eye(N*(p-1)),zeros(N*(p-1),N)]; % fill in lower part of companion form

    elb.B            = zeros(K,N);

    % encode "prior" over initial conditions (which are fixed)
    elb.X0  = X(elbT0+1,:)'; % time zero values of state. Recall that X already contains lagged values

    % actualrate = actualrateData(p+elbT0+1:end,:)'; % used for checks later
end

%% -----------------Prior hyperparameters for bvar model

% Prior on conditional mean coefficients, use Minnesota setup
% EM - NOTE: using actual fed funds rate data for Minnesota AR regressions ...
ARresid=NaN(T-1,N);
for i=1:N
    yt_0=[ones(T-1,1) Y(1:end-1,i)];
    yt_1=Y(2:end,i);
    ARresid(:,i)=yt_1-yt_0*(yt_0\yt_1);
end
AR_s2= diag(diag(ARresid'*ARresid))./(T-2);

Pi_pm=zeros(N * Klagreg,1); Pi_pv=eye(N * Klagreg); co=0;
sigma_const = NaN(1,N);
for i=1:N
    sigma_const(i)=AR_s2(i,i)*theta(3);  % this sets the prior variance on the intercept
    for l=1:p; %#ok<*NOSEL>
        for j=1:N
            co=co+1;
            if (i==j)
                if l==1
                    Pi_pm(co)=minnesotaPriorMean(i); % this sets the prior means for the first own lag coefficients.
                end
                Pi_pv(co,co)=theta(1)/(l^theta(4)); % prior variance, own lags
            else
                Pi_pv(co,co)=(AR_s2(i,i)/AR_s2(j,j)*theta(1)*theta(2)/(l^theta(4))); % prior variance, cross-lags
            end
        end
    end
end

% Pai~N(vec(MU_pai),OMEGA_pai)
OMEGA_pai   = diag(vec([sigma_const;reshape(diag(Pi_pv),Klagreg,N)])); % prior variance of Pai
MU_pai      = [zeros(1,N);reshape(Pi_pm,Klagreg,N)];                   % prior mean of Pai

% RHO prio
hRHO_mean    = 0.8 .* ones(N,1);
hRHO_V0i     = (1 / 0.2^2) * eye(N);

% A~N(MU_A,inv(OMEGA_A_inv))
Nareg        = N-1;
MU_A         = NaN(Nareg,N);
OMEGA_A_inv  = NaN(Nareg,Nareg,N);
for ii = 2 : N
    thisNareg = (ii-1);
    % Iim1      = eye(ii-1);
    MU_A(1:thisNareg,ii)                       = 0;
    OMEGA_A_inv(1:thisNareg,1:thisNareg,ii)    = 0;
end;

% PHI~IW(s_PHI,d_PHI)
d_PHI = N+3;                 % prior dofs
s_PHI = d_PHI*(0.15*eye(N)) * 12 / np; % prior scale, where eye(N)=PHI_

% prior on initial states
Vol_0mean    = zeros(N,1);
Vol_0vcvsqrt = 10 * speye(N);

%% PREPARE KSC sampler

[gridKSC, gridKSCt, logy2offset] = getKSC7values(T, N);

%% >>>>>>>>>>>>>>>>>>>>>>>>>> Gibbs sampler <<<<<<<<<<<<<<<<<<<<<<<<<<<

% Storage arrays for posterior draws
PAI_all     = NaN(MCMCdraws,K,N);
PHI_all     = NaN(MCMCdraws,N*(N-1)/2+N);
hRHO_all     = NaN(MCMCdraws,N);
hBAR_all     = NaN(MCMCdraws,N);
invA_all    = NaN(MCMCdraws,N,N);
sqrtht_all  = NaN(MCMCdraws,T,N);

shadowrate_all  = NaN(MCMCdraws,Nshadowrates,elbT); % note: will be empty if elbT < 1
missingrate_all = NaN(MCMCdraws,Nshadowrates,elbT);

% define some useful matrices prior to the MCMC loop
% PAI  = zeros(K,N);                                          % pre-allocate space for PAI
comp = [eye(N*(p-1)),zeros(N*(p-1),N)];                       % companion form
iV   = diag(1./diag(OMEGA_pai)); iVb_prior=iV*vec(MU_pai);    % inverses of prior matrices
EYEn = eye(N);

%% start of MCMC loop

if doprogress
    progressbar(0);
end
m = 0;
maxShake             = 100; % AR-SV draws
countELBaccept       = 0;
countELBacceptBurnin = 0;
stackAccept = NaN(MCMCdraws,1);
while m < MCMCreps % using while, not for loop to allow going back in MCMC chain

    if m == 0
        % initializations
        PREVdraw.Y          = Y0;                                      % data matrix Y
        PREVdraw.X          = X0;                                      % data matrix X
        PREVdraw.A_         = eye(N);
        warning('off', 'MATLAB:rankDeficientMatrix')
        PREVdraw.PAI        = X0\Y0;
        warning( 'on', 'MATLAB:rankDeficientMatrix')
        PREVdraw.sqrtht     = sqrt([ARresid(1,:).^2; ARresid.^2]);     % Initialize sqrt_sqrtht
        PREVdraw.Vol_states = 2*log(PREVdraw.sqrtht);                  % Initialize states
        sqrtPHI_            = sqrt(0.0001)*speye(N);                   % Initialize PHI_, a draw from the covariance matrix W
        PREVdraw.sqrtPHI_   = sqrtPHI_;
        PREVdraw.hRHO       = zeros(N,1);

    end % m == 0
    m = m + 1;


    % init with previous draws values
    Y          = PREVdraw.Y;
    X          = PREVdraw.X;
    A_         = PREVdraw.A_;
    sqrtht     = PREVdraw.sqrtht;
    Vol_states = PREVdraw.Vol_states;
    PAI        = PREVdraw.PAI;
    sqrtPHI_   = PREVdraw.sqrtPHI_;
    hRHO       = PREVdraw.hRHO;


    %%  Draw from the conditional posterior of PAI
    stationary=0;
    while stationary==0;

        % EM: for now, initialize the MCMC with actual fed funds rate data

        % CCM: This is the only new step (triangular algorithm).
        % PAI=triang(Y,X,N,K,T,invA_,sqrtht,iV,iVb_prior,rndStream);

        XX(:,:,~actualrateBlock) = repmat(X, [1 1 NshadowBlocks]);

        PAI=CTAsys(Y,XX,N,K,T,A_,sqrtht,iV,iVb_prior,PAI,rndStream);


        if (check_stationarity==0 || max(abs(eig([PAI(ndxKlagreg,:)' ; comp]))) < 1); stationary = 1; end;
    end
    RESID = NaN(T,N);
    for jj = 1 : N
        RESID(:,jj) = Y(:,jj) - XX(:,:,jj) * PAI(:,jj);
    end

    %%  Draw the covariances
    for ii = 2 : N

        % weighted regression to get Z'Z and Z'z (in Cogley-Sargent 2005 notation)
        y_spread_adj = RESID(:,ii)./sqrtht(:,ii);
        X_spread_adj = RESID(:,1 : ii - 1) ./ sqrtht(:,ii); % note: use of implicit vector expansion
        
        ZZ=X_spread_adj'*X_spread_adj; 
        Zz=X_spread_adj'*y_spread_adj;
        
        % computing posteriors moments
        thisNareg         = (ii-1);
        iValpha_post      = ZZ + OMEGA_A_inv(1:thisNareg,1:thisNareg,ii);
        sqrtiVAlpha_post  = chol(iValpha_post);
        tildealpha_post   = transpose(sqrtiVAlpha_post) \ (Zz + OMEGA_A_inv(1:thisNareg,1:thisNareg,ii) * MU_A(1:thisNareg,ii));
        % draw and store
        alphadraw         =  sqrtiVAlpha_post \ (tildealpha_post + randn(rndStream,thisNareg,1));
        A_(ii,1:ii-1)     = -alphadraw';
    end
    invA_=A_\EYEn; % compute implied draw from A^-1, needed in step 2b.

    %% SV step: Draw mixture states and then volatility states


    logy2 = log((RESID*A_').^2 + logy2offset);

    [Vol_states, hBAR, Vol_shocks, Vol_statesdemeaned] = StochVolKSCcorrsqrtAR1(logy2', Vol_states', hRHO, sqrtPHI_, Vol_0mean, Vol_0vcvsqrt, gridKSC, gridKSCt, N, T, rndStream);
    Vol_states = Vol_states';
    eta        = Vol_shocks';

    sqrtht     = exp(Vol_states/2); %compute sqrtht^0.5 from volatility states, needed in step 2b


    %% Draw volatility variances
    Zdraw     = randn(rndStream, N, T + d_PHI);

    sqrtPHIpost = chol(s_PHI + eta'*eta, 'lower');
    sqrtZZ      = chol(Zdraw * Zdraw'); % note: right uppper choleski
    sqrtPHI_    = sqrtPHIpost / sqrtZZ; % just a square root, not choleski
    PHI_        = sqrtPHI_ * sqrtPHI_'; % derive posterior draw from PHI

    %% DRAW SV persistence
    rhodraws = bayesAR1SURdraw(Vol_statesdemeaned(:,2:end)', Vol_statesdemeaned(:,1:end-1)', PHI_, hRHO_mean, hRHO_V0i, maxShake, rndStream);
    shake    = 0;
    thisOK   = false;
    while ~thisOK && shake < maxShake
        shake     = shake + 1;
        thisOK    = all(abs(rhodraws(:,shake)) < 1);
    end
    if thisOK
        hRHO = rhodraws(:,shake);
    else
        notOKtxt = 'unstableSV';
        error('SV persistence draw failed: %s', notOKtxt)
    end

    %%  shadow rate sampling
    if hasELBdata

        % always apply Geweke Gibbs

        % update state space objects
        elb.Y               = Y(elbT0+1:end,:)';

        % adjust obs for lagged actual rates
        PAIactual                     = PAI(ndxSHADOWRATELAGS,:);
        PAIactual(:,~actualrateBlock) = 0;
        Yhatactual = transpose(Xactual(elbT0+1:end,ndxSHADOWRATELAGS) * PAIactual);

        % update VAR companion form
        PAIshadow                                    = PAI;
        PAIshadow(ndxSHADOWRATELAGS,actualrateBlock) = 0;
        elb.A(elb.ndxY, :)                           = PAIshadow';

        elb.B(elb.ndxY,:)   = invA_;


        % map SV into sqrtSigma
        elb.sqrtSigma = sqrtht(elbT0+1:end,:)';


        elbYgibbs           = elb.Y;
        elb.Y(elb.yNaN)     = NaN; % missing values

        %% PS setup
        pai0    = transpose(PAIshadow(1,:)) + Yhatactual;
        pai3    = reshape(transpose(PAIshadow(2:end,:)), N, N, p);
        invbbb  = A_ ./ permute(elb.sqrtSigma, [1 3 2]);
        elbY0   = reshape(elb.X0(2:end), N, p);
        if m == 1
            [~, CC, QQ, RR1, arows, acols, asortndx, brows, bcols, bsortndx] = VARTVPSVprecisionsamplerNaN(pai3,invbbb,elb.Y,elb.yNaN,elbY0,pai0,rndStream);
        end

        % b) reconstruct Y and X
        shadowYdata = Ydata;
        if doELBsampling

            if m < MCMCburnin * .5
                shadowrate = gibbsdrawShadowrates(elbYgibbs, elb.X0, Yhatactual, elb.ndxS, elb.sNaN, p, elb.A, elb.B, elb.sqrtSigma, ...
                    elb.bound, 1, elb.gibbsburn, rndStream);
            else % try acceptance sampling
                YYdraws = VARTVPSVprecisionsamplerNaN(pai3,invbbb,elb.Y,elb.yNaN,...
                    elbY0,pai0,rndStream, ...
                    CC, QQ, RR1, arows, acols, asortndx, brows, bcols, bsortndx, elb.Nproposals);
                YYdraws = reshape(YYdraws, N, elbT, elb.Nproposals);

                shadowrateProposals = YYdraws(ndxSHADOWRATE,:,:);

                Ok = false;
                ndxAccept = 0;
                while ~Ok && (ndxAccept < elb.Nproposals)
                    ndxAccept = ndxAccept + 1;
                    thisProposal = shadowrateProposals(:,:,ndxAccept);
                    Ok = all(thisProposal(elb.sNaN) < ELBbound);
                end
                if Ok
                    shadowrate     = shadowrateProposals(:,:,ndxAccept);
                    if m > MCMCburnin
                        countELBaccept = countELBaccept + 1;
                        stackAccept(m-MCMCburnin) = ndxAccept;
                    else
                        countELBacceptBurnin = countELBacceptBurnin + 1;
                    end
                else
                    shadowrate = gibbsdrawShadowrates(elbYgibbs, elb.X0, Yhatactual, elb.ndxS, elb.sNaN, p, elb.A, elb.B, elb.sqrtSigma, ...
                        elb.bound, 1, elb.gibbsburn, rndStream);
                    % fprintf('%d, none accepted\n', m);
                end
            end



            if doELBsampleAlternate
                YYdraws = VARTVPSVprecisionsamplerNaN(pai3,invbbb,elb.Y,elb.yNaN,...
                    elbY0,pai0,rndStream, ...
                    CC, QQ, RR1, arows, acols, asortndx, brows, bcols, bsortndx);
                YYdraws = reshape(YYdraws, N, elbT);
                missingrate = YYdraws(ndxSHADOWRATE,:);
            else
                missingrate = NaN;
            end

            % do some checks
            % checkdiff(actualrate(~elb.sNaN), shadowrate(~elb.sNaN), [], 'shadowrate sampling'); % sample shadowrate matrix should collapse to actualrate when above the ELB
            % shadowrate(~elb.sNaN) = actualrate(~elb.sNaN); % @ if desired, this line sets shadowrate equal to actual rate (above ELB)

            shadowYdata(p+elbT0+1:end,ndxSHADOWRATE) = shadowrate';
        else


            YYdraws = VARTVPSVprecisionsamplerNaN(pai3,invbbb,elb.Y,elb.yNaN,...
                elbY0,pai0,rndStream, ...
                CC, QQ, RR1, arows, acols, asortndx, brows, bcols, bsortndx);
            YYdraws = reshape(YYdraws, N, elbT);
            missingrate = YYdraws(ndxSHADOWRATE,:);



            if doELBsampleAlternate
                shadowrate = NaN;
            else
                shadowrate = gibbsdrawShadowrates(elb.Y, elb.X0, Yhatactual, elb.ndxS, elb.sNaN, p, elb.A, elb.B, elb.sqrtSigma, ...
                    elb.bound, 1, elb.gibbsburn, rndStream);
                % do some checks
                % checkdiff(actualrate(~elb.sNaN), shadowrate(~elb.sNaN)); % sample shadowrate matrix should collapse to actualrate when above the ELB
                % shadowrate(~elb.sNaN) = actualrate(~elb.sNaN); % if desired, this line sets shadowrate equal to actual rate (above ELB)

            end

            shadowYdata(p+elbT0+1:end,ndxSHADOWRATE) = missingrate';
        end

        % matrix X
        lags=zeros(Nobs,N*p);
        for l=1:p
            lags(p+1:Nobs,(N*(l-1)+1):N*l) = shadowYdata(p+1-l:Nobs-l,1:N);
        end
        X = [ones(Nobs-p,1) lags(p+1:Nobs,:)];

        % trim Y
        Y = shadowYdata(p+1:end,:);

        % generate state vector for forecast jumpoff
        Xjumpoff     = zeros(Nstates,1);
        Xjumpoff(1) = 1;
        for l=1:p
            Xjumpoff(1+(l-1)*N+(1:N)) = Y(end-(l-1),1:N);
        end
        % add lagged actual rates to jumpoff
        for l=1:p
            Xjumpoff(K+(l-1)*Nyields+(1:Nyields)) = data(Nobs-(l-1),ndxYIELDS);
        end


    end % hasELBdata

    %% post burnin: store draws and draw from oos-predictive density
    if m > MCMCburnin

        thiniter = thiniter + 1;

        if thiniter == MCMCthinstride

            thiniter     = 0;

            thisMCMCdraw = thisMCMCdraw + 1;

            % STORE DRAWS
            PAI_all(thisMCMCdraw,:,:)        = PAI;
            PHI_all(thisMCMCdraw,:)          = PHI_((tril(PHI_))~=0);
            invA_all(thisMCMCdraw,:,:)       = invA_;
            sqrtht_all(thisMCMCdraw,:,:)     = sqrtht;
            hRHO_all(thisMCMCdraw,:)       = hRHO;
            hBAR_all(thisMCMCdraw,:)       = hBAR;
            if hasELBdata
                shadowrate_all(thisMCMCdraw,:,:)  = shadowrate;
                missingrate_all(thisMCMCdraw,:,:) = missingrate;
            end


            %% compute OOS draws
            

            if doPredictiveDensity

                % draw and scale SV shocks
                logSV0      = Vol_states(end,:)'; % Note: Vol_states record logs of *variances*
                [fcstSVdraws, logSVdraws] = simulateSVar1(hRHO, hBAR, logSV0, sqrtPHI_, Ndraws, fcstNhorizons,rndStream);
                fcstSVdraws = permute(fcstSVdraws, [1,3,2]); % N x fcstNhorizons x Ndraws


                % draw random numbers and scale by SV
                nushocks    = fcstSVdraws .* randn(rndStream, N, fcstNhorizons, Ndraws);

                for nn = 1 : Ndraws

            
                    %% update VAR companion form
                    PAIactual                     = PAI(ndxYIELDLAGS,:);
                    PAIactual(:,~actualrateBlock) = 0;

                    PAIshadow = PAI;
                    PAIshadow(ndxYIELDLAGS,actualrateBlock) = 0;

                    fcstA(ndxfcstY, 1:K)         = PAIshadow';
                    fcstA(ndxfcstY, K+1:Nstates) = PAIactual';

                    %% d) predictive logscore (one step ahead)
                    if doLogscores

                        % init X0
                        fcstX0              = Xjumpoff;

                        % predictive one-step-ahead moments 
                        muX          = fcstA * fcstX0;
                        muY          = muX(ndxfcstY);
                        sqrtOmegaY   = invA_ * diag(fcstSVdraws(:,1,nn));

                        % full vector Y
                        if yNatELB > 0
                            fcstLogscoreDraws(nn,thisMCMCdraw) = logscoreGaussianCensored(muY, sqrtOmegaY, yrealized(:,1), ELBbound, ndxYIELDS);
                        else
                            logdetOmegaY = sum(logSVdraws(:,1,nn));
                            fcstLogscoreDraws(nn,thisMCMCdraw) = logscoreGaussian(muY, sqrtOmegaY, yrealized(:,1),logdetOmegaY);
                        end
                        % separate scores for X and I (macro and yields)
                        % X
                        muYx          = muY(ndxYx);
                        sqrtOmegaYx   = chol(sqrtOmegaY(ndxYx,:) * sqrtOmegaY(ndxYx,:)', 'lower');
                        logdetOmegaYx = 2 * sum(log(diag(sqrtOmegaYx)));
                        fcstLogscoreXdraws(nn,thisMCMCdraw) = logscoreGaussian(muYx, sqrtOmegaYx, Yxrealized, logdetOmegaYx);
                        % I
                        muYi          = muY(ndxYi);
                        sqrtOmegaYi   = chol(sqrtOmegaY(ndxYi,:) * sqrtOmegaY(ndxYi,:)', 'lower');
                        if yNatELB > 0
                            fcstLogscoreIdraws(nn,thisMCMCdraw) = logscoreGaussianCensored(muYi, sqrtOmegaYi, Yirealized, ELBbound);
                        else
                            fcstLogscoreIdraws(nn,thisMCMCdraw) = logscoreGaussian(muYi, sqrtOmegaYi, Yirealized);
                        end
                    end % doLogscores

                    %% forecast simulation
                    fcstX0      = Xjumpoff;
                    theseShocks = invA_ * nushocks(:,:,nn);

                    for hh = 1 : fcstNhorizons
                        fcstXdraw    = fcstA * fcstX0 + fcstB * theseShocks(:,hh);
                        ydraw        = fcstXdraw(ndxfcstY);

                        % collect draw
                        fcstYdraws(:,hh,nn,thisMCMCdraw) = ydraw;

                        % prepare next iteration
                        fcstX0                   = fcstXdraw;
                        fcstX0(ndxfcstActual)    = max(fcstX0(ndxfcstShadow), ELBbound);

                    end

                    if doIRF1


                        %% plus simulation
                        nushocks(1,1,nn) = IRF1scale;
                        theseShocks      = invA_ * nushocks(:,:,nn);
                        fcstX0           = Xjumpoff;

                        for hh = 1 : fcstNhorizons
                            fcstXdraw    = fcstA * fcstX0 + fcstB * theseShocks(:,hh);
                            ydraw        = fcstXdraw(ndxfcstY);

                            % collect draw
                            fcstYdraws1plus(:,hh,nn,thisMCMCdraw) = ydraw;

                            % prepare next iteration
                            fcstX0                   = fcstXdraw;
                            fcstX0(ndxfcstActual)    = max(fcstX0(ndxfcstShadow), ELBbound);
                        end

                        %% minus simulation
                        nushocks(1,1,nn) = -1 * IRF1scale;
                        theseShocks      = invA_ * nushocks(:,:,nn);
                        fcstX0           = Xjumpoff;


                        for hh = 1 : fcstNhorizons
                            fcstXdraw    = fcstA * fcstX0 + fcstB * theseShocks(:,hh);
                            ydraw        = fcstXdraw(ndxfcstY);

                            % collect draw
                            fcstYdraws1minus(:,hh,nn,thisMCMCdraw) = ydraw;

                            % prepare next iteration
                            fcstX0                   = fcstXdraw;
                            fcstX0(ndxfcstActual)    = max(fcstX0(ndxfcstShadow), ELBbound);
                        end
                    end

                end % nn


            end % doPredictiveDensity
        end % thiniter
    end % if > burnin


    %% store current draw into PREVdraw
    PREVdraw.Y          = Y;
    PREVdraw.X          = X;
    PREVdraw.X0         = Xjumpoff;
    PREVdraw.A_         = A_;
    PREVdraw.sqrtht     = sqrtht;
    PREVdraw.Vol_states = Vol_states;
    PREVdraw.PAI        = PAI;
    PREVdraw.sqrtPHI_   = sqrtPHI_;
    PREVdraw.hRHO       = hRHO;

    if doprogress
        progressbar(m / MCMCreps)
    end

end %end of the Gibbs sampler


% uncensored forecasts
fcstShadowrateDraws = fcstYdraws(ndxYIELDS,:,:);
fcstShadowrateHat   = mean(fcstShadowrateDraws,3);

% censor the forecast draws
yieldDraws                   = fcstYdraws(ndxYIELDS,:,:,:);
ndx                          = yieldDraws < ELBbound;
yieldDraws(ndx)              = ELBbound;
fcstYdraws(ndxYIELDS,:,:,:)  = yieldDraws;

if doIRF1
    %% censor IRFsims
    % plus
    yieldDraws                        = fcstYdraws1plus(ndxYIELDS,:,:,:);
    ndx                               = yieldDraws < ELBbound;
    yieldDraws(ndx)                   = ELBbound;
    fcstYdraws1plus(ndxYIELDS,:,:,:)  = yieldDraws;

    irfdraws                          = fcstYdraws1plus - fcstYdraws;
    irfdraws(IRFcumcode,:,:,:)        = cumsum(irfdraws(IRFcumcode,:,:,:), 2);
    % note: fcstYdraws is not cumulated inside this function (for legacy reasons)

    IRF1drawsPlus                     = mean(irfdraws, 3);
    IRF1plus                          = mean(IRF1drawsPlus, 4);
    IRF1drawsPlus                     = permute(IRF1drawsPlus, [1 2 4 3]);

    fcstYdraws1plus(IRFcumcode,:,:,:) = cumsum(fcstYdraws1plus(IRFcumcode,:,:,:), 2);
    fcstY1hatPlus                     = mean(fcstYdraws1plus, [3 4]);

    % minus
    yieldDraws                        = fcstYdraws1minus(ndxYIELDS,:,:,:);
    ndx                               = yieldDraws < ELBbound;
    yieldDraws(ndx)                   = ELBbound;
    fcstYdraws1minus(ndxYIELDS,:,:,:) = yieldDraws;

    irfdraws                          = fcstYdraws1minus - fcstYdraws;
    irfdraws(IRFcumcode,:,:,:)        = cumsum(irfdraws(IRFcumcode,:,:,:), 2);
    
    IRF1drawsMinus                    = mean(irfdraws, 3);
    IRF1minus                         = mean(IRF1drawsMinus, 4);
    IRF1drawsMinus                    = permute(IRF1drawsMinus, [1 2 4 3]);

    fcstYdraws1minus(IRFcumcode,:,:,:) = cumsum(fcstYdraws1minus(IRFcumcode,:,:,:), 2);
    fcstY1hatMinus                     = mean(fcstYdraws1minus, [3 4]);

end

if ~isempty(fcstYdraws)
    fcstYdraws = reshape(fcstYdraws, N, fcstNhorizons, fcstNdraws);
end
fcstYhat   = mean(fcstYdraws,3);

if doLogscores
    fcstLogscoreDraws          = reshape(fcstLogscoreDraws, fcstNdraws, 1);
    fcstLogscoreXdraws         = reshape(fcstLogscoreXdraws, fcstNdraws, 1);
    fcstLogscoreIdraws         = reshape(fcstLogscoreIdraws, fcstNdraws, 1);
end

if doELBsampling
    fprintf('DONE with thisT %d, TID %d (countELBaccept = %d / %d, countELBacceptBurnin = %d / %d ) \n', thisT, TID, countELBaccept, MCMCdraws, countELBacceptBurnin, MCMCburnin)
else
    fprintf('DONE with thisT %d, TID %d \n', thisT, TID)
end

return

