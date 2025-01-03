function [PAI_all,  hRHO_all, hBAR_all, PHI_all, invA_all, sqrtht_all, shadowrate_all, missingrate_all, ...
    BETA_all, ...
    fcstYdraws, fcstYhat, ...
    fcstShadowrateDraws, fcstShadowrateHat ...
    ] = mcmcVARshadowrateGeneralAR1SV(thisT, MCMCdraws, ...
    p, np, data0, ydates0, ...
    ~, ...
    minnesotaPriorMean, doRATSprior, ...
    ndxSHADOWRATE, ndxOTHERYIELDS, doELBsampling, doELBsampleAlternate, ELBbound, elbT0, check_stationarity, ...
    ~, ~, ...
    ~, fcstNdraws, fcstNhorizons, rndStream, doprogress)

% empty argument is for actualrateWeight

if nargin < 23
    doprogress = false;
end

if check_stationarity
    warning('no stationarity check for general model')
end

%% get TID
% used to provide context for warning messages
TID   = parid;

doPredictiveDensity = nargout > 9;

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
theta5 = 1; % hyperparameter of prior on shadow rates

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

Nshadowrates  = length(ndxSHADOWRATE);
FFRdata       = data(:,ndxSHADOWRATE);

FFRlags=zeros(Nobs,p * Nshadowrates);
for l=1:p
    FFRlags(p+1:Nobs,(l-1)*Nshadowrates + (1 : Nshadowrates)) = FFRdata(p+1-l:Nobs-l,:);
end
Xffrlags = FFRlags(p+1:Nobs,:);
Xffrlags(Xffrlags < ELBbound) = ELBbound;

Xstar = zeros(Nobs-p,Nshadowrates);
X     = [ones(Nobs-p,1) lags(p+1:Nobs,:) Xffrlags];

% trim Y
Y = data(p+1:end,:);
% update pointers
[T,K]=size(X);
Kshadow       = K - p * Nshadowrates; % Number of standard BVAR regressors
Kshadowlag    = Kshadow - 1 ;         % number of standard (endogenous) regressors (i.e. w/o intercept)
% ndxKshadowlag = 1 + (1 : Kshadowlag); % location of the Klag regressors in X, only used to assess stability of partitioned VAR transition


% store original data matrices
Y0 = Y;
X0 = X;
% generate state vector for forecast jumpoff
ndxYIELDS     = union(ndxSHADOWRATE, ndxOTHERYIELDS);
% ndxYIELDLAGS  = cat(2, false, repmat(ismember(1:N, ndxYIELDS), 1, p)); % prepend by false for CONST
% ndxSHADOWRATELAGS = cat(2, false, repmat(ismember(1:N, ndxSHADOWRATE), 1, p)); % used for shadow-rate sampling
Nyields       = length(ndxYIELDS);

% check that shadow rates are ordered on top:
if any(ndxSHADOWRATE > Nshadowrates)
    error('shadow rates not ordered on top')
end
if any(ndxYIELDS > Nshadowrates + Nyields)
    error('yields should be ordered after shadow rates')
end

Nmacro        = N - Nyields;
ndxMACRO      = Nyields + 1 : N;

% Number of states for the forecasting state space
% in principle, Nstates should be larger than K since forecast system needs to track all yields
% (not only those constrained inside the estimation sample)
% however, absent identification of coefficients on those yields, will only track actual rates for which we have ELB data in estimation sample
Nstates       = 1 + (N + Nshadowrates) * p;
if Nstates ~= K
    error('something off with indices')
end
Xjumpoff      = NaN(Nstates,1);
Xjumpoff(1)   = 1;
for l=1:p
    Xjumpoff(1+(l-1)*N+(1:N)) = data(Nobs-(l-1),1:N);
end
% add lagged actual rates to jumpoff
for l=1:p
    Xjumpoff(Kshadow+(l-1)*Nshadowrates+(1:Nshadowrates)) = data(Nobs-(l-1),ndxSHADOWRATE);
end
XjumpoffActualYieldLags                                     = Xjumpoff(Kshadow+1:end);
XjumpoffActualYieldLags(XjumpoffActualYieldLags < ELBbound) = ELBbound;

%% prepare some objects for logscore evalation
doLogscores = false; % doPredictiveDensity && ~isempty(yrealized);
if doLogscores
    
    error('logscore evaluation not implemented for general model')
    
    % macro and yields part of Y vector
    ndxYx      = ~ismember(1:N,ndxYIELDS);
    ndxYi      = ~ndxYx;
    Yxrealized = yrealized(ndxYx,1);
    Yirealized = yrealized(ndxYi,1);
    
    yNatELB    = sum(Yirealized <= ELBbound);
end

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


%% prepare state space for forecasting
% note: there are some redundacies between the forecast state space and the elb state space
% the redundancies allow, however, for easier comparison and maintenance of codes for different model version
fcstA                  = zeros(Nstates,Nstates);
fcstA(1,1)             = 1; % unit root for the constant
% companion for conventional lags TODO CHECK
fcstA(1+N+1:Kshadow,2:Kshadow)     = [eye(N*(p-1)),zeros(N*(p-1),N)];
% companion for actual rates
fcstA(Kshadow+Nshadowrates+(1:Nshadowrates*(p-1)),Kshadow+(1:Nshadowrates*(p-1))) = eye(Nshadowrates*(p-1));
ndxfcstActual          = Kshadow+(1:Nshadowrates);
ndxfcstShadow          = 1+ndxSHADOWRATE;
ndxfcstMacro           = 1+ndxMACRO;

ndxfcstY          = 1+(1:N);
fcstB             = zeros(Nstates,N);
fcstB(ndxfcstY,:) = eye(N);


%% EM, new: prepare ELB data
Ydata                     = data;
shadowrateData     = data(:,ndxSHADOWRATE);
ndxELB                   = shadowrateData <= ELBbound;
shadowrateData(ndxELB)   = NaN;
Ydata(:,ndxSHADOWRATE)   = shadowrateData;

% yNaN
YdataNaN                  = isnan(Ydata);
Ydata(YdataNaN)       = 0;
yNaN                          = YdataNaN(p+1:end,:);
ELBdummy                = double(ndxELB(p+1:end,:));

% EM: prepare state space for shadowrate sampling
% state space with K states, N observables and N shocks
% transition matrix elb.A,  shock loadings elb.B (ex SV)
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
    elb.A                  = zeros(Kshadow,Kshadow,elbT);
    elb.A(1,1,:)             = 1; % unit root for the constant
    for tt = 1 : elbT
        elb.A(1+N+(1:N*(p-1)),1+(1:N*(p-1)),tt) = eye(N*(p-1)); % fill in lower part of companion form
    end
    
    elb.B            = zeros(Kshadow,N,elbT);
    
    % fixed initial conditions for ELB sampler
    elb.X0  = X(elbT0+1,1:Kshadow)'; % time zero values of state.
    elbY0   = reshape(elb.X0(2:end), N, p);
    
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

% shadowrateWeight   = 1 - actualrateWeight;

Pi_pm=zeros(N * Kshadowlag,1); Pi_pv=eye(N * Kshadowlag); co=0;
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
            
            % adjust SHADOWRATE prior for presence of actual rate
            % if any(ismember(ndxSHADOWRATE, j))
            %     Pi_pv(co,co) = shadowrateWeight(j) * Pi_pv(co,co); % scaling variance by half works since prior sees no correlation between coefficients
            %     Pi_pm(co)    = shadowrateWeight(j) * Pi_pm(co);
            % end
        end
    end
end

sigma_FFRlags = NaN(p*Nshadowrates,N);
mean_FFRlags  = zeros(p*Nshadowrates,N);
for i=1:N
    for l=1:p; %#ok<*NOSEL>
        for s=1:Nshadowrates
            j = ndxSHADOWRATE(s);
            if (i==j)
                % if (l==1)
                %     mean_FFRlags((l-1)*Nshadowrates+s,i) = minnesotaPriorMean(i);
                % end
                sigma_FFRlags((l-1)*Nshadowrates+s,i) = theta(1)/(l^theta(4));
            else
                sigma_FFRlags((l-1)*Nshadowrates+s,i) = (AR_s2(i,i)/AR_s2(j,j))*theta(1)*theta(2)/(l^theta(4));
            end
        end
    end
end
% again, scale by weight
% mean_FFRlags  = mean_FFRlags  .* actualrateWeight';
% sigma_FFRlags = sigma_FFRlags .* actualrateWeight';

% Pai~N(vec(MU_pai),OMEGA_pai)
OMEGA_pai   = diag(vec([sigma_const; reshape(diag(Pi_pv),Kshadowlag,N); sigma_FFRlags])); % prior variance of Pai
MU_pai      = [zeros(1,N);reshape(Pi_pm,Kshadowlag,N);mean_FFRlags];                   % prior mean of Pai

% RHO prio
hRHO_mean    = 0.8 .* ones(N,1);
hRHO_V0i     = (1 / 0.2^2) * eye(N);

Nbeta       = Nshadowrates * Nmacro;
OMEGA_beta  = NaN(Nshadowrates, Nmacro);
for i = 1 : Nmacro
    for j = 1 : Nshadowrates
        OMEGA_beta(j,i) = theta(1) * theta(2) * theta5 * AR_s2(Nyields + i, Nyields +i) / AR_s2(j,j);
    end
end
OMEGA_beta  = diag(OMEGA_beta); % prior variance of beta
MU_beta     = zeros(Nbeta,1); % prior mean of beta

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
BETA_all    = NaN(MCMCdraws,Nshadowrates,Nmacro);

shadowrate_all  = NaN(MCMCdraws,Nshadowrates,elbT); % note: will be empty if elbT < 1
missingrate_all = NaN(MCMCdraws,Nshadowrates,elbT);

% define some useful matrices prior to the MCMC loop
% PAI  = zeros(K,N);                                          % pre-allocate space for PAI
% comp = [eye(N*(p-1)),zeros(N*(p-1),N)];                       % companion form
iV     = diag(1./diag(OMEGA_pai)); iVb_prior=iV*vec(MU_pai);    % inverses of prior matrices
iVbeta = diag(1./diag(OMEGA_beta)); iVbetamean=iVbeta*MU_beta; % inverses of prior matrices
EYEn = eye(N);

%% start of MCMC loop

if doprogress
    progressbar(0);
end
m = 0;
maxShake               = 100; % AR-SV draws
countELBaccept         = 0;
% countELBacceptBurnin = 0;
% stackAccept = NaN(MCMCdraws,1);
while m < MCMCreps % using while, not for loop to allow going back in MCMC chain
    
    if m == 0
        % initializations
        PREVdraw.Y          = Y0;                                      % data matrix Y
        PREVdraw.X          = X0;                                      % data matrix X
        PREVdraw.Xstar      = Xstar;
        PREVdraw.A_         = eye(N);
        warning('off', 'MATLAB:rankDeficientMatrix')
        PREVdraw.PAI        = X0\Y0;
        warning( 'on', 'MATLAB:rankDeficientMatrix')
        PREVdraw.BETA       = zeros(Nshadowrates,Nmacro);
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
    Xstar      = PREVdraw.Xstar;
    A_         = PREVdraw.A_;
    sqrtht     = PREVdraw.sqrtht;
    Vol_states = PREVdraw.Vol_states;
    PAI        = PREVdraw.PAI;
    sqrtPHI_   = PREVdraw.sqrtPHI_;
    hRHO       = PREVdraw.hRHO;
    BETA       = PREVdraw.BETA;
    
    
    %%  Draw from the conditional posterior of PAI
    [PAI,BETA] = CTAshadowratebeta(Y,X,Xstar,N,K,T,Nshadowrates,Nmacro,A_,sqrtht,iV,iVb_prior,iVbeta,iVbetamean,PAI,BETA,rndStream);
    RESID = Y - X*PAI; % compute the new residuals
    % account for BETA
    RESID(:,Nyields+1:end) = RESID(:,Nyields+1:end) - Xstar * BETA;
    
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
        PAIactual  = PAI(Kshadow+1:end,:);
        % Yhatactual = X(elbT0+1:end,Kshadow+1:end) * PAIactual;
        Yhatactual = transpose(Xffrlags(elbT0+1:end,:) * PAIactual);
        
        
        
        % map SV into sqrtSigma
        elb.sqrtSigma = sqrtht(elbT0+1:end,:)';
        
        
        elbYgibbs           = elb.Y;
        elb.Y(elb.yNaN)     = NaN; % missing values
        
        %% PS setup
        PAIshadow  = PAI(1:Kshadow,:);
        pai0       = transpose(PAIshadow(1,:)) + Yhatactual;
        pai3       = reshape(transpose(PAIshadow(2:end,:)), N, N, p);
        invbbb     = A_ ./ permute(elb.sqrtSigma, [1 3 2]);
        
        %% construct beta adjustment
        BETAt = zeros(Nmacro,Nshadowrates,elbT);
        for tt = 1 : elbT
            BETAt(:,:,tt) = transpose(BETA) .* ELBdummy(tt,:);
        end
        GAMMA = zeros(N,N,elbT);
        for nn = 1 : N
            GAMMA(nn,nn,:) = 1;
        end
        for tt = 1 : elbT
            GAMMA(ndxMACRO,ndxSHADOWRATE,tt) = BETAt(:,:,tt);
        end
        PAI3  = pagemtimes(reshape(GAMMA, N, N, 1, elbT),pai3);
        for tt = 1 : elbT
            pai0(ndxMACRO,tt) = pai0(ndxMACRO,tt) + BETAt(:,:,tt) * (pai0(ndxSHADOWRATE,tt) - ELBbound);
        end
        
        % invGAMMA
        invGAMMA = GAMMA;
        for tt = 1 : elbT
            invGAMMA(ndxMACRO,ndxSHADOWRATE,tt) = -1 .* BETAt(:,:,tt);
        end
        invBBB = pagemtimes(invbbb,invGAMMA);
        
        % elbA and elbB for Geweke-Sampler
        for tt = 1 : elbT
            elb.A(elb.ndxY,:,tt) = GAMMA(:,:,tt) * transpose(PAI(1:Kshadow,:));
            elb.B(elb.ndxY,:,tt) = GAMMA(:,:,tt) * invA_;
        end
        
        if m == 1
            [~, CC, QQ, RR1, arows, acols, asortndx, brows, bcols, bsortndx] = VARTVPSVprecisionsamplerNaN(PAI3,invBBB,elb.Y,elb.yNaN,elbY0,pai0,rndStream);
        end
        
        % b) reconstruct Y and X
        shadowYdata = Ydata;
        if doELBsampling
            
            
            YYdraws = VARTVPSVprecisionsamplerNaN(PAI3,invBBB,elb.Y,elb.yNaN,...
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
                countELBaccept = countELBaccept + 1;
                % stackAccept(m-MCMCburnin) = ndxAccept;
            else
                % warning('m=%d: proposal did not work, resorting to Gibbs', m)
                shadowrate = gibbsdrawShadowratesTVP(elbYgibbs, elb.X0, Yhatactual, elb.ndxS, elb.sNaN, p, elb.A, elb.B, elb.sqrtSigma, ...
                    elb.bound, 1, elb.gibbsburn, rndStream);
            end
            
            
            missingrate = shadowrateProposals(:,:,1);
            
            % if doELBsampleAlternate
            %     YYdraws = VARTVPSVprecisionsamplerNaN(pai3,invbbb,elb.Y,elb.yNaN,...
            %         elbY0,pai0,rndStream, ...
            %         CC, QQ, RR1, arows, acols, asortndx, brows, bcols, bsortndx);
            %     YYdraws = reshape(YYdraws, N, elbT);
            %     missingrate = YYdraws(ndxSHADOWRATE,:);
            % else
            %     missingrate = NaN;
            % end
            
            shadowYdata(p+elbT0+1:end,ndxSHADOWRATE) = shadowrate';
        else
            
            
            YYdraws = VARTVPSVprecisionsamplerNaN(PAI3,invBBB,elb.Y,elb.yNaN,...
                elbY0,pai0,rndStream, ...
                CC, QQ, RR1, arows, acols, asortndx, brows, bcols, bsortndx);
            YYdraws = reshape(YYdraws, N, elbT);
            missingrate = YYdraws(ndxSHADOWRATE,:);
            
            
            
            if doELBsampleAlternate
                shadowrate = NaN;
            else
                shadowrate = gibbsdrawShadowratesTVP(elb.Y, elb.X0, Yhatactual, elb.ndxS, elb.sNaN, p, elb.A, elb.B, elb.sqrtSigma, ...
                    elb.bound, 1, elb.gibbsburn, rndStream);
            end
            
            shadowYdata(p+elbT0+1:end,ndxSHADOWRATE) = missingrate';
        end
        
        % matrix X
        lags=zeros(Nobs,N*p);
        for l=1:p
            lags(p+1:Nobs,(N*(l-1)+1):N*l) = shadowYdata(p+1-l:Nobs-l,1:N);
        end
        X = [ones(Nobs-p,1) lags(p+1:Nobs,:) Xffrlags];
        
        % trim Y
        Y     = shadowYdata(p+1:end,:);
        Xstar = ELBdummy .* (Y(:,ndxSHADOWRATE) - ELBbound);
        
        % generate state vector for forecast jumpoff
        Xjumpoff      = NaN(Nstates,1);
        Xjumpoff(1)   = 1;
        for l=1:p
            Xjumpoff(1+(l-1)*N+(1:N)) = Y(end-(l-1),1:N);
        end
        Xjumpoff(Kshadow+1:end) = XjumpoffActualYieldLags;
        
        
    end % hasELBdata
    
    %% post burnin: store draws and draw from oos-predictive density
    if m > MCMCburnin
        
        thiniter = thiniter + 1;
        
        if thiniter == MCMCthinstride
            
            thiniter     = 0;
            
            thisMCMCdraw = thisMCMCdraw + 1;
            
            % STORE DRAWS
            PAI_all(thisMCMCdraw,:,:)        = PAI;
            BETA_all(thisMCMCdraw,:,:)       = BETA;
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
                    
                    % note: tracks only actual rates with ELB data in estimation window
                    fcstA(ndxfcstY, :)         = PAI';
                    
                    %% d) predictive logscore (one step ahead)
                    % if doLogscores
                    %
                    %     % % init X0
                    %     % fcstX0              = Xjumpoff;
                    %     %
                    %     % % predictive one-step-ahead moments
                    %     % muX          = fcstA * fcstX0;
                    %     % muY          = muX(ndxfcstY);
                    %     % sqrtOmegaY   = invA_ * diag(fcstSVdraws(:,1,nn));
                    %     %
                    %     % % full vector Y
                    %     % if yNatELB > 0
                    %     %     fcstLogscoreDraws(nn,thisMCMCdraw) = logscoreGaussianCensored(muY, sqrtOmegaY, yrealized(:,1), ELBbound, ndxYIELDS);
                    %     % else
                    %     %     logdetOmegaY = sum(logSV(:,1,nn));
                    %     %     fcstLogscoreDraws(nn,thisMCMCdraw) = logscoreGaussian(muY, sqrtOmegaY, yrealized(:,1),logdetOmegaY);
                    %     % end
                    %     % % separate scores for X and I (macro and yields)
                    %     % % X
                    %     % muYx          = muY(ndxYx);
                    %     % sqrtOmegaYx   = chol(sqrtOmegaY(ndxYx,:) * sqrtOmegaY(ndxYx,:)', 'lower');
                    %     % logdetOmegaYx = 2 * sum(log(diag(sqrtOmegaYx)));
                    %     % fcstLogscoreXdraws(nn,thisMCMCdraw) = logscoreGaussian(muYx, sqrtOmegaYx, Yxrealized, logdetOmegaYx);
                    %     % % I
                    %     % muYi          = muY(ndxYi);
                    %     % sqrtOmegaYi   = chol(sqrtOmegaY(ndxYi,:) * sqrtOmegaY(ndxYi,:)', 'lower');
                    %     % if yNatELB > 0
                    %     %     fcstLogscoreIdraws(nn,thisMCMCdraw) = logscoreGaussianCensored(muYi, sqrtOmegaYi, Yirealized, ELBbound);
                    %     % else
                    %     %     fcstLogscoreIdraws(nn,thisMCMCdraw) = logscoreGaussian(muYi, sqrtOmegaYi, Yirealized);
                    %     % end
                    % end % doLogscores
                    
                    %% forecast simulation
                    fcstX0      = Xjumpoff;
                    theseShocks = invA_ * nushocks(:,:,nn);
                    
                    for hh = 1 : fcstNhorizons
                        fcstXdraw               = fcstA * fcstX0 + fcstB * theseShocks(:,hh);
                        shadowStar              = fcstXdraw(ndxSHADOWRATE) - ELBbound;
                        atELB                   = shadowStar < 0;
                        shadowStar(~atELB)      = 0;
                        fcstXdraw(ndxfcstMacro) = fcstXdraw(ndxfcstMacro) + transpose(BETA) * shadowStar;
                        
                        ydraw        = fcstXdraw(ndxfcstY);
                        
                        % collect draw
                        fcstYdraws(:,hh,nn,thisMCMCdraw) = ydraw;
                        
                        % prepare next iteration
                        fcstX0                   = fcstXdraw;
                        fcstX0(ndxfcstActual)    = max(fcstX0(ndxfcstShadow), ELBbound);
                        
                    end
                    
                end % nn
                
                
            end % doPredictiveDensity
        end % thiniter
    end % if > burnin
    
    
    %% store current draw into PREVdraw
    PREVdraw.Y          = Y;
    PREVdraw.X          = X;
    PREVdraw.Xstar      = Xstar;
    PREVdraw.X0         = Xjumpoff;
    PREVdraw.A_         = A_;
    PREVdraw.sqrtht     = sqrtht;
    PREVdraw.Vol_states = Vol_states;
    PREVdraw.PAI        = PAI;
    PREVdraw.BETA       = BETA;
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

if ~isempty(fcstYdraws)
    fcstYdraws = reshape(fcstYdraws, N, fcstNhorizons, fcstNdraws);
end
fcstYhat   = mean(fcstYdraws,3);

% if doLogscores
%     fcstLogscoreDraws          = reshape(fcstLogscoreDraws, fcstNdraws, 1);
%     fcstLogscoreXdraws         = reshape(fcstLogscoreXdraws, fcstNdraws, 1);
%     fcstLogscoreIdraws         = reshape(fcstLogscoreIdraws, fcstNdraws, 1);
% end

% if doELBsampling
%     fprintf('DONE with thisT %d, TID %d (countELBaccept = %d / %d, countELBacceptBurnin = %d / %d ) \n', thisT, TID, countELBaccept, MCMCdraws, countELBacceptBurnin, MCMCburnin)
% else
%     fprintf('DONE with thisT %d, TID %d \n', thisT, TID)
% end
fprintf('DONE with thisT %d, TID %d \n', thisT, TID)

return

