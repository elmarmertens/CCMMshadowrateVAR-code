function [PAI_all, PHI_all, invA_all, sqrtht_all, shadowrate_all, missingrate_all, ...
    fcstYdraws, fcstYhat, ...
    fcstShadowrateDraws, fcstShadowrateHat, ...
    fcstLogscoreDraws, ...
    fcstLogscoreXdraws, fcstLogscoreIdraws, ...
    IRF1plus, IRF1minus, ...
    IRF1drawsPlus, IRF1drawsMinus, fcstY1hatPlus, fcstY1hatMinus ...
    ] = mcmcVARhybridGibbs(thisT, MCMCdraws, ...
    p, np, data0, ydates0, ...
    ~, ...
    minnesotaPriorMean, doRATSprior, ...
    ndxSHADOWRATE, ndxOTHERYIELDS, doELBsampling, doELBsampleAlternate, ELBbound, elbT0, check_stationarity, ...
    IRF1scale, IRFcumcode, ...
    yrealized, fcstNdraws, fcstNhorizons, rndStream, doprogress)

    % stackAccept, ...
    % empty argument is for actualrateWeight
    
if nargin < 23
    doprogress = false;
end

if check_stationarity
    warning('no stationarity check for hybrid model')
end

%% get TID
% used to provide context for warning messages
TID   = parid;

doPredictiveDensity = nargout > 6;
doIRF1              = nargout > 13;

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

MCMCburnin        = MCMCdraws;
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

X = [ones(Nobs-p,1) lags(p+1:Nobs,:) Xffrlags];

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
% Nyields       = length(ndxYIELDS);
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
doLogscores = doPredictiveDensity && ~isempty(yrealized);
if doLogscores
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
% companion for conventional lags TODO CHECK
fcstA(1+N+1:Kshadow,2:Kshadow)     = [eye(N*(p-1)),zeros(N*(p-1),N)];
% companion for actual rates
fcstA(Kshadow+Nshadowrates+(1:Nshadowrates*(p-1)),Kshadow+(1:Nshadowrates*(p-1))) = eye(Nshadowrates*(p-1)); 
ndxfcstActual          = Kshadow+(1:Nshadowrates);
ndxfcstShadow          = 1+ndxSHADOWRATE;

ndxfcstY          = 1+(1:N);
fcstB             = zeros(Nstates,N);
fcstB(ndxfcstY,:) = eye(N);


%% EM, new: prepare ELB data
Ydata                    = data;
shadowrateData           = data(:,ndxSHADOWRATE);
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
    elb.gibbsburn      = 1e1;
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
    elb.A                  = zeros(Kshadow,Kshadow);
    elb.A(1,1)             = 1; % unit root for the constant
    elb.A(1+N+1:end,2:end) = [eye(N*(p-1)),zeros(N*(p-1),N)]; % fill in lower part of companion form
    
    elb.B            = zeros(Kshadow,N);
    
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
Vol_0mean     = zeros(N,1);   % time 0 states prior mean
Vol_0vcvsqrt  = 10*speye(N); % chol(Vol_0var)';

%% PREPARE KSC sampler

[gridKSC, gridKSCt, logy2offset] = getKSC7values(T, N);

%% >>>>>>>>>>>>>>>>>>>>>>>>>> Gibbs sampler <<<<<<<<<<<<<<<<<<<<<<<<<<<

% Storage arrays for posterior draws
PAI_all     = NaN(MCMCdraws,K,N);
PHI_all     = NaN(MCMCdraws,N*(N-1)/2+N);
invA_all    = NaN(MCMCdraws,N,N);
sqrtht_all  = NaN(MCMCdraws,T,N);

shadowrate_all  = NaN(MCMCdraws,Nshadowrates,elbT); % note: will be empty if elbT < 1
missingrate_all = NaN(MCMCdraws,Nshadowrates,elbT);

% define some useful matrices prior to the MCMC loop
% PAI  = zeros(K,N);                                          % pre-allocate space for PAI
% comp = [eye(N*(p-1)),zeros(N*(p-1),N)];                       % companion form
iV   = diag(1./diag(OMEGA_pai)); iVb_prior=iV*vec(MU_pai);    % inverses of prior matrices
EYEn = eye(N);

%% start of MCMC loop

if doprogress
    progressbar(0);
end
m = 0;
% countELBaccept       = 0;
% countELBacceptBurnin = 0;
% stackAccept = NaN(MCMCdraws,1);
while m < MCMCreps % using while, not for loop to allow going back in MCMC chain

    if m == 0
        % initializations
        PREVdraw.Y          = Y0;                                      % data matrix Y
        PREVdraw.X          = X0;                                      % data matrix X
        PREVdraw.A_         = eye(N);
        PREVdraw.PAI        = X0\Y0;
        PREVdraw.sqrtht     = sqrt([ARresid(1,:).^2; ARresid.^2]);     % Initialize sqrt_sqrtht
        PREVdraw.Vol_states = 2*log(PREVdraw.sqrtht);                  % Initialize states
        PREVdraw.sqrtPHI_   = sqrt(0.0001)*speye(N);                   % Initialize PHI_, a draw from the covariance matrix W

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


    %% STEP 2b: Draw from the conditional posterior of PAI
    PAI   = CTA(Y,X,N,K,T,A_,sqrtht,iV,iVb_prior,PAI,rndStream);
    RESID = Y - X*PAI; % compute the new residuals
    
    %% STEP 2c: Draw the covariances
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

    %% STEP 2d and STEP 1: Draw mixture states and then volatility states


    logy2 = log((RESID*A_').^2 + logy2offset);

    [Vol_states, ~, Vol_shocks] = StochVolKSCcorrsqrt(logy2', Vol_states', sqrtPHI_, Vol_0mean, Vol_0vcvsqrt, gridKSC, gridKSCt, N, T, rndStream);
    Vol_states = Vol_states';
    eta        = Vol_shocks';

    sqrtht     = exp(Vol_states/2); %compute sqrtht^0.5 from volatility states, needed in step 2b

    %% STEP 2a: Draw volatility variances
    Zdraw     = randn(rndStream, N, T + d_PHI);

    sqrtPHIpost = chol(s_PHI + eta'*eta, 'lower');
    sqrtZZ      = chol(Zdraw * Zdraw'); % note: right uppper choleski
    sqrtPHI_    = sqrtPHIpost / sqrtZZ; % just a square root, not choleski
    PHI_        = sqrtPHI_ * sqrtPHI_'; % derive posterior draw from PHI
    sqrtPHI_    = sparse(chol(PHI_, 'lower'));  % now choleski, and sparse (helps with performance of precision-based sampler)

    %% EM NEW STEP: shadow rate sampling
    if hasELBdata

        % always apply Geweke Gibbs

        % update state space objects
        elb.Y               = Y(elbT0+1:end,:)';

        % adjust obs for lagged actual rates
        PAIactual  = PAI(Kshadow+1:end,:);
        % Yhatactual = X(elbT0+1:end,Kshadow+1:end) * PAIactual;
        Yhatactual = transpose(Xffrlags(elbT0+1:end,:) * PAIactual);


        % update VAR companion form
        PAIshadow           = PAI(1:Kshadow,:);
        elb.A(elb.ndxY, :)  = PAIshadow';
        elb.B(elb.ndxY,:)   = invA_;


        % map SV into sqrtSigma
        elb.sqrtSigma = sqrtht(elbT0+1:end,:)';


        elbYgibbs           = elb.Y;
        elb.Y(elb.yNaN)     = NaN; % missing values

        %% PS setup
        pai0    = transpose(PAIshadow(1,:)) + Yhatactual;
        pai3    = reshape(transpose(PAIshadow(2:end,:)), N, N, p);
        invbbb  = A_ ./ permute(elb.sqrtSigma, [1 3 2]);
        
        if m == 1
            [~, CC, QQ, RR1, arows, acols, asortndx, brows, bcols, bsortndx] = VARTVPSVprecisionsamplerNaN(pai3,invbbb,elb.Y,elb.yNaN,elbY0,pai0,rndStream);
        end

        % b) reconstruct Y and X
        shadowYdata = Ydata;
        if doELBsampling

            shadowrate = gibbsdrawShadowrates(elbYgibbs, elb.X0, Yhatactual, elb.ndxS, elb.sNaN, p, elb.A, elb.B, elb.sqrtSigma, ...
                elb.bound, 1, elb.gibbsburn, rndStream);

            if doELBsampleAlternate
                YYdraws = VARTVPSVprecisionsamplerNaN(pai3,invbbb,elb.Y,elb.yNaN,...
                    elbY0,pai0,rndStream, ...
                    CC, QQ, RR1, arows, acols, asortndx, brows, bcols, bsortndx);
                YYdraws = reshape(YYdraws, N, elbT);
                missingrate = YYdraws(ndxSHADOWRATE,:);
            else
                missingrate = NaN;
            end

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
        Y = shadowYdata(p+1:end,:);

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
            PHI_all(thisMCMCdraw,:)          = PHI_((tril(PHI_))~=0);
            invA_all(thisMCMCdraw,:,:)       = invA_;
            sqrtht_all(thisMCMCdraw,:,:)     = sqrtht;
            if hasELBdata
                shadowrate_all(thisMCMCdraw,:,:)  = shadowrate;
                missingrate_all(thisMCMCdraw,:,:) = missingrate;
            end


            %% compute OOS draws
            

            if doPredictiveDensity

                % draw and scale SV shocks
                logSV0      = Vol_states(end,:)'; % Note: Vol_states record logs of *variances*
                logSVshocks = sqrtPHI_ * randn(rndStream, N, fcstNhorizons * Ndraws);
                logSVshocks = reshape(logSVshocks, N, fcstNhorizons, Ndraws);

                logSV       = logSV0 + cumsum(logSVshocks,2);
                fcstSVdraws = exp(logSV * 0.5);

                % draw random numbers and scale by SV
                nushocks    = fcstSVdraws .* randn(rndStream, N, fcstNhorizons, Ndraws);

                for nn = 1 : Ndraws

            
                    %% update VAR companion form

                    % note: tracks only actual rates with ELB data in estimation window
                    fcstA(ndxfcstY, :)         = PAI';

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
                            logdetOmegaY = sum(logSV(:,1,nn));
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
    %     PREVdraw.PHI_       = PHI_;
    PREVdraw.sqrtPHI_   = sqrtPHI_;

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

% if doELBsampling
%     fprintf('DONE with thisT %d, TID %d (countELBaccept = %d / %d, countELBacceptBurnin = %d / %d ) \n', thisT, TID, countELBaccept, MCMCdraws, countELBacceptBurnin, MCMCburnin)
% else
%     fprintf('DONE with thisT %d, TID %d \n', thisT, TID)
% end
fprintf('DONE with thisT %d, TID %d \n', thisT, TID)

return

