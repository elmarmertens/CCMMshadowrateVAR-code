function [PAI_all, PHI_all, invA_all, sqrtht_all, shadowrate_all, missingrate_all, ...
    fcstYdraws, fcstYhat, ...
    fcstShadowrateDraws, fcstShadowrateHat ] = ...
    mcmcVARhybridshadowrate(thisT, MCMCdraws, ...
    p, np, data0, ydates0, ...
    actualrateBlock, ...
    minnesotaPriorMean, doTightPrior, ...
    ndxSHADOWRATE, ndxOTHERYIELDS, doELBsampling, doSampleBoth, ELBbound, elbT0, check_stationarity, ...
    fcstNdraws, fcstNhorizons, rndStream, doprogress)


if nargin < 20
    doprogress = false;
end

%% get TID
% used to provide context for warning messages
TID   = parid;

doPredictiveDensity = nargout > 6;

%% truncate sample
samEnd = ydates0(thisT);
ndx    = ydates0 <= samEnd;

data   = data0(ndx,:);

%% --------------------------- OPTIONS ----------------------------------
if doTightPrior
    theta=[0.05 0.5 100 2];       % hyperparameters of Minnesota prior:
else
    theta=[0.1 0.5 100 2];       % hyperparameters of Minnesota prior:
end
% [lambda1 lambda2 int lambda3], int is the
% prior on the intercept. lambda1, lambda2
% and lambda3 are as in equation (42) with
% lambda1 the overall shrinkage, lambda2 the
% cross shrinkage and lambda 3 the lag decay
% (quadratic if =2). Note lambda2~=1 implies
% the prior becomes asymmetric across eqation,
% so this would not be implementable in the
% standard conjugate setup.
% Minn_pmean = 0;              % Prior mean of the 1-st own lag for each
% equation. For nonstationary variables, this
% is usually set to 1. For transformed
% stationary variables this is set to 0.

burnin            = ceil(0.5 * MCMCdraws);    % burn in
MCMCreps          = MCMCdraws + burnin;        % total MCMC draws


%% -------------------------Create data matrices-------------------------
% pointers
[Nobs,N]                 = size(data);
Nshadowrates             = length(ndxSHADOWRATE);
Ydata                    = data;
shadowrateData           = data(:,ndxSHADOWRATE);
ndxELB                   = shadowrateData <= ELBbound;
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
Nstates       = K + Nshadowrates * p; % track lagged actual rates as well
Xjumpoff      = zeros(K,1);
Xjumpoff(1)   = 1;
for l=1:p
    Xjumpoff(1+(l-1)*N+(1:N)) = data(Nobs-(l-1),1:N);
end
% add lagged actual rates to jumpoff
for l=1:p
    Xjumpoff(K+(l-1)*Nshadowrates+(1:Nshadowrates)) = data(Nobs-(l-1),ndxSHADOWRATE);
end

ndxSHADOWRATELAGS = cat(2, false, repmat(ismember(1:N, ndxSHADOWRATE), 1, p)); % preprend by false for CONST

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
fcstYdraws = NaN(N,fcstNhorizons, Ndraws, MCMCdraws); % see reshape near end of script

%% prepare state space for forecasting
% note: there are some redundacies between the forecast state space and the elb state space
% the redundancies allow, however, for easier comparison and maintenance of codes for different model version
fcstA                  = zeros(Nstates,Nstates);
fcstA(1,1)             = 1; % unit root for the constant
fcstA(1+N+1:K,2:K)     = [eye(N*(p-1)),zeros(N*(p-1),N)]; % fill in lower part of companion form
fcstA(K+Nshadowrates+(1:Nshadowrates*(p-1)),K+(1:Nshadowrates*(p-1))) = eye(Nshadowrates*(p-1)); % companion for actual rates
ndxfcstActual          = K+(1:Nshadowrates);
ndxfcstShadow          = 1+ndxSHADOWRATE;

ndxfcstY          = 1+(1:N);
fcstB             = zeros(Nstates,N);
fcstB(ndxfcstY,:) = eye(N);


%% EM, new: prepare ELB data
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

if hasELBdata
    
    elb.gibbsburn      = 1e2;
    elb.ndxS           = ismember(1:N, ndxSHADOWRATE);
    
    elb.ndxY           = 1+(1:N);
    elb.T0             = elbT0;
    elb.T              = T - elb.T0;
    elb.yNaN           = yNaN(elb.T0+1:end,:)'; % note the transpose, needed for call of sampler later
    elb.sNaN           = elb.yNaN(elb.ndxS,:);
        
    if any(any(yNaN(1:elb.T0,ndxSHADOWRATE)))
        error('something off about elbT0')
    end
    
    if ~any(any(yNaN(elb.T0+1,ndxSHADOWRATE)))
        warning('elb.T0 + 1 should be a missing obs, but it is not ...')
    end
    
    elb.bound = ELBbound;
    
    
    % init transition matrix
    elb.A                = zeros(K,K);
    elb.A(1,1)           = 1; % unit root for the constant
    elb.A(1+N+1:end,2:end) = [eye(N*(p-1)),zeros(N*(p-1),N)]; % fill in lower part of companion form
    
    elb.B            = zeros(K,N);
    
    % encode "prior" over initial conditions (which are fixed)
    elb.X0  = X(elb.T0+1,:)'; % time zero values of state. Recall that X already contains lagged values
    
    % construct state vector for ELB state space
    dummy           = Ydata;
    dummy(YdataNaN) = NaN;
    dummy(:,ndxSHADOWRATE) = NaN; % to ignore all information about a actual or shadow rate
    elb.X   = ones(1+N*p,elb.T);
    for l=0:p-1
        elb.X(1+N*l+(1:N),:) = dummy(p+elb.T0-l+(1:elb.T),1:N)'; % note: T0 indexes into data after cutting out p lags
    end
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
% todo: Pi_pv could be encoded as vector, only using diagonal elements anyhow
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

% Pai~N(vec(MU_pai),OMEGA_pai), equation 7.
OMEGA_pai   = diag(vec([sigma_const;reshape(diag(Pi_pv),Klagreg,N)])); % prior variance of Pai
MU_pai      = [zeros(1,N);reshape(Pi_pm,Klagreg,N)];                   % prior mean of Pai

% A~N(MU_A,inv(OMEGA_A_inv)), equation 8.
MU_A = NaN(N-1,N);
OMEGA_A_inv = NaN(N-1,N-1,N);
for i = 2:N;
    MU_A(1:i-1,i) = zeros(i-1,1);             % prior mean of A
    OMEGA_A_inv(1:i-1,1:i-1,i) = 0*eye(i-1);  % prior precision of A
end;

% PHI~IW(s_PHI,d_PHI), equation 9.
d_PHI = N+3;                 % prior dofs
s_PHI = d_PHI*(0.15*eye(N)) * 12 / np; % prior scale, where eye(N)=PHI_ in equation (9)

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

if hasELBdata
    shadowrate_all  = NaN(MCMCdraws,Nshadowrates,elb.T);
    missingrate_all = NaN(MCMCdraws,Nshadowrates,elb.T);
    shadowrate      = NaN(Nshadowrates,elb.T);
    missingrate     = NaN(Nshadowrates,elb.T);
else
    shadowrate_all  = [];
    missingrate_all = [];
    shadowrate      = [];
    missingrate     = [];
end

% define some useful matrices prior to the MCMC loop
% PAI  = zeros(K,N);                                          % pre-allocate space for PAI
comp = [eye(N*(p-1)),zeros(N*(p-1),N)];                       % companion form
iV   = diag(1./diag(OMEGA_pai)); iVb_prior=iV*vec(MU_pai);    % inverses of prior matrices
EYEn = eye(N);


%% start of MCMC loop
% the algorithm is as described in page 15, but it starts from step 2b,
% this is the same as starting from step 1 but is more convenient as
% it requires less initializations (one can think of steps 2b to 2d in
% repetition 1 as an initialization).

if doprogress
    progressbar(0);
end
m = 0;

while m < MCMCreps % using while, not for loop to allow going back in MCMC chain
    
    if m == 0
        % initializations
        PREVdraw.Y          = Y0;                                      % data matrix Y
        PREVdraw.X          = X0;                                      % data matrix X
        A_                  = eye(N);                                  % initialize A matrix
        PREVdraw.A_         = A_;
        PREVdraw.PAI        = X0\Y0;
        PREVdraw.sqrtht    = sqrt([ARresid(1,:).^2; ARresid.^2]);     % Initialize sqrt_sqrtht
        PREVdraw.Vol_states = 2*log(PREVdraw.sqrtht);                 % Initialize states
        %         PREVdraw.PHI_       = 0.0001*eye(N);                           % Initialize PHI_, a draw from the covariance matrix W
        PREVdraw.sqrtPHI_   = sqrt(0.0001)*speye(N);                           % Initialize PHI_, a draw from the covariance matrix W
        
    end % m == 0
    m = m + 1;
    
    
    % init with previous draws values
    Y          = PREVdraw.Y;
    X          = PREVdraw.X;
    A_         = PREVdraw.A_;
    sqrtht    = PREVdraw.sqrtht;
    Vol_states = PREVdraw.Vol_states;
    PAI        = PREVdraw.PAI;
    sqrtPHI_   = PREVdraw.sqrtPHI_;
    
    % if mod(m,10) == 0; clc; disp(['percentage completed:' num2str(100*m/MCMCreps) '%']); toc; end
    
    %% STEP 2b: Draw from the conditional posterior of PAI, equation 10.
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
    
    %% STEP 2c: Draw the covariances, equation 11.
    for ii = 2:N
        % weighted regression to get Z'Z and Z'z (in Cogley-Sargent 2005 notation)
        y_spread_adj=RESID(:,ii)./sqrtht(:,ii);
        %         X_spread_adj=[]; for vv=1:ii-1;  X_spread_adj=[X_spread_adj RESID(:,vv)./sqrtht(:,ii)]; end  %#ok<AGROW>
        X_spread_adj = RESID(:,1 : ii - 1) ./ sqrtht(:,ii); % note: use of implicit vector expansion
        ZZ=X_spread_adj'*X_spread_adj; Zz=X_spread_adj'*y_spread_adj;
        % computing posteriors moments
        Valpha_post = (ZZ + OMEGA_A_inv(1:ii-1,1:ii-1,ii))\eye(ii-1);
        alpha_post  = Valpha_post*(Zz + OMEGA_A_inv(1:ii-1,1:ii-1,ii)*MU_A(1:ii-1,ii));
        % draw and store
        alphadraw   = alpha_post+chol(Valpha_post,'lower')*randn(rndStream,ii-1,1);
        A_(ii,1:ii-1)= -1*alphadraw';
        % EM: Note: need to init A_ at least once to set A_(:,end)
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
    PHI_        = sqrtPHI_ * sqrtPHI_'; % derive posterior draw from PHI, equation 12.
    sqrtPHI_    = sparse(chol(PHI_, 'lower'));  % now choleski, and sparse (helps with performance of precision-based sampler)
    
    %% EM NEW STEP: shadow rate sampling
    if hasELBdata
        
        % always apply Geweke Gibbs
        
        % update state space objects
        elb.Y               = Y(elb.T0+1:end,:)';
        
        % adjust obs for lagged actual rates
        PAIactual                     = PAI(ndxSHADOWRATELAGS,:);
        PAIactual(:,~actualrateBlock) = 0;
        Yhatactual = Xactual(elb.T0+1:end,ndxSHADOWRATELAGS) * PAIactual;
        elb.Y      = elb.Y - Yhatactual';
        %         checkdiff(Yhatactual(:,ndxSHADOWRATE));
        
        % update VAR companion form
        PAIshadow                                    = PAI;
        PAIshadow(ndxSHADOWRATELAGS,actualrateBlock) = 0;
        elb.A(elb.ndxY, :)                           = PAIshadow';
        
        
        elb.B(elb.ndxY,:)   = invA_;
        
        % map SV into sqrtSigma
        elb.sqrtSigma = sqrtht(elb.T0+1:end,:)';

        if doSampleBoth
            shadowrate = gibbsdrawShadowrates(elb.Y, elb.X0, elb.ndxS, elb.sNaN, p, elb.A, elb.B, elb.sqrtSigma, ...
                elb.bound, 1, elb.gibbsburn, rndStream);
            missingrate = gibbsdrawShadowrates(elb.Y, elb.X0, elb.ndxS, elb.sNaN, p, elb.A, elb.B, elb.sqrtSigma, ...
                [], 1, elb.gibbsburn, rndStream);
        else
            if doELBsampling
                shadowrate = gibbsdrawShadowrates(elb.Y, elb.X0, elb.ndxS, elb.sNaN, p, elb.A, elb.B, elb.sqrtSigma, ...
                    elb.bound, 1, elb.gibbsburn, rndStream);
            else
                missingrate = gibbsdrawShadowrates(elb.Y, elb.X0, elb.ndxS, elb.sNaN, p, elb.A, elb.B, elb.sqrtSigma, ...
                    [], 1, elb.gibbsburn, rndStream);
            end
        end
        
        % b) reconstruct Y and X
        shadowYdata = Ydata;
        if doELBsampling
            shadowYdata(p+elb.T0+1:end,ndxSHADOWRATE) = shadowrate';
        else
            shadowYdata(p+elb.T0+1:end,ndxSHADOWRATE) = missingrate';
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
            Xjumpoff(K+(l-1)*Nshadowrates+(1:Nshadowrates)) = data(Nobs-(l-1),ndxSHADOWRATE);
        end
        
        
    end % hasELBdata
    
    %% post burnin: store draws and draw from oos-predictive density
    if m > burnin;
        
        thisdraw = m-burnin;
        % STORE DRAWS
        PAI_all(thisdraw,:,:)        = PAI;
        PHI_all(thisdraw,:)          = PHI_((tril(PHI_))~=0);
        invA_all(thisdraw,:,:)       = invA_;
        sqrtht_all(thisdraw,:,:)     = sqrtht;
        if hasELBdata
            shadowrate_all(thisdraw,:,:)  = shadowrate;
            missingrate_all(thisdraw,:,:) = missingrate;
        end
        
        
        %% compute OOS draws
        %
        
        if doPredictiveDensity
		
            % draw and scale SV shocks
            logSV0      = Vol_states(end,:)'; % Note: Vol_states record logs of *variances*
            logSVshocks = sqrtPHI_ * randn(rndStream, N, fcstNhorizons * Ndraws);
            logSVshocks = reshape(logSVshocks, N, fcstNhorizons, Ndraws);
            
            % draw random numbers
            zdraws  = randn(rndStream, N, fcstNhorizons, Ndraws);
            
            for nn = 1 : Ndraws
                
                logSV       = bsxfun(@plus, logSV0, cumsum(logSVshocks(:,:,nn),2));
                
                
                fcstSVdraws         = exp(logSV * 0.5);
                nushocks            = zeros(N, fcstNhorizons+1); % padding with a line of zeros for use with ltitr
                nushocks(:,1:end-1) = invA_ * (fcstSVdraws .* zdraws(:,:,nn));
                
                %% update VAR companion form
                PAIactual                     = PAI(ndxSHADOWRATELAGS,:);
                PAIactual(:,~actualrateBlock) = 0;
                
                PAIshadow = PAI;
                PAIshadow(ndxSHADOWRATELAGS,actualrateBlock) = 0;
                
                fcstA(ndxfcstY, 1:K)         = PAIshadow';
                fcstA(ndxfcstY, K+1:Nstates) = PAIactual';
                
                %% forecast simulation
                fcstX0 = Xjumpoff;
                
                for n = 1 : fcstNhorizons
                    fcstXdraw    = fcstA * fcstX0 + fcstB * nushocks(:,n);
                    ydraw        = fcstXdraw(ndxfcstY);
                    
                    % collect draw
                    fcstYdraws(:,n,nn,thisdraw) = ydraw;
                    % prepare next iteration
                    fcstX0                   = fcstXdraw;
                    fcstX0(ndxfcstActual)    = max(fcstX0(ndxfcstShadow), ELBbound);
                    
                end
            end % nn
            
        end
    end
    
    
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

% censor shadowrate draws as needed
ndxYIELDS = union(ndxSHADOWRATE, ndxOTHERYIELDS);

% uncensored forecasts
fcstShadowrateDraws = fcstYdraws(ndxYIELDS,:,:);
fcstShadowrateHat   = mean(fcstShadowrateDraws,3);

% collect draws
fcstYdraws = reshape(fcstYdraws, N, fcstNhorizons, fcstNdraws);

% censor shadow rates
yieldDraws                 = fcstYdraws(ndxYIELDS,:,:);
ndx                        = yieldDraws < ELBbound;
yieldDraws(ndx)            = ELBbound;
fcstYdraws(ndxYIELDS,:,:)  = yieldDraws;
% fcstYhat(ndxYIELDS,:)      = mean(fcstYdraws(ndxYIELDS,:,:),3);

fcstYhat   = mean(fcstYdraws,3);



fprintf('DONE with thisT %d, TID %d \n', thisT, TID)

return

