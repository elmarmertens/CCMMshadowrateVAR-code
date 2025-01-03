function [PAI,BETA]=CTAshadowratebeta(Y,X,Xstar,N,K,~,Nshadowrates,Nmacro,A_,sqrtht,iVpai,iVpaimean_prior,iVbeta,iVbetamean_prior,PAI,BETA,rndStream)

% =========================================================================
% Performs a draw from the conditional posterior of the VAR conditional
% mean coefficients by using the triangular algorithm. The
% triangularization achieves computation gains of order N^2 where N is the
% number of variables in the VAR. Carriero, Clark and Marcellino (2015),
% Large Vector Autoregressions with stochastic volatility and flexible
% priors. 
%
% The model is:
%
%     Y(t) = Pai(L)Y(t-1) + v(t); Y(t) is Nx1; t=1,...,T, L=1,...,p.
%     v(t) = inv(A)*(LAMBDA(t)^0.5)*e(t); e(t) ~ N(0,I);
%                _                                         _
%               |    1          0       0       ...      0  |
%               |  A(2,1)       1       0       ...      0  |
%      A =      |  A(3,1)     A(3,2)    1       ...      0  |
%               |   ...        ...     ...      ...     ... |
%               |_ A(N,1)      ...     ...   A(N,N-1)    1 _|
%
%     Lambda(t)^0.5 = diag[sqrt_h(1,t)  , .... , sqrt_h(N,t)];
%
% INPUTS
% Data and pointers:
% Y     = (TxN) matrix of data appearing on the LHS of the VAR
% X     = (TxK) matrix of data appearing on the RHS of the VAR
% N     = scalar, #of variables in VAR 
% K     = scalar, #of regressors (=N*p+1)  
% T     = scalar, #of observations
% The matrix X needs to be ordered as: [1, y(t-1), y(t-2),..., y(t-p)]
% 
% Error variance stuff:
% invA_   = (NxN) inverse of lower triangular covariance matrix A
% sqrtht = (TxN) time series of diagonal elements of volatility matrix 
% For a homosckedastic system, with Sigma the error variance, one can
% perform the LDL decomposition (command [L,D]=LDL(Sigma)) and set inv_A=L
% and sqrtht=repmat(sqrt(diag(D)'),T,1). 
%
% Priors:
% iV          = (NKxNK) precision matrix for VAR coefficients 
% iVB_prior   = (NKx1) (prior precision)*(prior mean)
% Note 1:iV is the inverse of the prior matrix and iVB_prior is the product
% of iV and the prior mean vector, which both need to be computed only once,
% before the start of the main MCMC loop.
% Note 2:in this code, iV is assumed block-diagonal. This implies that the
% prior is independent across equations. This includes most of the priors 
% usually considered, including the Minnesota one.  To use a non-block
% diagonal iV one needs to modify the code using the recursions illustrated 
% in equations (37) and (38).  
%
% OUTPUT
% One draw from (PAI|A,Lambda,data)
% PAI=[Pai(0), Pai(1), ..., Pai(p)].
% =========================================================================

% y_til=Y*A_';

Kmax                    = K + Nshadowrates; % not all equations contain beta, but let's draw rnds for all
zdraws                  = randn(rndStream,Kmax,N);
BBB                     = zeros(Nshadowrates,N);
BBB(:,end-Nmacro+1:end) = BETA;
otherX                  = [X Xstar];

for j=1:N


    
    if j > N - Nmacro
        indexPAI             = K*(j-1)+(1:K);
        indexBETA            = Nshadowrates*(j-1-(N-Nmacro))+(1:Nshadowrates);
        thisPriorVi          = blkdiag(iVpai(indexPAI,indexPAI), iVbeta(indexBETA,indexBETA));
        thisPriorViMean      = cat(1, iVpaimean_prior(indexPAI), iVbetamean_prior(indexBETA));
        thisX                = otherX;
        thisK                = Kmax;
    else
        indexPAI        = K*(j-1)+1:(K*(j-1)+K);
        thisPriorVi     = iVpai(indexPAI,indexPAI);
        thisPriorViMean = iVpaimean_prior(indexPAI);
        thisX           = X;
        thisK           = K;
    end

    otherSLOPES          = cat(1, PAI, BBB);

    otherSLOPES(:,j) = 0; % to ignore slopes associated with current equation

    % build model
    lambda = vec(sqrtht(:,j:N));
    Y_j    = vec((Y - otherX * otherSLOPES) * A_(j:N,:)')./lambda; 
    X_j    = kron(A_(j:N,j),thisX)./lambda;
    
    % posterior moments
    iV_post = thisPriorVi  + X_j'*X_j;
    [iVchol_post, ispd]  = chol(iV_post, 'lower');
    if ispd == 0

        Vchol_post      = (iVchol_post \ eye(thisK))';  % note: Vchol_post used below for constructing draws; hence the transpose
        V_post          = Vchol_post * Vchol_post'; 
    
    else % proceed via QR
        
        warning('switching to QR routine')
        
        iVchol = chol(thisPriorVi)'; % could do this once and for all, not for every call of triang
        % set up Kailath's fast array matrix
        qrM = [iVchol, X_j'];
        R = triu(qr(qrM',0))';
    
        iVchol_post = R(1:thisK,1:thisK);
        Vchol_post  = (iVchol_post \ eye(thisK))'; % upper triangular but OK
        V_post      = Vchol_post * Vchol_post';
    end
    
    % posterior draw
    b_post      = V_post*(thisPriorViMean + X_j'*Y_j);
    theseSLOPES = b_post + Vchol_post * zdraws(1:thisK,j);
    PAI(:,j)    = theseSLOPES(1:K);
    if j > N - Nmacro
        BBB(:,j) = theseSLOPES(K+1:end);
    end
end
BETA = BBB(:,end-Nmacro+1:end);                    
