function [PAI,BETA]=CTAshadowratebeta(Y,X,Xstar,N,K,~,Nshadowrates,Nmacro,A_,sqrtht,iVpai,iVpaimean_prior,iVbeta,iVbetamean_prior,PAI,BETA,rndStream)

% adapated from CTA.m to handle general shadow-rate VAR

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
