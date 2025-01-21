function PAI=CTAsys(Y,X,N,K,T,A_,sqrtht,iV,iVb_prior,PAI,rndStream)

% adapted from CTA.m to handle regressor switches in restricted non-structural shadow-rate VAR

zdraws = randn(rndStream,K,N);
XPAI   = NaN(T,N);

Ik       = eye(K);
for j=1:N
    
    % select coefficients of equation j to remove from the LHS
    %     PAI(:,j)=zeros(K,1);

    for jj = 1 : N
        if jj == j
            XPAI(:,jj) = 0;
        else
            XPAI(:,jj) = X(:,:,jj) * PAI(:,jj);
        end
    end
    
    % build model
    lambda=vec(sqrtht(:,j:N));
    Y_j=vec((Y-XPAI)*A_(j:N,:)')./lambda; 
    
    X_j=kron(A_(j:N,j),X(:,:,j))./lambda;
    
    % posterior moments
    index=K*(j-1)+1:(K*(j-1)+K);
    iV_post = iV(index,index)  + X_j'*X_j;
    [iVchol_post, ispd]  = chol(iV_post, 'lower');
    if ispd == 0

        Vchol_post      = (iVchol_post \ Ik)';  % note: Vchol_post used below for constructing draws; hence the transpose
        V_post          = Vchol_post * Vchol_post'; 
    
    else % proceed via QR
        
        warning('switching to QR routine')
        
        iVchol = chol(iV(index,index))'; % could do this once and for all, not for every call of triang
        % set up Kailath's fast array matrix
        qrM = [iVchol, X_j'];
        R = triu(qr(qrM',0))';
    
        iVchol_post = R(1:K,1:K);
        Vchol_post  = (iVchol_post \ Ik)'; % upper triangular but OK
        V_post      = Vchol_post * Vchol_post';
    end
    
    % posterior draw
    b_post   = V_post*(iVb_prior(index) + X_j'*Y_j);
    PAI(:,j) = b_post + Vchol_post * zdraws(:,j);
    
end
                    
  
