function missingrateDraw = drawMissingrates(yData, yNaN, p, XX, XX0, A, B, sqrtSigma, rndStream)
% DRAWMISSINGRATES ...
%
%   ...

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 15-Feb-2020 13:53:49 $
% $Revision : 1.00 $
% DEVELOPED : 9.7.0.1296695 (R2019b) Update 4
% FILENAME  : drawMissingrates.m

%% Notation for setup
% y is VAR vector
% y = [x;s] where s are potentially missing variables
% Y, X, S are their companion vector analogues. Ny, Nx, and Nss are fixed
% observed data (signal) is Z = C Y

%% process arguments

[Ny,T] = size(yData);

% setup indices into y and Y
ndxS = any(yNaN,2);
Ns   = sum(ndxS);

% note: state vector contains unity as first element
ndxXX      = cat(1, true,  repmat(~ndxS, p, 1));
ndxSS      = cat(1, false, repmat(ndxS, p, 1));

Nstate = length(ndxXX);
if Nstate ~= Ny * p + 1
    error('dimension mismatch')
end

Nxx  = sum(ndxXX);
Nss  = sum(ndxSS);
if Nss ~= p * Ns
    error('dimension mismatch, Nss should be equal to p, but we have p=%d and Nss=%d', p, Nss)
end

% define inno vector *relative to reshuffled companion system*
Nx = Ny - Ns;
ndxINNOx = cat(1, false, true(Nx,1), false(Ns + Ny * (p-1),1));
% ndxINNOs = cat(1, false(Nxx,1), true(Ns,1), false(Ns * (p - 1), 1));
ndxSIGNAL  = ndxINNOx;
ndxSIGNAL(end-(Nss-1):end) = true;
% Nsignal  = sum(ndxSIGNAL);


Cy             = zeros(Ny,Nstate);
Cy(:,1+(1:Ny)) = eye(Ny);
Cy             = [Cy(:,ndxXX) Cy(:,ndxSS)]; % reorder states as in companion form

S0    = XX0(ndxSS);

X     = XX(ndxXX,:);
Xlag  = [XX0(ndxXX) XX(ndxXX,1:end-1)];

%% analyze time structure of missing data
anymiss = any(yNaN,1);
T0      = find(anymiss, 1, 'first'); % the first missing observation
Tstar   = find(anymiss, 1, 'last');  % the last missing observation

%% companion form VAR with partition into X and S

Axx = A(ndxXX,ndxXX);
Axs = A(ndxXX,ndxSS);
Asx = A(ndxSS,ndxXX);
Ass = A(ndxSS,ndxSS);
% A   = [Axx Axs; Asx Ass];
As  = [Axs; Ass];

Nw  = size(B,2);

Bx  = NaN(Nxx, Nw, T);
Bs  = NaN(Nss, Nw, T);

for t = 1 : T
    Bx(:,:,t) = B(ndxXX,:) * diag(sqrtSigma(:,t));
    Bs(:,:,t) = B(ndxSS,:) * diag(sqrtSigma(:,t));
end
BSV = cat(1, Bx, Bs);


%% Kalman filter
% shat stores prior/posterior mean
% Sigma stores posterior variance of S
% Omega stores prior variance of Y = C * X

Shat       = NaN(Nss, T);
Shattm1    = NaN(Nss,T);
Xhattm1    = NaN(Nxx,T);
sqrtSigmaS = NaN(Nss,Nss,T);


for t = 1 : T
    
    
    Zt = yData(~yNaN(:,t),t);
    C  = Cy(~yNaN(:,t),:);
    
    CC = C * As;
    DD = C * BSV(:,:,t);
    
    % compute Kalman gain
    if t == 1
        sqrtSigmaSlag = zeros(Nss);
        Shatlag       = S0;
    else
        sqrtSigmaSlag = sqrtSigmaS(:,:,t-1);
        Shatlag       = Shat(:,t-1);
    end
    
    % QR Kalman
    [~, K, sqrtSigmaS(:,:,t)] = kalmanQRupdateABCD(Ass, Bs(:,:,t), CC, DD, sqrtSigmaSlag);
    
    % update shat
    Shattm1(:,t) = Asx * Xlag(:,t) + Ass * Shatlag;
    Xhattm1(:,t) = Axx * Xlag(:,t) + Axs * Shatlag;

    STATEttm1    = [Xhattm1(:,t); Shattm1(:,t)];
    ztilde       = Zt - C * STATEttm1;
    Shat(:,t)    = Shattm1(:,t) + K * ztilde; % = Asx * Xlag(:,t) + Ass * Shatlag + K * ztilde;
    
end


%% Smoothing

missingrateDraw = NaN(Ns,T);

% values from Tstar + 1 : T are known, need to start drawing only as of
% Tstar + p

if T < (Tstar + p)

    % draw from end of filter posterior
    

    % catch cases where Tstar < T < (Tstar + p)

    Sdraw   = NaN(Nss,1);
    these   = T-(0:p-1);
    thisNdx = these <= Tstar;
    
    thisShat         = Shat(:,T);
    Sdraw(~thisNdx)  = thisShat(~thisNdx);
    
    % checkdiff(sqrtSigmaS(~thisNdx,:,T)); % check that values are really known
    
    % compute choleski; note: fails when T < p
    %     thisSqrtSigmaS   = chol(sqrtSigmaS(thisNdx,:,T) * sqrtSigmaS(thisNdx,:,T)')';
    %     checkdiff(thisSqrtSigmaS, cholqr(sqrtSigmaS(thisNdx,:,T)));
    thisSqrtSigmaS = cholqr(sqrtSigmaS(thisNdx,:,T));
    
    Sdraw(thisNdx) = thisShat(thisNdx) + thisSqrtSigmaS * randn(rndStream, sum(thisNdx), 1);
    
    % fill Sdraw into missingrateDraw
    % If T < p, Sdraw will contain lagged (actual) rates that need not be
    % filled into missingrateDraw; these(these > 0) tracks the relevant entries
    % into missingrateDraw; note further that these counts backwards
    %     these                              = T-(0:p-1);
    missingrateDraw(:,these(these > 0)) = Sdraw(1 : sum(these > 0));

else % T >= Tstar + p
    % fill in known values before handing over to the smoother
    for t = T : -1 : (Tstar + 1)
        missingrateDraw(:,t) = yData(ndxS,t); 
    end
    % note: values are filled in through Tstar+1, 
    % ... but jump off for smoother is Tstar + p - 1
    % ... this ensures that STATEdraw (created below) has no NaN entries
end

thisT = min(T-1, Tstar + p - 1);

drawNdx = thisT : -1 : T0 + (p - 1);
zdraws = randn(rndStream, Ns, length(drawNdx));

for t = drawNdx

    % t indexes date for which statevector is drawn, 
    % ... missingrate will be drawn for t-(p-1)
    
    % use QR
    AA            = As(ndxSIGNAL,:);
    BB            = BSV(ndxSIGNAL,:,t+1); 
    thisSqrtSigma = sqrtSigmaS(:,:,t);
    [J, sqrtSigmaMissingrate] = statesmootherQRupdate(AA, BB, thisSqrtSigma, Ns, t > Tstar);
    
    
    % construct SdrawHat 
    Sdraw          = missingrateDraw(:,t+1-(0:p-1))';
    STATEdraw      = cat(1, X(:,t+1), Sdraw);
    STATEhat       = cat(1, Xhattm1(:,t+1), Shattm1(:,t+1));
    SIGNAL         = STATEdraw(ndxSIGNAL) - STATEhat(ndxSIGNAL);
    missingrateHat  = Shat(end-(Ns-1):end,t) + J * SIGNAL;

    missingrateDraw(:,t-(p-1)) = missingrateHat + sqrtSigmaMissingrate * zdraws(:,drawNdx == t);

end

if any(isnan(missingrateDraw(:)))
    error('NaN missingrate draws')
end

for t = T0 - 1 : -1 : 1 % T0 should be equal to one anyway
    missingrateDraw(:,t) = yData(ndxS,t); % just fill in data
end

function [J, sqrtSigmaPosterior,sqrtOmega] = statesmootherQRupdate(A, B, sqrtSigmaPrior, Ns, doPinv)
% QR decomposition for shadow-rate state smoother, with extended signal vector
%

if nargin < 4 || isempty(Ns)
    Ns = 1;
end

if nargin < 5
    doPinv = true;
end

if Ns ~= 1
    error('code not vetted for Ns > 1')
end

[Nz,Nw]    = size(B);
[~,Nstate] = size(A);



if Nstate + Nw < Nz + Ns
    warning('Should have Nstate + Nw > Nz + Ns but get %d+%d < %d+Ns', Nstate, Nw, Nz + Ns) 
end


M = zeros(Nz + Ns, Nw + Nstate);

M(1:Nz, 1:Nw)                 = B;

M(1:Nz, Nw + (1:Nstate))      = A * sqrtSigmaPrior;

H                             = [zeros(Ns,Nstate-Ns), eye(Ns)];
M(Nz+(1:Ns), Nw+(1:Nstate))   = H * sqrtSigmaPrior;

[~,R] = qr(M');

R = R';

sqrtOmega          = R(1:Nz,1:Nz);

if doPinv
    J  = R(Nz+(1:Ns),1:Nz) * pinv(sqrtOmega);
else
    if rcond(sqrtOmega) > 1e-10
        J  = R(Nz+(1:Ns),1:Nz) / sqrtOmega;
    else
        % warning('ill conditioned sqrtOmega, using pinv')
        J  = R(Nz+(1:Ns),1:Nz) * pinv(sqrtOmega);
    end
end

sqrtSigmaPosterior = R(Nz+(1:Ns), Nz+(1:Ns));




