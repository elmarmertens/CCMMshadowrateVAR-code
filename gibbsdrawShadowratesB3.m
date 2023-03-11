function shadowrateDraws = gibbsdrawShadowratesB3(Y, STATE0, ndxS, sNaN, p, A, B, SVol, elbBound, Ndraws, burnin, rndStream)
% GIBBSDRAWSHADOWRATESB3 ...
%
%   ...


% B3 allows for 3D matrix B

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 05-Mar-2021 14:01:43 $
% $Revision : 1.00 $
% DEVELOPED : 9.9.0.1592791 (R2020b) Update 5
% FILENAME  : gibbsdrawShadowrates.m



%% process arguments


if nargin < 10 || isempty(Ndraws)
    Ndraws = 1;
end
if nargin < 11 || isempty(burnin)
    burnin = 0;
end
if nargin < 12
    rndStream = getDefaultStream;
end

[Ny,T] = size(Y);


S      = Y(ndxS,:);
Ns     = sum(ndxS);

shadowrateDraws = NaN(Ns,T,Ndraws);

%% prepare state space

Nstate  = Ny * p;     % without intercept
NNstate = Nstate + 1; % with    intercept

Nx            = Ny - Ns;
ndxX          = ~ndxS;
C             = zeros(Ns + Nx,NNstate);
C(:,1+(1:Ny)) = eye(Ny);

if ismatrix(B) % 2D
    B = repmat(B, [1 1 T]);
end

Nw  = size(B,2);
if Nw ~= Ny
    error('dimension mismatch');
end

Bsvol = zeros(Ny, Nw, T);
bb    = B(1+(1:Ny),:,:); % w/o constant
for t = 1 : T
    Bsvol(:,:,t) = bb(:,:,t) * diag(SVol(:,t));
end


%% prepare smoothing weights
J                  = NaN(Ns, Nstate + Nx, T);
sqrtOmegaPosterior = NaN(Ns,Ns,T);

aa             = A(2:end,2:end); % drop constant
Apowerp        = NaN(Nstate, Nstate, p+1);
Apowerp(:,:,1) = eye(Nstate);
for k = 1 : p
    Apowerp(:,:,k+1) = aa * Apowerp(:,:,k);
end
Apowerp = Apowerp(:,1:Ny,:);


% smoothing weights
for t = 1 : T-p
    
    if any(sNaN(:,t))
        % prepare M
        M                                    = zeros(Nstate + Nw, Nw * (p + 1));
        for j = 0 : p
            M(1 : Nstate, Nw * j + (1 : Nw)) = Apowerp(:,:,p+1-j) * Bsvol(:,:,t+j);
        end
        M(Nstate + (1 : Nx), 1 : Nw)         = Bsvol(ndxX,:,t);
        M(Nstate + Nx + 1 : end, 1 : Nw)     = Bsvol(ndxS,:,t);
        
        
        % perform qr
        [~,R] = qr(M');
        R     = R';
        
        sqrtSigma                  = R(1:Nstate+Nx,1:Nstate+Nx);
        J(:,:,t)                   = R(Nstate+Nx+(1:Ns),1:Nstate+Nx) / sqrtSigma;
        sqrtOmegaPosterior(:,:,t)  = R(Nstate+Nx+(1:Ns), Nstate+Nx+(1:Ns));
    end
    
end

if isempty(t) % if T-p<1
    t = 0;
end

while t < T
    
    t = t + 1;
    
    if any(sNaN(:,t))
        
        k = T - t;
        Nsignal = Ny * k + Nx;
        
        thisM                               = zeros(Nsignal + Ns, Nsignal + Ns);
        for j = 0 : k
            thisM(1 : Nsignal, Nw * j + (1 : Nw)) = Apowerp(1:Nsignal,:,k+1-j) * Bsvol(:,:,t+j);
        end
        thisM(k * Ny + (1 : Nx), 1 : Nw)    = Bsvol(ndxX,:,t);
        thisM(Nsignal + 1 : end, 1 : Nw)    = Bsvol(ndxS,:,t);
        
        [~,R] = qr(thisM');
        R     = R';
        
        sqrtSigma                      = R(1 : Nsignal, 1 : Nsignal);
        J(:,:,t)                       = 0;
        J(:, Ny * (p - k) + 1 : end,t) = R(Nsignal + (1:Ns), 1 : Ny * k + Nx) / sqrtSigma; % "live values" at bottom of state vector
        sqrtOmegaPosterior(:,:,t)      = R(Nsignal + (1:Ns), Nsignal + (1:Ns));
        
        
    end
end

%% weights for truncMV in case of Ns > 1
if Ns > 1
    sqrtOmega1 = NaN(Ns, T);
    beta1      = NaN(Ns, Ns-1, T);
    for t = 1 : T
        if any(sNaN(:,t))
            vcv = sqrtOmegaPosterior(:,:,t) * sqrtOmegaPosterior(:,:,t)';
            for s = 1 : Ns
                ndxOther = (1:Ns ~= s);
                beta1(s,:,t)    = vcv(s,ndxOther) / vcv(ndxOther,ndxOther);
                sqrtOmega1(s,t) = sqrt(vcv(s,s) - vcv(s,ndxOther) / vcv(ndxOther,ndxOther) * vcv(ndxOther,s));
            end
        end
    end
    
    % checkdiff(sqrtOmega1(end,:), abs(sqrtOmegaPosterior(end,end,:)));
end

%% prepare further state-space objects
CAA   = C * A;
AApp1 = A^(p+1); % including constant
AApp1 = AApp1(2:end,:); % dropping constant

%% Gibbs draws
totalNdraws = burnin + Ndraws;

if isempty(elbBound)
    zdraws = randn(rndStream, Ns, T, totalNdraws);
else
    udraws = rand(rndStream, Ns, T, totalNdraws);
end


% progressbar(0)
for n = 1 : totalNdraws
    
    % init t = 1
    STATElag        = STATE0;
    YY              = cat(2, Y, zeros(Ny, p)); % padded (dummy) values to compute STATEfuture below
    
    for t = 1 : T
        
        if any(sNaN(:,t))
            
            Yhat           = CAA * STATElag;
            Xresid         = Y(ndxX,t) - Yhat(ndxX);
            Shat           = Yhat(ndxS);
            
            STATEfuturehat = AApp1 * STATElag;
            STATEfuture    = YY(:,t+p:-1:t+1); % using zeros as dummy values for t + k > T
            STATEtilde     = STATEfuture(:) - STATEfuturehat;
            Sposterior     = Shat + J(:,:,t) * [STATEtilde; Xresid]; % note: could adapt to t > T - p
            
            if isempty(elbBound)
                S(:,t) = Sposterior + sqrtOmegaPosterior(:,:,t) * zdraws(:,t,n);
            else
                if Ns == 1
                    
                    S(:,t) = drawTruncNormal(Sposterior, sqrtOmegaPosterior(:,:,t), elbBound, udraws(:,t,n));
                    
                else
                    
                    for s = find(sNaN(:,t)')
                        ndxOther = 1 : Ns ~= s;
                        thisMu  = Sposterior(s) + beta1(s,:,t) * (S(ndxOther,t) - Sposterior(ndxOther));
                        thisSig = sqrtOmega1(s,t);
                        S(s,t)  = drawTruncNormal(thisMu, thisSig, elbBound, udraws(s,t,n));
                    end
                    
                end
            end
            
            
            Y(ndxS,t) = S(:,t); % note: future values in YY need
                                % not be updated here
            
        end % any(sNaN(:,t))
        
        % update (even if not(sNaN))
        if t >= p
            % this code is quicker than the implementation for ...
            % t < p, used in the else clause below
            this            = Y(:,t:-1:t-p+1);
            STATElag(2:end) = this(:);
        else
            statelag = STATElag(1 + (1 : Ny * (p - 1)));
            STATElag = cat(1, 1, Y(:,t), statelag);
        end
    end % for t
    
    if n > burnin
        shadowrateDraws(:,:,n-burnin) = S;
    end
    
    %     progressbar(n / (totalNdraws))
end


