function [Xdraw, Xhat, P] = rwnoisePrecisionBasedSampler(Y, Ny, T, volRW, volNOISE, X0, sqrtV0, Ndraws, rndStream)
% abcdPrecisionBasedSampler computes smoothed kalman states using the stacked approach of 
% Chan and Jeliazkov
%  
%   ... 

% assumes volRW is Nx x 1 vector, and volNOISE is (Ny x T) x 1 vectors (i.e. no correlation within X and Y)


if nargin < 9
    Ndraws = 1;
end
if nargin < 10
    rndStream = getDefaultStream;
end

%% read parameters
Y        = Y(:);
NyT      = Ny * T;
Nx       = Ny;
NxT      = Nx * T; % obsolete but for better readability
NxTm1    = Nx * (T - 1); % obsolete but for better readability

%% construct stacked system
XX0 = sparse(1:Nx, 1, X0, Nx * T, 1);

rowndx = [1 : NxT, Nx + 1 : NxT];
colndx = [1 : NxT, 1 : NxTm1];
values = [ones(1,NxT), -ones(1, NxTm1)];
AA  = sparse(rowndx, colndx, values);

CC  = speye(NyT, NyT);


% sqrtSIGMA = sparse(1:NxT, 1:NxT, [sqrtV0(:)', repmat(volRW(:)', 1, T-1)]);
ITm1        = speye(T-1);
sqrtSIGMA   = blkdiag(sqrtV0,  kron(ITm1, volRW));

sqrtOMEGA   = sparse(1:NyT, 1:NyT, volNOISE);

%% set up  stacked system


AAtilde            = sqrtSIGMA \ AA;
XX0tilde           = sqrtSIGMA \ XX0;

CCtilde            = sqrtOMEGA \ CC;
Ytilde             = sqrtOMEGA \ Y;

P                   = AAtilde' * AAtilde + (CCtilde' * CCtilde);
[sqrtP, flag]       = chol(P, 'lower');

if flag > 0
    error('P not posdf, using QR instead')
    % via qr -- much slower
    M = [AAtilde; CCtilde]; %#ok<UNRCH>
    m = size(M,2);
    [~, R] = qr(M);
    sqrtP = R(1:m,1:m)';
    % checkdiff(sqrtP * sqrtP', sqrtP2 * sqrtP2');
end

sqrtPXhat   = sqrtP \ (AAtilde' * XX0tilde + CCtilde' * Ytilde); 

Zdraw        = randn(rndStream, Nx * T, Ndraws);
Xdraw        = (sqrtP') \ (sqrtPXhat + Zdraw);
Xdraw        = reshape(Xdraw, Nx, T, Ndraws);

if nargout > 1
    Xhat        = (sqrtP') \ sqrtPXhat;
    Xhat        = reshape(Xhat, Nx, T);
end

