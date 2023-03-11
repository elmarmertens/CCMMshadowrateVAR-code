function llf = logscoreGaussianCensored(mu,sqrtOmega,y,ELBbound,ndxCensored)
% LOGSCOREGAUSSIANCENSORED ...
%
%   ...

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 03-Jan-2023 10:22:35 $
% $Revision : 1.00 $
% DEVELOPED : 9.12.0.2039608 (R2022a) Update 5
% FILENAME  : logscoreGaussianCensored.m

N = length(y);

if nargin < 5 || isempty(ndxCensored)
    ndxCensored = true(N,1);
elseif ~islogical(ndxCensored)
    ndxCensored = ismember(1:N, ndxCensored(:));
end

% ensure column vectors
ndxCensored = ndxCensored(:); 
y           = y(:);

ndxAtELB   = (y <= ELBbound) & ndxCensored;
ndxOffELB  = ~ndxAtELB;
NatELB     = sum(ndxAtELB);
% checkdiff(NatELB, sum(y(ndxCensored) <= ELBbound));
NoffELB    = N - NatELB;


if NatELB == 0

    llf = logscoreGaussian(muY, sqrtOmegaY, y);
    % this could also be caught at level of calling function to avoid the overhead

else % at least one element of y at ELB

    ydev     = y - mu;
    logtwopi = log(2 * pi);

    yAtELB   = y(ndxAtELB);
    muAtELB  = mu(ndxAtELB);

    %% reorder with unconstrained on top
    ndx       = 1:N;
    ndxOrder  = [ndx(ndxOffELB), ndx(ndxAtELB)];
    ydev      = ydev(ndxOrder);
    ndxOffELB = 1:NoffELB;
    ndxAtELB  = NoffELB+1:N;
    sqrtOmega = sqrtOmega(ndxOrder,:);
    Omega     = sqrtOmega * sqrtOmega';
    sqrtOmega = chol(Omega, 'lower');

    % some checks
    % mu        = mu(ndxOrder);
    % checkdiff(muAtELB, mu(ndxAtELB));
    % y         = y(ndxOrder);
    % checkdiff(yAtELB, y(ndxAtELB));

    %% compute likelihood of obs away from ELB
    if NoffELB > 1
        sqrtOmega11   = sqrtOmega(ndxOffELB,ndxOffELB);
        logdetOmega11 = 2 * sum(log(diag(sqrtOmega11)));
        z1dev         = sqrtOmega11 \ (ydev(ndxOffELB));
        llf1          = -.5 * (NoffELB * logtwopi + logdetOmega11 + sum(z1dev.^2));

        % prepare moments for computing likelihood at ELB
        sqrtOmega21 = sqrtOmega(ndxAtELB,ndxOffELB);
        sqrtOmega22 = sqrtOmega(ndxAtELB,ndxAtELB);
        y21         = muAtELB + sqrtOmega21 * z1dev;

    else
        llf1        = 0;
        sqrtOmega22 = sqrtOmega(ndxAtELB,ndxAtELB);
        y21         = muAtELB;
    end

    %% compute likleihood of obs at ELB
    if NatELB == 1
        llf2 = log(normcdf(yAtELB, y21, sqrtOmega22)); % note: normcdf takes vol as input
    else
        llf2 = log(mvncdf(yAtELB', y21', sqrtOmega22 * sqrtOmega22')); % note: mvncdf expects row vectors and VCV matrix (though results appear robust when inputting two column vectors)
        %         checkdiff(llf2, log(mvncdf(yAtELB(:), y21(:), sqrtOmega22 * sqrtOmega22')));
    end

    %% sum both together
    llf = llf1 + llf2;
end