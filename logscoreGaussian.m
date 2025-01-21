function llf = logscoreGaussian(mu,sqrtOmega,y,logdetOmega)
% LOGSCOREGAUSSIAN ... 
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 


if nargin < 4 || isempty(logdetOmega)
    logdetOmega = 2 * sum(log(diag(sqrtOmega)));
end

N = length(y);

ydev = sqrtOmega \ (y - mu);
llf  = -.5 * (N * log(2 * pi) + logdetOmega + sum(ydev.^2));
