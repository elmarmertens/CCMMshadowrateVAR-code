function llf = logscoreGaussian(mu,sqrtOmega,y,logdetOmega)
% LOGSCOREGAUSSIAN ... 
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 03-Jan-2023 10:22:35 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.12.0.2039608 (R2022a) Update 5 
% FILENAME  : logscoreGaussian.m 


if nargin < 4 || isempty(logdetOmega)
    logdetOmega = 2 * sum(log(diag(sqrtOmega)));
end

N = length(y);

ydev = sqrtOmega \ (y - mu);
llf  = -.5 * (N * log(2 * pi) + logdetOmega + sum(ydev.^2));
