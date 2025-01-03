function pcode = fredMDveryshortlabel(ncode)
% FREDMDPRETTYLABEL ...
%
%   ...

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 09-Dec-2019 18:55:22 $
% $Revision : 1.00 $
% DEVELOPED : 9.7.0.1247435 (R2019b) Update 2
% FILENAME  : fredmdPrettylabel.m



pcode = cell(size(ncode));

for n = 1 : length(ncode)
    switch ncode{n}
        case 'RPI'
            pcode{n} = 'Income';
        case 'DPCERA3M086SBEA'
            pcode{n} = 'Cons';
        case 'CMRMTSPLx'
            pcode{n} = 'Sales';
        case 'INDPRO'
            pcode{n} = 'IP';
        case 'CUMFNS'
            pcode{n} = 'CapUtil';
        case 'UNRATE'
            pcode{n} = 'Unemp';
        case 'PAYEMS'
            pcode{n} = 'Nfm Pyr';
        case 'CES0600000007'
            pcode{n} = 'Hours';
        case 'CES0600000008'
            pcode{n} = 'H. Earn';
        case 'WPSFD49207'
            pcode{n} = 'PPI Fin';
        case 'PPICMM'
            pcode{n} = 'PPI Met';
        case 'DSERRA3M086SBEA'
            pcode{n} = 'PCE Services Cons.';
        case 'DSERRG3M086SBEA'
            pcode{n} = 'PCE Services Prices';
        case 'PCEPI'
            pcode{n} = 'PCE Inf';
        case 'PCEPILFE'
            pcode{n} = 'PCE Core Prices';
        case 'CPIAUCSL'
            pcode{n} = 'CPI';
        case 'FEDFUNDS'
            pcode{n} = 'FFR';
        case {'WUXIASHADOWRATE', 'KRIPPNERSHADOWRATE'}
            pcode{n} = 'Policy Rate';
        case 'TB3MS'
            pcode{n} = '3m T-Bill';
        case 'HOUST'
            pcode{n} = 'Hsng St';
        case {'S_P500', 'SP500'}
            pcode{n} = 'S\&P500';
        case 'EXUSUKx'
            pcode{n} = 'FX \$/\pounds';
        case 'TB6MS'
            pcode{n} = '6m Tsy';
        case 'BAA'
            pcode{n} = 'BAA';
        case 'GS1'
            pcode{n} = '1y Tsy';
        case 'GS5'
            pcode{n} = '5y Tsy';
        case 'GS10'
            pcode{n} = '10y Tsy';
        case 'GS20'
            pcode{n} = '20y Tsy';
        case 'BAAFFM'
            pcode{n} = 'Baa Sprd';
        case 'EXCESSBONDPREMIUM'
            pcode{n} = 'EBP';

        otherwise
            pcode{n}= ncode{n};
    end % switch
end % for n
