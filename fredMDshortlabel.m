function pcode = fredMDshortlabel(ncode)
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
            pcode{n} = 'Consumption';
        case 'CMRMTSPLx'
            pcode{n} = 'Sales';
        case 'INDPRO'
            pcode{n} = 'IP';
        case 'CUMFNS'
            pcode{n} = 'Cap. Util.';
        case 'UNRATE'
            pcode{n} = 'Unemp.';
        case 'PAYEMS'
            pcode{n} = 'Nfm Pyrlls';
        case 'CES0600000007'
            pcode{n} = 'Hours';
        case 'CES0600000008'
            pcode{n} = 'H. Earnings';
        case 'WPSFD49207'
            pcode{n} = 'PPI (Fin.)';
        case 'PPICMM'
            pcode{n} = 'PPI (Metals)';
        case 'PCEPI'
            pcode{n} = 'PCE Prices';
        case 'CPIAUCSL'
            pcode{n} = 'CPI';
        case 'FEDFUNDS'
            pcode{n} = 'FFR';
        case {'WUXIASHADOWRATE', 'KRIPPNERSHADOWRATE'}
            pcode{n} = 'Policy Rate';
        case 'TB3MS'
            pcode{n} = '3m T-Bill';
        case 'HOUST'
            pcode{n} = 'Hsng Strts';
        case {'S_P500', 'SP500'}
            pcode{n} = 'S\&P 500';
        case 'EXUSUKx'
            pcode{n} = 'USD / GBP';
            case 'TB6MS'
                pcode{n} = '6m Tbill';
            case 'BAA'
                pcode{n} = 'BAA Yld';
        case 'GS1'
            pcode{n} = '1y Trsy';
        case 'GS5'
            pcode{n} = '5y Trsy';
        case 'GS10'
            pcode{n} = '10y Trsy';
        case 'GS20'
            pcode{n} = '20y Trsy';
        case 'BAAFFM'
            pcode{n} = 'Baa Sprd';
        otherwise
            pcode{n}= ncode{n};
    end % switch
end % for n
