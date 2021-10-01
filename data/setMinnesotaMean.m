%% mean for Minnesota prior: zero (diff) or RW (level)

minnesotaPriorMean = NaN(N,1);

for n = 1 : N
    switch ncode{n}
        case {'CUMFNS', 'UNRATE', ...
                'WPSFD49207',  'PPICMM', 'PCEPI', ...
                'HOUST',  'BAAFFM', ...
                'BAA10Y', 'BAA', ...
                'FEDFUNDS', 'TB3MS', 'TB6MS', 'GS1', 'GS5', 'GS10', 'GS20', ...
                'WUXIASHADOWRATE', 'KRIPPNERSHADOWRATE'}
            minnesotaPriorMean(n) = 1;
        otherwise
            minnesotaPriorMean(n) = 0;
    end
end