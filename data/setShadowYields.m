if exist('ELBbound', 'var') && ELBbound > 0.25
    ndxSHADOWRATE = find(ismember(ncode, {'FEDFUNDS', 'TB3MS', 'TB6MS', 'GS1', ...
        'WUXIASHADOWRATE', 'KRIPPNERSHADOWRATE','GS5'}));
    % define index of yields that need to obey ELB (at least out of sample)
    ndxOTHERYIELDS = find(ismember(ncode, {'GS10', 'GS20', 'BAA'}));
else
    ndxSHADOWRATE = find(ismember(ncode, {'FEDFUNDS', 'TB3MS', 'TB6MS', 'GS1', ...
        'WUXIASHADOWRATE', 'KRIPPNERSHADOWRATE'}));
    % define index of yields that need to obey ELB (at least out of sample)
    ndxOTHERYIELDS = find(ismember(ncode, {'GS5', 'GS10', 'GS20', 'BAA'}));
end

ndxYIELDS = union(ndxSHADOWRATE, ndxOTHERYIELDS);

ndxFINANCIALS  = union(ndxYIELDS, find(ismember(ncode, {'SP500', 'EXUSUKx'})));