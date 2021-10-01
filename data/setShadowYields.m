ndxSHADOWRATE = find(ismember(ncode, {'FEDFUNDS', 'TB3MS', 'TB6MS', 'GS1', ...
       'WUXIASHADOWRATE', 'KRIPPNERSHADOWRATE'}));
% define index of yields that need to obey ELB (at least out of sample)
ndxOTHERYIELDS = find(ismember(ncode, {'GS5', 'GS10', 'GS20', 'BAA'}));
