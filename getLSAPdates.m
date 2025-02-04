function [lsapDates, lsapLabels] = getLSAPdates()
% GETLSAPDATES ... 
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 

lsapDates = cat(1, datetime(2008,11, 25), datetime(2010,11,3), datetime(2011, 9, 21), datetime(2012,9,13), ...
    datetime(2013,5,22), datetime(2013,12,18), datetime(2014,10,29), ...
    datetime(2020,3,20), datetime(2020,12,16), datetime(2022,3,16));
% lsapLabels = {'LSAP1', 'LSAP2', 'MEP', 'LSAP3', ...
%     'taper tantrum', 'taper begins', 'LSAP end', ...
%     'Pandemic QE', 'QE slows', 'QE stops'};
lsapLabels = {'QE1', 'QE2', 'Operation Twist', 'QE3', ...
    'Taper Tantrum', 'Taper begins', 'QE ends', ...
    'QE begins', 'QE slows', 'QE ends'};

