% generate_freddata.m
% =========================================================================
% DESCRIPTION:
% This script loads in raw data from a monthly database CSV file,
% transforms each series based on transformation code using
% prepare_missing.m, and removes outliers from the transformed data using
% remove_outliers.m.
%
% NOTE:
% The default CSV file read by this code is 2015-04.csv, which is the April
% 2015 version of the dataset. If using a different version, make sure to
% change the variable "csv_in" on line 26 to match the name of the relevant
% CSV file.
%
% =========================================================================
% CLEAR:
clear
close all
clc

%#ok<*DATNM>
%#ok<*DATST>

%% load em toolboxes
path(pathdef)

addpath ../matlabtoolbox/emtools/
addpath ../matlabtoolbox/emtexbox/
addpath ../matlabtoolbox/emgibbsbox/
addpath ../matlabtoolbox/emeconometrics/
addpath ../matlabtoolbox/emstatespace/
addpath ..

doPlots = true;

%% settings

for this = {'fredMD20VXO', 'fredMD20VXOexYield', 'fredMD20EBP', 'fredMD20EBPexYield'}
%     {'fredMD20', 'fredMD20exYield', ...
%     'fredMD14longyields', 'fredMD15plus6M', 'fredMD15plus1Y', 'fredMD15plus5Y', 'fredMD15plus10Y', ...
%         'fredMD15plusShortYields',  'fredMD15plusLongYields', ...
%         'fredMD15plusBAA', 'fredMD15plusLongYieldsBAA'}

    clear rawdata tcode names tabledata

    datalabel = this{:};

    vintage = '2022-09';

    outputlabel = strcat(datalabel, '-', vintage);
    doQuarterly = false;

    %#ok<*UNRCH>

    % =========================================================================
    % PARAMETER TO BE CHANGED:
    % Update the .csv filename to match the desired version

    % CSV file name

    csv_in= strcat(vintage, '.csv');

    % =========================================================================
    % LOAD AND LABEL DATA:
    % Load data from CSV file
    dum=importdata(csv_in,',');

    % Variable names
    names = dum.textdata(1,2:end);

    % Transformation numbers
    tcode = dum.data(1,:);

    % Raw data
    rawdata = dum.data(2:end,:); % first row contains tcodes

    % Month of final observation
    % final_month=month(dum.textdata(end,1));

    % Year of final observation
    % final_year=year(dum.textdata(end,1));

    % =========================================================================
    % SET UP DATES:
    % Dates (monthly) are of the form YEAR+MONTH/12
    % e.g. March 1970 is represented as 1970+3/12
    % Dates go from 1959:01 to final_year:final_month (see above)
    % dates = (1959+1/12:1/12:final_year+final_month/12)';

    dates = datenum(dum.textdata(3:end,1), 'mm/dd/yy', 1950); % using matlab dates

    % T = number of months in sample
    T       = size(dates,1);
    rawdata = rawdata(1:T,:);

    %% patch in BAA10
    ndxBAA    = strcmp(names, 'BAA');
    ndx10Y    = strcmp(names, 'GS10');

    BAA10Y      = rawdata(:,ndxBAA) - rawdata(:,ndx10Y);
    ndxBAAFFM   = strcmp(names, 'BAAFFM');
    BAAFFM      = rawdata(:,ndxBAAFFM);

    if doPlots
        figure
        h = plot(dates, [BAA10Y BAAFFM]);
        nbershades(dates)
        legend(h, 'BAA10Y', 'BAAFM')
    end

    % ndxFFR      = strcmp(names, 'FEDFUNDS');
    % spreadcheck = rawdata(:,ndxBAA) - rawdata(:,ndxFFR);
    % checkdiff(spreadcheck, rawdata(:,ndxBAAFFM));

    % augment FREDMD
    rawdata = [rawdata, BAA10Y]; %#ok<AGROW>
    names   = cat(2, names, 'BAA10Y');
    tcode   = [tcode, 1];  %#ok<AGROW>

    %% do not first difference the fedfunds rate or yields
    ndxFFR        = strcmp(names, 'FEDFUNDS');
    tcode(ndxFFR) = 1;


    tcode(strcmp(names, 'TB3MS'))   = 1;
    tcode(strcmp(names, 'TB6MS'))   = 1;
    tcode(strcmp(names, 'GS1'))     = 1;
    tcode(strcmp(names, 'GS10'))    = 1;
    tcode(strcmp(names, 'GS5'))     = 1;
    tcode(strcmp(names, 'PCEPI'))   = 5;
    tcode(strcmp(names, 'UNRATE'))  = 1;

    % do not double difference inflation, not difference unrate
    tcode(tcode == 2) = 1; % takes also care of the interest rates listed above
    tcode(tcode == 6) = 5;

    % other
    if contains(datalabel, 'levels')
        tcode(tcode == 5) = 4;
    end



    %% patch in Krippner
    krippner         = importdata('SSR_Estimates_2022March.xlsx');
    krippnerDatesEOM = xls2mdate(krippner.data.D0x2EMonthlyAverageSSRSeries(:,2));
    ndxSSR = 3; % 3rd data column, 2nd text column
    if ~strcmpi('US SSR', krippner.textdata.D0x2EMonthlyAverageSSRSeries(20,ndxSSR - 1))
        error('check out Krippner XLS sheet')
    end
    krippnerRate = krippner.data.D0x2EMonthlyAverageSSRSeries(:,ndxSSR);

    [y, m] = datevec(krippnerDatesEOM);
    krippnerDates = datenum(y,m,1);

    ndxTB3MS   = strcmpi('tb3ms', names);
    if tcode(ndxTB3MS) ~= 1
        error('TB3MS not collected in levels?')
    end


    if doPlots
        figure
        hold on
        plot(dates, rawdata(:,ndxTB3MS), 'k-')
        plot(krippnerDates, krippnerRate, 'r--')
        xtickdates(dates)
    end

    shadowrate          = rawdata(:,ndxTB3MS);
    ndxKRIPPNERinFRED   = ismember(dates, krippnerDates);
    ndxKRIPPNER4FRED    = ismember(krippnerDates, dates);

    shadowrate(ndxKRIPPNERinFRED) = krippnerRate(ndxKRIPPNER4FRED);

    % augment FREDMD
    rawdata = [rawdata, shadowrate]; %#ok<AGROW>
    names   = cat(2, names, 'KRIPPNERSHADOWRATE');
    tcode   = [tcode, 1]; %#ok<AGROW>

    % store krippner/tb3 rate as separate csv for use in plotting
    thislabel = 'KRIPPNERSHADOWRATE';
    datatable   = array2table(shadowrate,  'VariableNames', {thislabel});
    output      = cat(2, table(dates, 'VariableNames', {'dates'}), datatable);
    writetable(output, sprintf('%s.csv', thislabel))

    %% patch in WuXia
    wuxia = importdata('shadowrate_US.xls');

    y = floor(wuxia(:,1) / 100);
    m = mod(wuxia(:,1),100);

    wuxiaDates = datenum(y,m,1);

    ndxTB3MS   = strcmpi('tb3ms', names);
    if tcode(ndxTB3MS) ~= 1
        error('TB3MS not collected in levels?')
    end


    if doPlots
        figure
        hold on
        plot(dates, rawdata(:,ndxTB3MS), 'k-')
        plot(wuxiaDates, wuxia(:,2), 'r--')
        xtickdates(dates)
    end

    shadowrate       = rawdata(:,ndxTB3MS);
    ndxWUXIAinFRED   = ismember(dates, wuxiaDates);
    ndxWUXIA4FRED    = ismember(wuxiaDates, dates);

    shadowrate(ndxWUXIAinFRED) = wuxia(ndxWUXIA4FRED,2);

    % augment FREDMD
    rawdata = [rawdata, shadowrate]; %#ok<AGROW>
    names   = cat(2, names, 'WUXIASHADOWRATE');
    tcode   = [tcode, 1]; %#ok<AGROW>

    thislabel = 'WUXIASHADOWRATE';
    datatable   = array2table(shadowrate,  'VariableNames', {thislabel});
    output      = cat(2, table(dates, 'VariableNames', {'dates'}), datatable);
    writetable(output, sprintf('%s.csv', thislabel))

    %% patch in G-Z EBP
    ebpfile = importdata('ebp_csv.csv');

    ebpdates = datenum(ebpfile.textdata(2:end,1));
    
    ebp        = NaN(size(rawdata,1), 1);
    [~, ndx1, ndx2] = intersect(dates, ebpdates);
    ebp(ndx1,1) = ebpfile.data(ndx2,2); 

    % augment FREDMD
    rawdata = [rawdata, ebp]; %#ok<AGROW>
    names   = cat(2, names, 'EXCESSBONDPREMIUM');
    tcode   = [tcode, 1]; %#ok<AGROW>

    thislabel   = 'EXCESSBONDPREMIUM';
    datatable   = array2table(ebp,  'VariableNames', {thislabel});
    output      = cat(2, table(dates, 'VariableNames', {'dates'}), datatable);
    writetable(output, sprintf('%s.csv', thislabel))

    %% transform into quarterly data if needed
    if doQuarterly

        quarters = quarterlydates(dates);



        qdata = grpstats(rawdata, quarters);
        qdates = unique(quarters);

        % drop incomplete quarters
        qcount = arrayfun(@(x) sum(quarters == x), qdates);
        ndx    = qcount == 3;

        rawdata = qdata(ndx,:);
        dates   = qdates(ndx,:);

        T       = size(dates,1);

        if any(~ndx)
            warning('dropped %d monthly data points that belong to incomplete quarters', sum(~ndx))
        end

    end

    %% TRANSFORM RAW DATA INTO STATIONARY FORM:
    % Use function prepare_missing.m
    %   Output yt: matrix containing data after transformation
    %
    %     case 1, % Level (i.e. no transformation): x(t)
    %     case 2, % First difference: x(t)-x(t-1)
    %     case 3, % Second difference: (x(t)-x(t-1))-(x(t-1)-x(t-2))
    %     case 4, % Natural log: ln(x)
    %     case 5, % First difference of natural log: ln(x)-ln(x-1)
    %     case 6, % Second difference of natural log: (ln(x)-ln(x-1))-(ln(x-1)-ln(x-2))
    %     case 7, % First difference of percent change: (x(t)/x(t-1)-1)-(x(t-1)/x(t-2)-1)



    yt  = prepare_missing(rawdata,tcode);
    ndx = tcode == 5;
    if doQuarterly
        yt(:,ndx) = yt(:,ndx) * 400;
    else
        yt(:,ndx) = yt(:,ndx) * 1200;
    end

    % =========================================================================
    % REDUCE SAMPLE TO USABLE DATES:
    % Remove first two months because some series have been second differenced
    yt=yt(3:T,:);
    dates=dates(3:T,:);
    % T = length(dates);

    % =========================================================================
    % REMOVE OUTLIERS: EM dropped
    % Use function remove_outliers.m (see for definition of outliers)
    %   Output data: matrix containing transformed series after removal of outliers
    %   Output n: matrix containing number of outliers removed for each series
    % [~,n]=remove_outliers(yt);
    %
    % disp('there are quite a few outliers:')
    % display(n)

    data = yt;

    % =========================================================================
    % SELECT 20 VARIABLES AND STORE IN CSV

    switch datalabel
        case {'fredMD3'}
            codeVariableSelection = {'UNRATE', 'PCEPI', 'FEDFUNDS'};
        case {'fredMD3cumYield'}
            codeVariableSelection = {'UNRATE', 'PCEPI', 'FEDFUNDS', 'TB6MS', 'GS1', 'GS5', 'GS10'};
        case {'fredMD3krippner'}
            codeVariableSelection = {'UNRATE', 'PCEPI', 'KRIPPNERSHADOWRATE'};
        case {'fredMD3wuxia'}
            codeVariableSelection = {'UNRATE', 'PCEPI', 'WUXIASHADOWRATE'};
        case {'fredMD20'}
            codeVariableSelection = {'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
                'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
                'PPICMM', 'PCEPI', 'FEDFUNDS', 'HOUST', 'S&P 500', 'EXUSUKx', ...
                'TB6MS', 'GS1', ...
                'GS5', 'GS10', ...
                'BAA'};
        case {'fredMD20exYield'}
            codeVariableSelection = {'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
                'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
                'PPICMM', 'PCEPI', 'FEDFUNDS', 'HOUST', 'S&P 500', 'EXUSUKx'};
        case {'fredMD20VXO'}
            codeVariableSelection = {'VIXCLSx', 'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
                'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
                'PPICMM', 'PCEPI', 'FEDFUNDS', 'HOUST', 'S&P 500', 'EXUSUKx', ...
                'TB6MS', 'GS1', ...
                'GS5', 'GS10', ...
                'BAA'};
        case {'fredMD20VXOexYield'}
            codeVariableSelection = {'VIXCLSx', 'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
                'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
                'PPICMM', 'PCEPI', 'FEDFUNDS', 'HOUST', 'S&P 500', 'EXUSUKx'};
        case {'fredMD20EBP'}
            codeVariableSelection = {'EXCESSBONDPREMIUM', 'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
                'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
                'PPICMM', 'PCEPI', 'FEDFUNDS', 'HOUST', 'S&P 500', 'EXUSUKx', ...
                'TB6MS', 'GS1', ...
                'GS5', 'GS10', ...
                'BAA'};
        case {'fredMD20EBPexYield'}
            codeVariableSelection = {'EXCESSBONDPREMIUM', 'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
                'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
                'PPICMM', 'PCEPI', 'FEDFUNDS', 'HOUST', 'S&P 500', 'EXUSUKx'};
        case {'fredMD14longyields'}
            codeVariableSelection = {'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
                'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
                'PPICMM', 'PCEPI', 'HOUST', 'S&P 500', 'EXUSUKx', ...
                'GS5', 'GS10', ...
                'BAA'};
        case {'fredMD15plus6M'}
            codeVariableSelection = {'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
                'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
                'PPICMM', 'PCEPI', 'FEDFUNDS', 'HOUST', 'S&P 500', 'EXUSUKx', ...
                'TB6MS'};
        case {'fredMD15plus1Y'}
            codeVariableSelection = {'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
                'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
                'PPICMM', 'PCEPI', 'FEDFUNDS', 'HOUST', 'S&P 500', 'EXUSUKx', ...
                'GS1'};
        case {'fredMD15plus5Y'}
            codeVariableSelection = {'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
                'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
                'PPICMM', 'PCEPI', 'FEDFUNDS', 'HOUST', 'S&P 500', 'EXUSUKx', ...
                'GS5'};
        case {'fredMD15plus10Y'}
            codeVariableSelection = {'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
                'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
                'PPICMM', 'PCEPI', 'FEDFUNDS', 'HOUST', 'S&P 500', 'EXUSUKx', ...
                'GS10'};
        case {'fredMD15plusShortYields'}
            codeVariableSelection = {'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
                'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
                'PPICMM', 'PCEPI', 'FEDFUNDS', 'HOUST', 'S&P 500', 'EXUSUKx', ...
                'TB6MS', 'GS1'};
        case {'fredMD15plusLongYields'}
            codeVariableSelection = {'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
                'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
                'PPICMM', 'PCEPI', 'FEDFUNDS', 'HOUST', 'S&P 500', 'EXUSUKx', ...
                'GS5', 'GS10'};
        case {'fredMD15plusBAA'}
            codeVariableSelection = {'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
                'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
                'PPICMM', 'PCEPI', 'FEDFUNDS', 'HOUST', 'S&P 500', 'EXUSUKx', ...
                'BAA'};
        case {'fredMD15plusLongYieldsBAA'}
            codeVariableSelection = {'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
                'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
                'PPICMM', 'PCEPI', 'FEDFUNDS', 'HOUST', 'S&P 500', 'EXUSUKx', ...
                'GS5', 'GS10', 'BAA'};
        case {'fredMDmacroFFRGS10'}
            codeVariableSelection = {'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
                'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
                'PPICMM', 'PCEPI', 'FEDFUNDS', 'HOUST', 'S&P 500', 'EXUSUKx', 'GS10'};
        case {'fredMD20krippner'}
            codeVariableSelection = {'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
                'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
                'PPICMM', 'PCEPI', ...
                'KRIPPNERSHADOWRATE', ...
                'HOUST', 'S&P 500', 'EXUSUKx', ...
                'TB6MS', 'GS1', ...
                'GS5', 'GS10', ...
                'BAA'};
        case {'fredMD20wuxia'}
            codeVariableSelection = {'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
                'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
                'PPICMM', 'PCEPI', ...
                'WUXIASHADOWRATE', ...
                'HOUST', 'S&P 500', 'EXUSUKx', ...
                'TB6MS', 'GS1', ...
                'GS5', 'GS10', ...
                'BAA'};
        case {'fredMD20exYieldKrippner'}
            codeVariableSelection = {'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
                'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
                'PPICMM', 'PCEPI', ...
                'KRIPPNERSHADOWRATE', ...
                'HOUST', 'S&P 500', 'EXUSUKx' ...
                };
        case {'fredMD20exYieldWuxia'}
            codeVariableSelection = {'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
                'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
                'PPICMM', 'PCEPI', ...
                'WUXIASHADOWRATE', ...
                'HOUST', 'S&P 500', 'EXUSUKx', ...
                };
        otherwise
            error('datalabel %s not recognized', datalabel)
    end



    % map FRED-MD into list
    [~,ndxVariableSelection] = ismember(codeVariableSelection, names);

    % collect data and names
    ncode  = names(ndxVariableSelection);
    tcode  = tcode(ndxVariableSelection);

    cumcode   = false(1,length(ndxVariableSelection));
    cumcode(tcode == 2) = true;
    % cumcode(tcode == 5) = true; per Todd's suggestion from Dec 9 2019
    cumcode(tcode == 6) = true;

    tabledata = data(:,ndxVariableSelection);



    ncode = matlab.lang.makeValidName(ncode, 'ReplacementStyle', 'delete'); % to get rid of '&' and other unwanted signs, seems necessary at least under linux

    %% check for outliers
    tableNoOutliers = remove_outliers(tabledata);
    ndx = isnan(tableNoOutliers);
    for n = 1 : size(tabledata,2)
        thisdata = tableNoOutliers(:,n);
        nanny = isnan(thisdata);
        % if any(nanny)
        theseOutliers = NaN(size(tabledata,1),1);
        theseOutliers(nanny) = tabledata(nanny,n);

        if doPlots
            figure
            hold on
            plot(dates, thisdata)
            if any(nanny)
                plot(dates, theseOutliers, 'rx', 'linewidth', 2)
                plot(dates, theseOutliers, 'ro', 'linewidth', 2)
            end
            xtickdates(dates)
            titletxt = sprintf('%s (%d)', ncode{n}, tcode(n));
            title(titletxt)
            set(gcf, 'name', titletxt)
        end
    end

    %% clean missing values
    nanny = any(isnan(tabledata), 2);

    if ~iscompact(nanny)
        error('missing data inside sample')
    end

    if any(nanny)
        nannyCell = arrayfun(@(x) datestr(x), dates(nanny), 'UniformOutput', false);
        warning('data is missing for %s\n', nannyCell{:})
        %     arrayfun(@(x) warning('missing data at %s', datestr(x)), dates(nanny));
        warning('data is missing for %s\n', ncode{any(isnan(tabledata), 1)})
    end

    if ~iscompact(nanny)
        error('date vector not compact after pruning NaN')
    end

    %% prepare table
    tabledata = tabledata(~nanny,:);
    dates     = dates(~nanny);
    % T         = length(dates);

    N     = length(ncode);

    %% construct yield code
    setShadowYields;
    ndxYIELDS = union(ndxSHADOWRATE, ndxOTHERYIELDS);
    ndxYIELDS = ismember(1:N,ndxYIELDS);

    %% prepend dates, tcode etc and store table
    if doQuarterly
        datalabel = strcat(datalabel, '-quarterly');
    end


    tabledata = [tcode; cumcode; tabledata]; %#ok<AGROW>

    datatable   = array2table(tabledata,  'VariableNames', ncode);

    % check
    if any(any(ismissing(datatable)))
        error('there are missing observations')
    end

    tabledates  = [NaN;NaN;dates]; % note: recycle the variable name dates
    output = cat(2, table(tabledates, 'VariableNames', {'dates'}), datatable);
    writetable(output, sprintf('%s.csv', outputlabel))

    %% define minnesota prior means
    setMinnesotaMean


    %% prepare fortran files
    % % output as txt file for fortran
    % mat2fortran(sprintf('%s.yData.txt', outputlabel), tabledata(3:end,:))
    % % mat2fortran(sprintf('%s.cumcode.txt', outputlabel), tabledata(2,:))
    % int2fortran(sprintf('%s.tcode.txt', outputlabel), tabledata(1,:))
    % mat2fortran(sprintf('%s.minnesotaPriorMean.txt', outputlabel), minnesotaPriorMean)
    % logical2fortran(sprintf('%s.ndxYIELDS.txt', outputlabel), ndxYIELDS)
    % mat2fortran(sprintf('%s.dates.txt', outputlabel), dates)
    %
    % % output ncode
    % cellstr2fortran(sprintf('%s.ynames.txt', outputlabel), ncode)
    % cellstr2fortran(sprintf('%s.ylabel.txt', outputlabel), fredMDprettylabel(ncode))

    %% generates table of variables
    varlabels = fredMDprettylabel(ncode);

    filename = sprintf('datalist-%s.tex', outputlabel);
    fid = fopen(filename, 'wt');

    fprintf(fid, '\\begin{center}\n');
    fprintf(fid, '\\begin{tabular}{lll}\n');
    fprintf(fid, '\\toprule\n');
    fprintf(fid, 'Variable & FRED-MD code & tcode ');
    fprintf(fid, '\\\\\n');
    fprintf(fid, '\\midrule\n');
    for n = 1 : N
        fprintf(fid, '%s ', varlabels{n});
        fprintf(fid, '& %s ', ncode{n});
        fprintf(fid, '& %d ', tcode(n));
        fprintf(fid, '\\\\\n');
    end
    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{center}\n');
    fprintf(fid, '\n');
    fprintf(fid, 'Note: ');
    fprintf(fid, 'Data obtained from the %s vintage of FRED-MD. ', strtok(csv_in, '.'));
    if doQuarterly
        fprintf(fid, 'Quarterly observations (constructed from monthly averages) ');
    else
        fprintf(fid, 'Monthly observations ');
    end
    % fprintf(fid, 'from %s to %s. ', datestryymm(dates(1)), datestryymm(dates(end)));
    fprintf(fid, 'from %s to %s. ', datestr(dates(1), 'yyyy:mm'), datestr(dates(end), 'yyyy:mm'));
    fprintf(fid, '%s \n', fredMDtcodeNote);
    fclose(fid);
    type(filename)

    %% generates alt table of variables: expressive transformation
    varlabels = fredMDprettylabel(ncode);
    N = length(varlabels);

    filename = sprintf('datalist2-%s.tex', outputlabel);
    fid = fopen(filename, 'wt');

    fprintf(fid, '\\begin{center}\n');
    fprintf(fid, '\\begin{tabular}{lll}\n');
    fprintf(fid, '\\toprule\n');
    fprintf(fid, 'Variable & FRED-MD code & transformation ');
    fprintf(fid, '\\\\\n');
    fprintf(fid, '\\midrule\n');
    for n = 1 : N
        fprintf(fid, '%s ', varlabels{n});
        fprintf(fid, '& %s ', ncode{n});

        switch tcode(n)
            case 1
                % % no transformation
            case 2
                fprintf(fid, '& %s\n', '\ensuremath{\Delta x_t}');
            case 4
                fprintf(fid, '& %s\n', '\ensuremath{\log(x_t)}');
            case 5
                if doQuarterly
                    fprintf(fid, '& %s\n', '\ensuremath{\Delta\log(x_t) \cdot 400}');
                else
                    fprintf(fid, '& %s\n', '\ensuremath{\Delta\log(x_t) \cdot 1200}');
                end
            otherwise
                error houston
        end

        fprintf(fid, '\\\\\n');
    end
    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{center}\n');
    fprintf(fid, '\n');
    fprintf(fid, 'Note: ');
    fprintf(fid, 'Data obtained from the %s vintage of FRED-MD. ', strtok(csv_in, '.'));
    if doQuarterly
        fprintf(fid, 'Quarterly observations (constructed from monthly averages) ');
    else
        fprintf(fid, 'Monthly observations ');
    end
    fprintf(fid, 'from %s to %s. ', datestr(dates(1), 'yyyy:mm'), datestr(dates(end), 'yyyy:mm'));
    % fprintf(fid, '%s \n', fredMDtcodeNote);
    fclose(fid);
    type(filename)

    %% generates alt table of variables: transformation and minnesota mean
    varlabels = fredMDprettylabel(ncode);
    N = length(varlabels);

    filename = sprintf('datalist3-%s.tex', outputlabel);
    fid = fopen(filename, 'wt');

    fprintf(fid, '\\begin{center}\n');
    fprintf(fid, '\\begin{tabular}{lllc}\n');
    fprintf(fid, '\\toprule\n');
    fprintf(fid, 'Variable & FRED-MD code & transformation & Minnesota prior');
    fprintf(fid, '\\\\\n');
    fprintf(fid, '\\midrule\n');
    for n = 1 : N
        fprintf(fid, '%s ', varlabels{n});
        fprintf(fid, '& %s ', ncode{n});

        switch tcode(n)
            case 1
                fprintf(fid, ' & ');
            case 2
                fprintf(fid, ' & %s', '\ensuremath{\Delta x_t}');
            case 4
                fprintf(fid, ' & %s', '\ensuremath{\log(x_t)}');
            case 5
                if doQuarterly
                    fprintf(fid, '& %s\n', '\ensuremath{\Delta\log(x_t) \cdot 400}');
                else
                    fprintf(fid, '& %s\n', '\ensuremath{\Delta\log(x_t) \cdot 1200}');
                end
            otherwise
                fprintf(fid, ' & ');
        end

        fprintf(fid, ' & %d ', minnesotaPriorMean(n));

        fprintf(fid, '\\\\\n');
    end
    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{center}\n');
    fprintf(fid, '\n');
    fprintf(fid, 'Note: ');
    fprintf(fid, 'Data obtained from the %s vintage of FRED-MD. ', strtok(csv_in, '.'));
    if doQuarterly
        fprintf(fid, 'Quarterly observations (constructed from monthly averages) ');
    else
        fprintf(fid, 'Monthly observations ');
    end
    fprintf(fid, 'from %s to %s.\n', datestr(dates(1), 'yyyy:mm'), datestr(dates(end), 'yyyy:mm'));
    fprintf(fid, 'Entries in the column ``Minnesota prior'''' report the prior mean on the first own-lag coefficient of the corresponding variable in each BVAR. Prior means on all other VAR coefficients are set to zero.\n');
    % fprintf(fid, '%s \n', fredMDtcodeNote);
    fclose(fid);
    type(filename)

    %% generates some stuff to collect latex tables
    varlabels = fredMDprettylabel(ncode);
    N = length(varlabels);

    filename = 'scratch.tex';
    fid = fopen(filename, 'wt');

    fprintf(fid, 'LABELS:\n');
    for n = 1 : N
        fprintf(fid, '%s\n', varlabels{n});
    end
    fprintf(fid, '\n');

    fprintf(fid, 'NAMES:\n');
    for n = 1 : N
        fprintf(fid, '%s\n', ncode{n});
    end

    fprintf(fid, 'FLOAT:\n');
    for n = 1 : N
        fprintf(fid, '\\subfloat[%s]{\\includegraphics[width=\\picwid]{pdf%s-\\jumpoff}\\label{subfig:%s-\\jumpoff}}\n', varlabels{n}, ncode{n}, ncode{n});
        if mod(n,3) == 0
            fprintf(fid, '\\\\\n');
        else
            fprintf(fid, '\\quad\n');
        end
    end




    fclose(fid);
    type(filename)

end

%% finish
dockAllFigures