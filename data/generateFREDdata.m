% generate_freddata.m
% =========================================================================
% DESCRIPTION:
% This script loads in raw data from a monthly database CSV file,
% transforms each series based on transformation code using
% prepare_missing.m, and removes outliers from the transformed data using
% remove_outliers.m.
%
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

for this = {...
        'fredsxMD20', 'fredsxMD20exYield', ...
        'fredsxMD14longyields', ...
        }

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

    data = yt;

    % =========================================================================
    % SELECT 20 VARIABLES AND STORE IN CSV

    switch datalabel
        case {'fredsxMD20'}
            codeVariableSelection = {'FEDFUNDS', ...
                'TB6MS', 'GS1', ...
                'GS5', 'GS10', ...
                'BAA', ...
                'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
                'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
                'PPICMM', 'PCEPI', 'HOUST', 'S&P 500', 'EXUSUKx'};
        case {'fredsxMD20exYield'}
            codeVariableSelection = {'FEDFUNDS', 'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
                'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
                'PPICMM', 'PCEPI', 'HOUST', 'S&P 500', 'EXUSUKx'};
        case {'fredsxMD14longyields'}
            codeVariableSelection = {'GS5', 'GS10', 'BAA', ...
                'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
                'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
                'PPICMM', 'PCEPI', 'HOUST', 'S&P 500', 'EXUSUKx', ...
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

    %% tabulate variable definitions
    varlabels = fredMDprettylabel(ncode);
    N = length(varlabels);

    filename = sprintf('datalist-%s.tex', outputlabel);
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

   
end

%% finish
dockAllFigures