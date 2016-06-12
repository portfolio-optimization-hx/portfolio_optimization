%% portfolio optmization example script
%  this is an example script portfolio optmization modelling using code in
%  this repo. the scripts downloads the data and runs three portfolio
%  optmizations:
%       
%       mean-variance optmization
%       multi-period mean-variance optimization
%       mean-cvar otpmization
%
%  after the optimization is completed, the portfolio with the highest return 
%  to risk measure is selected from the efficient frontier. the portfolio
%  is then benchmarked against an equally weighted portfolio. see comments
%  in script and function for additional insights, this codebase is a
%  work-in-progress. 

% settings
data_dname      = 'data';
dataset_name    = 'djia';

train_sdate     = '1986/01/01'; % yyyy/mm/dd train start date
train_edate     = '1995/12/31'; % yyyy/mm/dd train end date

validate_sdate  = '1996/01/01'; % yyyy/mm/dd validate / cross-validate start date
validate_edate  = '2005/12/31'; % yyyy/mm/dd validate / cross-validate end date

test_sdate      = '2006/01/01'; % yyyy/mm/dd test start date
test_edate      = '2015/12/31'; % yyyy/mm/dd test end date
% settings

%% setup
train_drange    = [datenum(train_sdate,'yyyy/mm/dd') datenum(train_edate,'yyyy/mm/dd')];
validate_drange = [datenum(validate_sdate,'yyyy/mm/dd') datenum(validate_edate,'yyyy/mm/dd')];
test_drange     = [datenum(test_sdate,'yyyy/mm/dd') datenum(test_edate,'yyyy/mm/dd')];

asdata_fpath    = [pwd filesep data_dname filesep dataset_name '_asdata.mat'];
irdata_fpath    = [pwd filesep data_dname filesep 'interest_rates.mat'];

% add functions directory to path
if isempty(strfind(path,[pwd filesep 'functions']))
    addpath([pwd filesep 'functions']);
end

%% load, create, format data 
% check if portfolio file exist, else create it
if exist(asdata_fpath,'file')
    asdata = load(asdata_fpath);
    asdata = asdata.asdata;
else
    fprintf('Missing portfolio data file, attempt to create from csv.\n');
    asdata_dname = [pwd filesep data_dname filesep dataset_name];
    
    % check if asset directory exist and has at least two files, else download it
    if ~exist(asdata_dname,'dir') || (numel(dir(asdata_dname)) < 4) % 4 since 2 listing will be '.' '..' directories
        fprintf('Downloading asset data ...\n\n');
        download_data(data_dname,dataset_name);
        fprintf('\n');
    end
    if ~exist(asdata_dname,'dir')
        fprintf('Asset data does not exist and failure to download data.\n');
        return;
    end
    
    % get csv file paths for ascsv_to_asdata
    aslist_fpaths = dir(asdata_dname);
    aslist_fpaths = cellfun(@(c1) ([asdata_dname filesep c1]),{aslist_fpaths(:).name},'uniformoutput',false);
    if numel(aslist_fpaths) < 4
        fprintf('''Portfolio'' optimization unnecessary for single asset.\n');
        return;
    end
    
    % read format csv file to asdata struct
    asdata = ascsv_to_asdata(aslist_fpaths);
    if ~isempty(asdata)
        save(asdata_fpath,'asdata','-v7.3');
    else
        fprintf('Assets list data does not exist and failure to create data.\n');
    end
end

% check if interest / risk-free rate data, else create zero return irdata
if exist(irdata_fpath,'file')
    irdata = load(irdata_fpath);
    irdata = irdata.irdata;
else
    fprintf('Interest / risk-free rate data does not exit, a copy of the file is available on Github.\n');
    irdata = struct();
    irdata.timestamp = asdata.timestamp;
    irdata.annual = zeros(size(irdata.timestamp,1),4);
    irdata.daily  = irdata.annual;
end

%% mean-variance portfolio optimization
drange    = train_drange;
ret_linma = price_to_ma_return(asdata,drange,'linear'); % linear mean annual return
ret_covm  = price_to_covmatrix(asdata,drange); % return covariance matrix
wpa_min   = 0.00; % weight per asset max
wpa_max   = 0.10; % weight per asset min
% when asset does not contain price data for the period replace NaN
% with following values and optimization will not include the asset(s)
ret_linma(isnan(ret_linma)) = 0;
ret_covm(isnan(ret_covm))   = 1;

% mean-variance optimization function with: return mean, return
% standard-deviation as inputs, and min and max weight per asset as
% optional inputs
[eff_frontier,eff_weights] = mean_variance_optimization(ret_linma,ret_covm,wpa_min,wpa_max);
[~,idx1]    = max(eff_frontier(:,2) ./ eff_frontier(:,1)); % select portfolio with highest mean / variance ratio
w1          = eff_weights(idx1,:);
w2          = zeros(size(w1)) + (1 / asdata.count);

% calculate optimal portfolio and benchmark portfolio metrics
drc1 = {train_drange validate_drange test_drange}; % date range cell
for i1 = 1:numel(drc1) % loop through all date range
    dr = drc1{i1};
    
    % create investment schedule matrix, invest $10,000 at start date
    inschedule = [dr(1) 1e4]; 
    
    % calc_portfolio_return takes two asset data, date range, asset
    % weights, and investment schedule as inputs, 'real' vs 'theoretical'
    % and interest rate data as optional inputs to calculate portfolio
    % values. the functions defaults to 'real' calculation method, if
    % interest rate data is not provided, interest rate defaults to zero
    pf1 = calc_portfolio_return(asdata,dr,w1,inschedule,'real',irdata);
    pf2 = calc_portfolio_return(asdata,dr,w2,inschedule,'real',irdata); 
    
    % print results
    fprintf('%d-%d PORTFOLIO A - SELECTED PORTFOLIO\n\n',year(dr));
    str1_title = sprintf(' %9s %9s %9s %9s %9s\n','YEAR','RFR','RET','STDEV','SHARPE');
    str2 = sprintf(' %9d %9.3f %9.3f %9.3f %9.3f\n',[year(pf1.metrics_yearly(:,1)) pf1.metrics_yearly(:,2:end)]');
    fprintf('%s%s\n\n',str1_title,str2);
    
    fprintf('%9s PORTFOLIO B - BENCHMARK PORTFOLIO\n\n','');
    str2 = sprintf(' %9d %9.3f %9.3f %9.3f %9.3f\n',[year(pf1.metrics_yearly(:,1)) pf1.metrics_yearly(:,2:end)]');
    fprintf('%s%s\n\n',str1_title,str2);
end 

%% multi-period mean-variance portfolio optimization 
% calculate multi-period mean annual return and covariance matrix
drange    = train_drange;

eday = drange(2); % same end date
sday = (drange(1):drange(2)); % different start dates
idx1 = floor(weeknum(sday)/2); % find start of every 2 weeks
idx1 = find(idx1(2:end) ~= idx1(1:end-1)) + 1;
sday = sday([1 idx1]);
sday(eday - sday < 252) = []; % remove date range less 252 days since value calculated would not be as accurate

ret_linma = zeros(1,asdata.count); % allocate mean annual return matrix
ret_covm  = zeros(asdata.count); % allocate covariance matrix matrix
fprintf('Calculating multi-period returns and covariance matrix.\n');
for i1 = 1:numel(sday)
    ret_linma = ret_linma + price_to_ma_return(asdata,[sday(i1) eday],'linear'); % calculate and sum mean annual return for date range
    ret_covm  = ret_covm  + price_to_covmatrix(asdata,[sday(i1) eday]); % calculate and sum covariance matrix for date range
end

ret_linma = ret_linma / numel(sday); % find average mean annual return of all date range
ret_covm  = ret_covm  / numel(sday); % find average covariance matrix of all date range
wpa_min   = 0.00; % weight per asset max
wpa_max   = 0.10; % weight per asset min

% when asset does not contain price data for the period replace NaN
% with following values and optimization will not include the asset(s)
ret_linma(isnan(ret_linma)) = 0;
ret_covm(isnan(ret_covm))   = 1;

% mea-variance optimization function with: return mean, return
% standard-deviation as inputs, and min and max weight per asset as
% optional inputs
[eff_frontier,eff_weights] = mean_variance_optimization(ret_linma,ret_covm,wpa_min,wpa_max);
[~,idx1]    = max(eff_frontier(:,2) ./ eff_frontier(:,1)); % select portfolio with highest mean / variance ratio
w1          = eff_weights(idx1,:);
w2          = zeros(size(w1)) + (1 / asdata.count); % benchmark equally weighted portfolio

% calculate optimal portfolio and benchmark portfolio metrics
drc1 = {train_drange validate_drange test_drange}; % date range cell
for i1 = 1:numel(drc1) % loop through all date range
    dr = drc1{i1};
    
    % create investment schedule matrix, invest $1,000 every 2 weeks
    inschedule  = (dr(1):dr(2))';
    tfi1 = floor(weeknum(inschedule)/2);
    tfi1 = find(tfi1(2:end) ~= tfi1(1:end-1)) + 1;
    inschedule  = inschedule([1; tfi1]);
    inschedule(weekday(inschedule) == 1) = inschedule(weekday(inschedule) == 1) + 1;
    inschedule(:,2) = 1e3;
    inschedule(inschedule(:,1) < asdata.timestamp(1),1)   = asdata.timestamp(1);
    inschedule(inschedule(:,1) > asdata.timestamp(end),:) = [];
    
    pf1 = calc_portfolio_return(asdata,dr,w1,inschedule,'real',irdata);
    pf2 = calc_portfolio_return(asdata,dr,w2,inschedule,'real',irdata); 
    
    % print results
    fprintf('%d-%d PORTFOLIO A - SELECTED PORTFOLIO\n\n',year(dr));
    str1_title = sprintf(' %9s %9s %9s %9s %9s\n','YEAR','RFR','RET','STDEV','SHARPE');
    str2 = sprintf(' %9d %9.3f %9.3f %9.3f %9.3f\n',[year(pf1.metrics_yearly(:,1)) pf1.metrics_yearly(:,2:end)]');
    fprintf('%s%s\n\n',str1_title,str2);
    
    fprintf('%9s PORTFOLIO B - BENCHMARK PORTFOLIO\n\n','');
    str2 = sprintf(' %9d %9.3f %9.3f %9.3f %9.3f\n',[year(pf1.metrics_yearly(:,1)) pf1.metrics_yearly(:,2:end)]');
    fprintf('%s%s\n\n',str1_title,str2);
end 

%% mean-cvar portfolio optimization
drange      = train_drange;

beta        = 0.99; % percent interval to calculate cVaR
wpa_min     = 0.00; % weight per asset max
wpa_max     = 0.10; % weight per asset min
samp_count  = 1000; % number of return samples

rng(0); % seed rng generator for return sampling
ret_linma   = price_to_ma_return(asdata,train_drange,'linear'); % linear mean annual return
ret_sample  = price_to_retsample(asdata,train_drange,samp_count,'linear'); % return covariance matrix

% when asset does not contain price data for the period replace NaN
% with following values and optimization will not include the asset(s)
ret_sample(:,isnan(ret_linma)) = -1;
ret_linma(isnan(ret_linma))    = 0;

% mean-cvar optimization function with: return mean, return
% samples, and beta as inputs, and min and max weight per asset as
% optional inputs. function returns efficient frontier in format 
%   [cVaR, return, whether data point is efficient relative to other points]
[eff_frontier,eff_weights] = mean_cvar_optimization(ret_linma,ret_sample,beta,wpa_min,wpa_max);
[~,idx1]    = max(eff_frontier(:,2) ./ eff_frontier(:,1)); % select portfolio with highest mean / variance ratio
w1          = eff_weights(idx1,:);
w2          = zeros(size(w1)) + (1 / asdata.count);

% calculate optimal portfolio and benchmark portfolio metrics
drc1 = {train_drange validate_drange test_drange}; % date range cell
for i1 = 1:numel(drc1) % loop through all date range
    dr = drc1{i1};
    
    % create investment schedule matrix, invest $10,000 at start date
    inschedule = [dr(1) 1e4]; 

    pf1 = calc_portfolio_return(asdata,dr,w1,inschedule,'real',irdata);
    pf2 = calc_portfolio_return(asdata,dr,w2,inschedule,'real',irdata); 
    
    % print results
    fprintf('%d-%d PORTFOLIO A - SELECTED PORTFOLIO\n\n',year(dr));
    str1_title = sprintf(' %9s %9s %9s %9s %9s\n','YEAR','RFR','RET','STDEV','SHARPE');
    str2 = sprintf(' %9d %9.3f %9.3f %9.3f %9.3f\n',[year(pf1.metrics_yearly(:,1)) pf1.metrics_yearly(:,2:end)]');
    fprintf('%s%s\n\n',str1_title,str2);
    
    fprintf('%9s PORTFOLIO B - BENCHMARK PORTFOLIO\n\n','');
    str2 = sprintf(' %9d %9.3f %9.3f %9.3f %9.3f\n',[year(pf1.metrics_yearly(:,1)) pf1.metrics_yearly(:,2:end)]');
    fprintf('%s%s\n\n',str1_title,str2);
end 
