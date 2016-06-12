function pfdata = calc_portfolio_return(asdata,date_range,w,inschedule,varargin)
    % calculate portfolio return, value given asset data, date range, asset
    % weights, and investment schedule
    %
    % investment schedule is a two column matrix with the first column
    % as investment date and second column investment amount
    %
    % two different methods for calculating returns, theoretical which
    % allows for fractional shares, and real (less theoretical) which will
    % hold cash until next investment date until allocation accumulates over
    % ones share. real will also account for any interest return on cash
    % holdings with interest rate struct from varargin, otherwise no
    % interest return
    %
    % varargin options:
    %   method: real, theoretical (default real)
    %   interest rate struct: with timestamp, annual, daily fields
    %       annual, daily contain risk 'free' rate in format
    %           [CAD TBill,CAD Savings,US TBill,US Savings]
    %       from http://www.bankofcanada.ca/, https://research.stlouisfed.org/
    %       daily calculated as trading days not calendar days from annual
    
    % verify date range
    if ~isnumeric(date_range) || numel(date_range) < 2
        date_range = asdata.timestamp([1 end])';
    else
        date_range = [max([date_range(1) asdata.timestamp(1)]) min([date_range(2) asdata.timestamp(end)])];
    end
        
    % check if [function, time frame] options set, if multiple value for same
    % option choose higher option
    vc1         = varargin(cellfun(@(c) (ischar(c)),varargin));
    cmethod     = max([1 find(ismember({'real','theoretical'},vc1))]); 
    
    idxas       = w > 0.0001;
    
    % setup timestmap price, ret, truncate to date range
    w1          = w(idxas);
    timestamp   = asdata.timestamp;
    price       = asdata.price(:,idxas); % number of shares calculated from time price
    ret         = asdata.price_adjusted(:,idxas); % return caculcated from time price_adjusted
    ret         = [ones(1,numel(w1)); ret(2:end,:) ./ ret(1:end-1,:)];
    
    idx1        = (timestamp(:,1)>=date_range(1)) & (timestamp(:,1)<=date_range(2));
    timestamp   = timestamp(idx1,:);
    price       = price(idx1,:);
    ret         = ret(idx1,:);
    ret(~isfinite(ret)) = 1;
    
    % setup investment schedule
    inschedule(inschedule(:,1) < date_range(1),1) = date_range(1); % investment before date range, set to date range start
    [~,idx1]    = sort(inschedule(:,1));
    inschedule  = inschedule(idx1,:);
    inschedule(:,3) = 0;
    
    % precompute / find investment schedule time in time series for ret, price 
    % look up etc, rather than find each time
    idx1        = 0;
    for i1 = 1:numel(timestamp)
        while (idx1 < size(inschedule,1)) && (inschedule(idx1+1,1) <= timestamp(i1))
            idx1 = idx1+1;
            inschedule(idx1,3) = i1;
        end
        if idx1 >= size(inschedule,1)
            break;
        end
    end
    inschedule(~inschedule(:,3),:) = []; % investment after date range, removed
    if isempty(inschedule)
        fprintf('Investment schedule falls outside of date range.\n');
    end
    
    % check if interest rate / risk-free data given otherwise irdata with zero return
    irdata      = [];
    idx1        = find(cellfun(@(c) (isstruct(c)),varargin));
    for i1 = idx1(:)'
        if ismember(fieldnames(varargin{i1}),{'timestamp','annual','daily'})
            irdata = varargin{i1};
            break;
        end
    end
    if isempty(irdata)
        irdata = struct();
        irdata.timestamp = asdata.timestamp;
        irdata.annual = zeros(size(irdata.timestamp,1),4);
        irdata.daily  = irdata.annual;
    end

    % match interest / risk-free rate to data timestamp
    [~,idx1]    = ismember(timestamp,asdata.timestamp);
    irret       = irdata.daily(idx1,:) + 1; 
        
    %% calculate
    % allocate asset value and remainder cash value accounts
    assets_value    = zeros(size(timestamp,1),numel(w1));
    cs_value        = zeros(size(timestamp,1),numel(w1)); % cash remainder from share rounding down in each asset 'account'
    
    % additional investment makes it difficult to calculate returns
    % from the assets_value, cs_value above because the jump in value 
    % are mistaken as return. allocate two new matrices two record accurate
    % last values for later returns calculations
    assets_valuel   = zeros(size(timestamp,1),numel(w1)); % record assets value last to accurately calculate returns, otherwise each additional investment from the schedule appears as asset return
    cs_valuel       = zeros(size(timestamp,1),numel(w1)); % record cash 'savings account' value last to accurately calculate returns
    if cmethod == 2 
        %% simple theoretical
        for i1 = 1:size(inschedule,1)
            idx1        = inschedule(i1,3); % precomputed corresponding idx on timestamp to lookup price, ret on date
            calloc      = inschedule(i1,2) * w1; % asset value at start, cash allocated to each asset
            assets_value(idx1:end,:) = assets_value(idx1:end,:) + ...
                bsxfun(@times,calloc,cumprod(ret(idx1:end,:),1)); % asset value at each time
        end
        
    else % default
        %% realistic, accounts for whole shares only, track cash and interest rate
        for i1 = 1:size(inschedule,1)
            idx1        = inschedule(i1,3); % precomputed corresponding idx on timestamp to lookup price, ret on date            
            calloc      = inschedule(i1,2) * w1; % asset value at start, cash allocated to each asset
            calloc      = calloc + cs_value(idx1,:); % add cash accumulated in asset 'account' previously
            [assets_value,cs_value,assets_valuel,cs_valuel] = calc_real( ...
                assets_value,cs_value,assets_valuel,cs_valuel,idx1,calloc,price,ret,irret(:,2)); % calculate asset and cash 'saving account' value using CAD savings rate
            
            % look bewteen this and next period (or end date range) for
            % opprotunity to deploy cash (when price has fallen, corporate
            % action cause single share price to be below cash, or cash
            % accumulated enough interest to buy a share)
            idx2 = size(cs_value,1);
            if i1 < size(inschedule,1)
                idx2 = inschedule(i1+1,3)-1;
            end
            for i2 = 1:size(cs_value,2) % iterate through assets individually
                idx3 = idx1-1 + find( (price(idx1:idx2,i2)~=0) & ...
                    (cs_value(idx1:idx2,i2) >= price(idx1:idx2,i2)) ,1,'first'); % find idx3 where price ~= 0 and cash >= price
                while ~isempty(idx3) % there may be more than one purchase opportunity, do while any cash >= price in period
                    % only calculate for single asset
                    calloc = cs_value(idx3,i2); % cash available for asset
                    [assets_value(:,i2),cs_value(:,i2),assets_valuel(:,i2),cs_valuel(:,i2)] = calc_real( ...
                        assets_value(:,i2),cs_value(:,i2),assets_valuel(:,i2),cs_valuel(:,i2),idx3,calloc,price(:,i2),ret(:,i2),irret(:,2));
                    idx3 = idx1-1 + find( (price(idx1:idx2,i2)~=0) & ...
                        (cs_value(idx1:idx2,i2) >= price(idx1:idx2,i2)) ,1,'first'); % find any other idx3 where price ~= 0 and cash >= price
                end
            end
        end
        
    end
    
    pfdata                  = struct();
    pfdata.assets_count     = sum(idxas);
    pfdata.assets_name      = asdata.name(idxas);
    pfdata.timestamp        = timestamp;
    pfdata.assets_value     = assets_value;
    pfdata.cs_value         = cs_value;
    pfdata.portfolio_value  = sum(assets_value + cs_value,2);
    pfdata.portfolio_dlret  = sum(assets_value + cs_value,2) ./ sum(assets_valuel + cs_valuel,2) - 1; % daily linear return
    pfdata.portfolio_dlret(~isfinite(pfdata.portfolio_dlret)) = 0; % correct divided by zero values (happens when asset has not been listed)
    
    pfdata.metrics          = pfmetrics(pfdata.timestamp,pfdata.portfolio_dlret,irret(:,1),'all'); % using CAD TBill as risk-free rate
    pfdata.metrics_annual   = []; % formatting struct field order
    pfdata.metrics_yearly   = pfmetrics(pfdata.timestamp,pfdata.portfolio_dlret,irret(:,1),'yearly'); % using CAD TBill as risk-free rate
    pfdata.metrics_annual   = [pfdata.metrics(:,1) mean(pfdata.metrics_yearly(:,2:end),1)];
    pfdata.metrics_monthly  = pfmetrics(pfdata.timestamp,pfdata.portfolio_dlret,irret(:,1),'monthly'); % using CAD TBill as risk-free rate
end

function [pfmet] = pfmetrics(timestamp,dlret,irret,timeunit)
    % portfolio metrics in given timeunit   
    switch timeunit
        case 'weekly'
            idx1 = weeknum(timestamp);
        case 'monthly'
            idx1 = month(timestamp);
        case 'yearly'
            idx1 = year(timestamp);
        case 'all'
%             % compounded annual, note pmet(5) in this case is geometric shape ratio
%             tu   = numel(unique(year(timestamp)));
%             pfmet = [0 prod(irret)^(1/tu)-1 prod(dlret+1)^(1/tu)-1 std(dlret)*sqrt(252)];
%             pfmet(1,5) = (pfmet(1,3) - pfmet(1,2)) / pfmet(1,4);
            idx1 = true(size(timestamp));
    end
    
    if numel(idx1) > 1
        idx1 = find(idx1(2:end)~=idx1(1:end-1));
        idx1 = [[1; idx1+1] [idx1; numel(timestamp)]]; % idx1+1 to not double count
        idx1(idx1(:,1)==idx1(:,2),:) = [];
    else
        idx1 = logical([idx1 idx1])*1;
    end
    
    pfmet  = zeros(size(idx1,1),5); % [timestamp, return, standard-deviation, risk-free return, sharpe ratio]
    for i1 = 1:size(idx1,1)
        pfmet(i1,1) = timestamp(idx1(i1,2)); % timestamp
        pfmet(i1,2) = prod(irret(idx1(i1,1):idx1(i1,2)))-1; % risk-free return
        pfmet(i1,3) = prod(dlret(idx1(i1,1):idx1(i1,2))+1)-1; % return
        pfmet(i1,4) = std(dlret(idx1(i1,1):idx1(i1,2))) * sqrt(idx1(i1,2)-idx1(i1,1)+1); % standard deviation from daily to time unit
        pfmet(i1,5) = (pfmet(i1,3) - pfmet(i1,2)) / pfmet(i1,4);
    end
end

function [assets_value,cs_value,assets_valuel,cs_valuel] = calc_real(assets_value,cs_value,assets_valuel,cs_valuel,idx1,calloc,price,ret,irret)
    ashare                      = floor(calloc ./ price(idx1,:)); % assets share = cash / price per share round down
    ashare(~price(idx1,:))      = 0; % if price is zero, asset hasn't been listed, buy 0 shares

    % calculate extra cash after buying max shares
    cs_value(idx1,:)            = calloc - (price(idx1,:) .* ashare); % after buying max shares, put cash in asset savings 'account'
    cs_value(idx1:end,:)        = bsxfun(@times,cs_value(idx1,:),cumprod(irret(idx1:end,1),1));  % accumulate interest forward if applicable
    cs_valuel(idx1,:)           = cs_value(idx1,:);
    cs_valuel(idx1+1:end,:)     = cs_value(idx1:end-1,:);

    % recalculate actual cash used after buying max shares
    calloc                      = price(idx1,:) .* ashare;
    assets_value(idx1:end,:)    = assets_value(idx1:end,:) + ...
        bsxfun(@times,calloc,cumprod(ret(idx1:end,:),1)); % asset value at each time
    assets_valuel(idx1,:)       = assets_valuel(idx1,:) + calloc;
    assets_valuel(idx1+1:end,:) = assets_value(idx1:end-1,:);
end