function mean_average_return = price_to_ma_return(asdata,date_range,varargin)
    % calculate mean annual return from price_adjusted
    % use price_adjusted to account for corporate action
    % mean annual return from annualized from base time frame, because
    % simply doing yearly can have incomplete data
    % varargin options:
    %   function:   linear, log (default linear)
    %   time frame: daily, weekly, monthly, yearly (default daily)
        
    % verify date range
    if ~isnumeric(date_range) || numel(date_range) < 2
        date_range = asdata.timestamp([1 end])';
    else
        date_range = [max([date_range(1) asdata.timestamp(1)]) min([date_range(2) asdata.timestamp(end)])];
    end
    
    % check if [function, time frame] options set, if multiple value for same
    % option choose higher option
    varargin    = varargin(cellfun(@(c) (ischar(c)),varargin));
    lin_log     = max([1 find(ismember({'linear','log'},varargin))]); 
    base_time   = max([1 find(ismember({'daily','weekly','monthly','yearly'},varargin))]);    
    
    mean_average_return = zeros(1,asdata.count);
    for i1 = 1:asdata.count
        timepricea = [asdata.timestamp asdata.price_adjusted(:,i1)];
        timepricea = timepricea( (timepricea(:,1)>=date_range(1)) & (timepricea(:,1)<=date_range(2)) ,:);
        timepricea(~timepricea(:,2),:) = []; % remove rows before equity started trading
        
        % find start-end price for given time frame
        % idx:  index of when timeframe starts
        % mdpp: minimum days per period, to remove non-standard data anomalies
        % aopy: average occurence per year, for annualization
        switch base_time
            case 1 % daily                
                idx1 = timepricea(:,1);
                mdpp = 1;
                aopy = 252; % average 252 trading days per year
            case 2 % weekly
                idx1 = weeknum(timepricea(:,1));
                mdpp = 3;
                aopy = 52; % 52 weeks per year
            case 3 % monthly
                idx1 = month(timepricea(:,1));
                mdpp = 19;
                aopy = 12; % 12 month per year
            case 4 % yearly
                idx1 = year(timepricea(:,1));
                mdpp = 248;
                aopy = 1; % 1 year per year
        end
        
        idx1 = find(idx1(2:end) ~= idx1(1:end-1));
        idx1 = [[1; idx1] [idx1; size(timepricea,1)]];
        idx1((idx1(:,2)-idx1(:,1)+1)<mdpp,:) = []; % remove period where number of days less than minimum required
        if lin_log == 2
            % log return
            ret = mean(log(timepricea(idx1(:,2),2) ./ timepricea(idx1(:,1),2)));
            ret = ret * aopy;
        else
            % linear return
            ret = mean(timepricea(idx1(:,2),2) ./ timepricea(idx1(:,1),2));
            ret = ret ^ aopy - 1;
        end
        mean_average_return(i1) = ret;
    end
end