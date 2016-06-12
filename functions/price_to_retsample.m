function retsample = price_to_retsample(asdata,date_range,n,varargin)
    % create n sample annual return for each asset from price_adjusted
    % use price_adjusted to account for corporate action
    % mean annual return from annualized from base time frame, because
    % simply doing yearly can have incomplete data
    % varargin options:
    %   function:   linear, log (default linear)
    %   sample method: monte carlo from historical daily (default monte carlo)
    %                  actual annual historical 
    % historical sampling methods should only be used when there is enough
    % sample data
        
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
    smethod     = max([1 find(ismember({'mc','hist'},varargin))]);
    
    
    % allocate return sample
    if smethod == 2
        retsample   = (year(date_range(1)):year(date_range(2)))';
        retsample   = [retsample zeros(numel(retsample),asdata.count)];
    else
        retsample   = zeros(n,asdata.count);
    end
    for i1 = 1:asdata.count
        timepricea = [asdata.timestamp asdata.price_adjusted(:,i1)];
        if numel(date_range) >= 2
            timepricea = timepricea( (timepricea(:,1)>=date_range(1)) & (timepricea(:,1)<=date_range(2)) ,:);
        end
        timepricea(~timepricea(:,2),:) = []; % remove rows before equity started trading
        
        if smethod == 2
            timepricea(:,1) = year(timepricea(:,1));
            idx1 = find(timepricea(2:end,1) ~= timepricea(1:end-1,1)) + 1;
            idx1 = [[1; idx1] [idx1; size(timepricea,1)]];
            idx1((idx1(:,2)-idx1(:,1)+1)<248,:) = []; % remove incomplete years
            ret  = [timepricea(idx1(:,1),1) timepricea(idx1(:,1),2) timepricea(idx1(:,2),2)];
            if lin_log == 2
                ret = [ret(:,1) log(ret(:,3) ./ ret(:,2))];
            else
                ret = [ret(:,1) -1+(ret(:,3) ./ ret(:,2))];
            end
            [~,idx1] = ismember(ret(:,1),retsample(:,1));
            retsample(idx1,1+i1) = ret(:,2);
        else
            % generate monte carlo only if there is enough daily return            
            if size(timepricea,1) > 252 * 1.5 % 252 trading days per year x 1.5 years
                % uncomment line 59-62 for random index that maintains 
                % historical average win/loss ratio
                oc      = timepricea(2:end,2) > timepricea(1:end-1,2);
                oc      = {find(~oc) find(oc)};
                lcount  = round((size(oc{1},1) / (size(timepricea,1)-1)) * 252); % loss count
                idxrn   = [oc{1}(randi(numel(oc{1}),n,lcount)) oc{2}(randi(numel(oc{2}),n,252-lcount))]; % generate random index for n samples x 252 days
                
                % uncomment line 61 to use alternative random index that will 
                % produce wider distribution, further tails, larger stdev
                idxrn   = [randi(size(timepricea,1)-1,n,252)]; % generate random index for n samples x 252 days
                if lin_log == 2
                    ret = log(timepricea(2:end,2) ./ timepricea(1:end-1,2)); % get daily log return
                    ret = ret(idxrn); % create random return based on random index
                    retsample(:,i1) = mean(ret,2).*252;
                    % retsample(:,i1) = sum(ret,2); % annual log return = sum(252 daily log returns)
                else
                    ret = (timepricea(2:end,2) ./ timepricea(1:end-1,2)); % get daily linear return
                    ret = ret(idxrn); % create random return matrix based on random index
                    retsample(:,i1) = mean(ret,2).^252 - 1;
                    % retsample(:,i1) = prod(ret,2) - 1; % annual linear return = prod(252 daily log returns) - 1
                end
            end
        end
    end
    
    if smethod == 2
        % for historical, if requested n sample is less than available
        % historical samples, randomly select n samples
        if n < size(retsample,1)
            retsample = retsample(randperm(size(retsample,1),n),:);
        end
        retsample(:,1) = []; % remove year column
    end
    retsample = sort(retsample,1);
end