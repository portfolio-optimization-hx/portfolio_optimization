function covmatrix = price_to_covmatrix(asdata,date_range,varargin)
    % calculate covariance matrix return from daily price_adjusted
    % varargin option for linear, log returns
    
    % verify date range
    if ~isnumeric(date_range) || numel(date_range) < 2
        date_range = asdata.timestamp([1 end])';
    else
        date_range = [max([date_range(1) asdata.timestamp(1)]) min([date_range(2) asdata.timestamp(end)])];
    end
    
    varargin  = varargin(cellfun(@(c) (ischar(c)),varargin));
    if sum(ismember('log',varargin))
        ret = log(asdata.price_adjusted(2:end,:) ./ prfdata.price_adjusted(1:end-1,:));
    else
        ret = asdata.price_adjusted(2:end,:) ./ asdata.price_adjusted(1:end-1,:) - 1;
    end
    ret = ret( (asdata.timestamp(2:end) >= date_range(1)) & (asdata.timestamp(2:end) <= date_range(2)) ,:);
    
    covmatrix = zeros(asdata.count,asdata.count);
    % use max number of data available, otherwise change nan,inf in ret o zero and use cov(ret)
    for i1 = 1:asdata.count
        for i2 = 1:asdata.count
            if i1 < i2
                ri1 = find(isfinite(ret(:,i1)) & isfinite(ret(:,i2)),1,'first'); % first row where both equity starts trading
                cmp = cov(ret(ri1:end,[i1 i2]));
                covmatrix(i1,i1) = cmp(1,1);
                covmatrix(i1,i2) = cmp(1,2);
                covmatrix(i2,i1) = cmp(1,2);
                covmatrix(i2,i2) = cmp(2,2);
            end
        end
    end
    
    covmatrix = covmatrix * 252; % daily covariance to annual covariance 252 trading days per year
end