function [eff_frontier,eff_weights] = mean_variance_optimization(ret_linma,ret_covm,varargin)
    %% mean-variance portfolio optimization
    
    assets_count        = size(ret_covm,2); % asset count
    min_weight_pa       = 0.00; % min weight per asset default
    max_weight_pa       = 1.00; % min weight per asset default
    if numel(varargin) >= 1 && isnumeric(varargin{1}) && ~isempty(varargin{1})
        min_weight_pa   = varargin{1};
    end
    if numel(varargin) >= 2 && isnumeric(varargin{2}) && ~isempty(varargin{2})
        max_weight_pa   = varargin{2};
    end

    % when asset does not contain price data for the period replace NaN
    % with following values and optimization will not include the asset(s)
    ret_linma(isnan(ret_linma)) = 0;
    ret_covm(isnan(ret_covm))   = 1;

    % eff_front       = (0.01:0.005:1.00)'; % test larger range for long and short
    eff_frontier    = [min(ret_linma):0.01:max(ret_linma) max(ret_linma)]';
    eff_frontier    = [zeros(size(eff_frontier,1),1) eff_frontier zeros(size(eff_frontier,1),1)];
    eff_weights     = zeros(size(eff_frontier,1),assets_count);

    % mean-variance Markowitz, solve using quadratic programming
    % quadprog setup, see doc quadprog
    qpo = optimset('Algorithm','interior-point-convex','Display','Off');
    H   = ret_covm;
    f   = zeros(assets_count,1);
    A   = -eye(assets_count);
    b   = zeros(assets_count,1)+1; % +1 to allow long and short, change to +0 to enforce long only
    Aeq = [ones(1,assets_count); ret_linma];
    beq = [1.00; 0.00];
    lb  = zeros(assets_count,1);
    ub  = zeros(assets_count,1);
    ub(1:assets_count) = max_weight_pa;
    lb(1:assets_count) = min_weight_pa;

    fprintf('Minimizing variance at given returns:\n\n');
    for i1 = 1:size(eff_frontier,1)
        fprintf('%10s',sprintf('%6.2f%%',eff_frontier(i1,2)*100));
        if ~mod(i1,8)
            fprintf('\n');
        end    
        
        beq(2) = eff_frontier(i1,2);
        [w,pvar,status] = quadprog(H,f,A,b,Aeq,beq,lb,ub,[],qpo); % weights, portfolio variance, status
        if status ~= 1
            eff_frontier(i1,1:2) = nan;
        else
            % portfolio variance to portfolio standard deviation
            % pvar * 2 since quadprog optimizes 1/2*x'*H*x, see doc quadprog
            eff_frontier(i1,1) = sqrt(pvar * 2);
            eff_weights(i1,:) = w;
        end
    end
    fprintf(repmat('\n',1,(mod(i1,8)~=0)+1));

    eff_weights(abs(eff_weights) < 0.0001) = 0;
    for i1 = 1:size(eff_weights,1)
    end
    if ~sum(isfinite(eff_frontier(:,1)))
        fprintf('No solution to efficient frontier, check weight limits, date range, underlying data\n');
    end

    eff_frontier(end,3) = isfinite(eff_frontier(end,1));
    for i1 = 1:size(eff_frontier,1)-1
        if isfinite(eff_frontier(i1,1))
            eff_frontier(i1,3) = ~sum((eff_frontier(i1,1) > eff_frontier(i1+1:end,1)));
        end
    end
end