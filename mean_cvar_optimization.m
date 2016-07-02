function [eff_frontier,eff_weights] = mean_cvar_optimization(ret_linma,ret_sample,beta,varargin)
    %% mean-cvar portfolio optimization
    % beta                    = 0.99;
    
    min_weight_per_asset    = 0.00; % default
    max_weight_per_asset    = 1.00; % default
    if numel(varargin) >= 1 && isnumeric(varargin{1}) && ~isempty(varargin{1})
        min_weight_per_asset = varargin{1};        
    end
    if numel(varargin) >= 2 && isnumeric(varargin{2}) && ~isempty(varargin{2})
        max_weight_per_asset = varargin{2};
    end
    
    assets_count = numel(ret_linma);
    samp_count   = size(ret_sample,1);


    % eff_frontier     = (0.01:0.005:1.00)'; % test larger range for long and short
    eff_frontier    = [min(ret_linma):0.01:max(ret_linma) max(ret_linma)]';
    eff_frontier    = eff_frontier((eff_frontier > 0) & (eff_frontier <= 2)); % limit return range
    eff_frontier    = [zeros(size(eff_frontier,1),1) eff_frontier zeros(size(eff_frontier,1),1) ];
    eff_weights     = zeros(size(eff_frontier,1),assets_count);

    % mean-CVaR, solve using linear programming described by R. T. Rockafellar and S. Uryasev
    % linprog setup, see doc linprog
    qpo = optimset('Display','Off');

    f = [zeros(1,assets_count) 1 zeros(1,samp_count) + ((1/(1-beta)) * (1/samp_count))];
    
    A = zeros(samp_count,1+assets_count+samp_count);
    A(1:end,1:assets_count) = -ret_sample;
    A(1:end,assets_count+1) = -1;
    A(1:end,assets_count+2:end) = -eye(samp_count);

    b = zeros(samp_count,1);

    Aeq = zeros(2,1+assets_count+samp_count);
    Aeq(1,1:assets_count) =  1.00;
    Aeq(2,1:assets_count) = -ret_linma';

    beq = [1; 0];

    lb = zeros(1+assets_count+samp_count,1);
    ub = zeros(1+assets_count+samp_count,1);
    lb(1:assets_count)     = min_weight_per_asset;
    ub(1:assets_count)     = max_weight_per_asset;
    ub(assets_count+1)     = inf;
    ub(assets_count+2:end) = inf;

    fprintf('Minimizing cVaR at given returns:\n\n');
    for i1 = 1:size(eff_frontier,1)
        fprintf('%10s',sprintf('%6.2f%%',eff_frontier(i1,2)*100));
        if ~mod(i1,8)
            fprintf('\n');
        end    
        
        beq(2) = -eff_frontier(i1,2);
        [w,cvar,status] = linprog(f,A,b,Aeq,beq,lb,ub,[],qpo); % weights, portfolio variance, status
        if status ~= 1
            eff_frontier(i1,1:2) = nan;
        else
            % portfolio variance to portfolio standard deviation
            % pvar * 2 since quadprog optimizes 1/2*x'*H*x, see doc quadprog
            eff_frontier(i1,1) = cvar;
            eff_weights(i1,:) = w(1:assets_count);
        end
    end
    fprintf(repmat('\n',1,(mod(i1,8)~=0)+1));

    eff_weights(abs(eff_weights) < 0.0001) = 0;
    if ~sum(isfinite(eff_frontier(:,1)))
        fprintf('No solution to efficient frontier, check weight limits, date range, underlying data\n');
    end

    for i1 = 1:size(eff_frontier,1)-1
        if isfinite(eff_frontier(i1,1))
            eff_frontier(i1,3) = ~sum((eff_frontier(i1,1) > eff_frontier(i1+1:end,1)));
        end
    end
end
