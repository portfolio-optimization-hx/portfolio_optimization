function asdata = ascsv_to_asdata(csv_fpaths)
    asdata = [];
    
    pa_ndays        = 255*2; % preallocate 255 trading days * 20 years
    timestamp       = (today-pa_ndays+1:today)'; % ascending order
    price           = zeros(pa_ndays,numel(csv_fpaths)); % close price 
    price_adjusted  = zeros(pa_ndays,numel(csv_fpaths)); % close price adjusted for corporate actions
    as_name         = cell(1,50); % portfolio asset names
    as_count        = 0; % numer of portfolio assets
    
    fprintf('\nProcessing %d files ...\n',numel(csv_fpaths));
    fprintf(['\n    (ensure csv dates are properly formatted and turn off' ...
             '\n     date format check [line 25-26] for better performance) \n\n']);
    for i1 = 1:numel(csv_fpaths)
        % read csv files only
        if (numel(csv_fpaths{i1}) > 4) && (strcmp(csv_fpaths{i1}(end-3:end),'.csv')) && (exist(csv_fpaths{i1},'file'))
            celldata1 = csv_to_cell(csv_fpaths{i1});
            
            % if cell data in proper format with at least 1 row, 7 columns
            if (size(celldata1,1) >= 1) && (size(celldata1,2) >= 7)
                % logical index to check string contain only valid date chars,
                % remove invalid rows. this check has high performance
                % penalty, can be removed if csv files guaranteed to be in
                % correct format. yahoo files can get away with just line 28
%                 li1 = cellfun(@(c) (sum(~ismember(c,['0' '1' '2' '3' '4' '5' '6' '7' '8' '9' '-']))),celldata1(:,1)) ~= 0; 
%                 celldata1(li1,:) = [];
                
                celldata1(1,:) = []; % for yahoo csv and comment out line 25-26

                % convert cell to temp matrix format 
                % [timestamp close-price close-price-adjusted]
                m1  = [datenum(celldata1(:,1),'yyyy-mm-dd') str2double(celldata1(:,5))  str2double(celldata1(:,7))];

                % check if date lies outside of existing timestamp range
                date_min = min(m1(:,1));
                if date_min < timestamp(1)
                    price            = [
                        zeros(timestamp(1)-date_min,size(price,2))
                        price
                        ];
                    price_adjusted   = [
                        zeros(timestamp(1)-date_min,size(price_adjusted,2))
                        price_adjusted
                        ];
                    timestamp        = [
                        (date_min:timestamp(1)-1)'
                        timestamp
                        ];
                end
                date_max = min(m1(:,1));
                if date_max > timestamp(end)
                    price            = [
                        price
                        zeros(date_max-timestamp(end),size(price,2))
                        ];
                    price_adjusted   = [
                        price_adjusted
                        zeros(date_max-timestamp(end),size(price_adjusted,2))
                        ];
                    timestamp        = [
                        timestamp
                        (timestamp(end)+1:date_max)'
                        ];
                end
                
                % match dates add to asdata
                [li1,idx1] = ismember(m1(:,1),timestamp);
                if sum(~li1)
                    % should not occur already handled in above code
                    fprintf('Date(s) missing in asdata, check ascsv_to_asdata.\n');
                end
                as_count = as_count+1;
                price(idx1,as_count)          = m1(li1,2);
                price_adjusted(idx1,as_count) = m1(li1,3);
                as_name{as_count}          = csv_fpaths{i1}(max([1 find(csv_fpaths{i1}==filesep,1,'last')+1]):end-4);
                
                fprintf('%10s',as_name{as_count});
                if ~mod(as_count,8)
                    fprintf('\n');
                end
            end
        end
    end
    fprintf(repmat('\n',1,(mod(i1,8)~=0)+1));
    fprintf('\nDone - processed %d csv files.\n',as_count);
    
    if as_count
        % remove empty rows, no price should be zero
        li1             = (sum(price == 0,2) < size(price,2)) | (sum(price_adjusted == 0,2) < size(price_adjusted,2));
        timestamp       = timestamp(li1,:);
        price           = price(li1,1:as_count);
        price_adjusted  = price_adjusted(li1,1:as_count);
        as_name       = as_name(:,1:as_count);
        
        % replace zeros with previous day price if previous day price is
        % valid / non-zero. this occurs on yahoo for Canadian equities
        % because yahoo doesn't recognize Canadian holidays
        price           = zerorep(price);
        price_adjusted  = zerorep(price_adjusted);
        
        % create return struct asdata
        asdata                  = struct();
        asdata.count            = as_count;
        asdata.name             = as_name;
        asdata.timestamp        = timestamp;
        asdata.price            = price;
        asdata.price_adjusted   = price_adjusted;
    end
end

function m1 = zerorep(m1)
    % replace zero with previous value if previous value is non-zero    
    for ci1 = 1:size(m1,2)
        for ri1 = 2:size(m1,1)
            if m1(ri1-1,ci1) && ~m1(ri1,ci1)
                m1(ri1,ci1) = m1(ri1-1,ci1);
            end
        end
    end
end