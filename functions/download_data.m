function download_data(datadir_name,dataset_name)% settings

    % dataset_name    = 'djia'; % dataset to download, 'sp500': S&P500, 'djia': Dow Jones Industrial Average
    % curl_path       = '\usr\bin\curl'; % path to curl for alt download method
    
    if isempty(strfind(datadir_name,pwd)) % change relative to absolute path
        datadir_name = [pwd filesep datadir_name];
    end
    if datadir_name(end) == filesep
        datadir_name(end) = [];
    end

    % check asset csv file, setup directory
    asset_csv_path   = [datadir_name filesep dataset_name '_assets.csv']; % asset csv path
    savedir_name    = [datadir_name filesep dataset_name];
    if ~exist(asset_csv_path,'file');
        fprintf('Missing dataset assets csv file.\n');
        return;
    end
    if ~exist(savedir_name,'file');
        mkdir(savedir_name);
    end

    % read assets csv, get first column must be ticker symbol
    tickers = csv_to_cell(asset_csv_path);
    if size(tickers,1) > 1
        tickers     = tickers(:,1);
    end
    
    for i1 = 1:numel(tickers)
        savefile_name   = [savedir_name filesep strrep(tickers{i1},'.','_') '.csv'];

        if ~exist(savefile_name,'file') % don't download twice
            tic;
            fprintf('Downloading %s data from Yahoo!: ',tickers{i1});

            yahoourl    = sprintf('http://ichart.finance.yahoo.com/table.csv?s=%s',tickers{i1});
            [~,~]       = urlwrite(yahoourl,savefile_name); % [~,~] to suppress error which is checked later

            % % use websave for MATLAB R2014b and above
            % [~]         = websave(yahoourl,savefile)
            % % use curl, curl_path point to the binary
            % curlstr1    = sprintf('"%s" -k -o "%s" "%s"',curl_path,savefile_name,yahoourl);
            % [~,~]       = system(curlstr1);

            if exist(savefile_name,'file')
                file_id = fopen(savefile_name,'r');
                str1    = fread(file_id,4,'*char');
                fclose(file_id);
                if ~strcmpi(str1(:)','Date')
                    % if wrong format delete file
                    fprintf(2,'INVALID FILE ');
                    delete(savefile_name);
                end
            else
                fprintf(2,'FAILED TO DOWNLOAD ');
            end
            toc;
        end
    end
end