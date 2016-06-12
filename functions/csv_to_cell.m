function csv_cell = csv_to_cell(file_path,varargin)
    % read csv file and convert to cell
    % take delimiter string as optional argument
    % read line by line to handle larger files 
    % data stored as string in cell, does not perform conversion
    
    csv_cell = {};
    if ~exist(file_path,'file')
        fprintf('%s does not exist.\n',file_path);
        return;
    end
    
    % default delimiter ',' accept optional delimiter from varargin
    delimiter_str = ',';
    if numel(varargin) && ischar(varargin{1})
        delimiter_str = varargin{1};
    end
    
    % allocate cell, approximate size
    file_id   = fopen(file_path,'r');
    fline     = fgetl(file_id);
    
    file_info = dir(file_path);
    
    csv_cell    = cell( ...
        ceil(file_info.bytes/numel(fline)*1.5), ...
        numel(strfind(fline,delimiter_str))+1 ...
        );
    csv_cell_ri = 0;
    
    while ischar(fline)
        % format line add to csv_cell
        csv_cell_ri = csv_cell_ri + 1;
        tcell   = textscan(fline,'%s','Delimiter',delimiter_str);
        csv_cell(csv_cell_ri,1:numel(tcell{:})) = tcell{:};
        
        % get next line
        fline   = fgetl(file_id);
    end
    csv_cell = csv_cell(1:csv_cell_ri,:);
    
    fclose(file_id);
end