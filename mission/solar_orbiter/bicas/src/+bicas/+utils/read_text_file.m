function rowStrList = read_text_file(path)
% Read text file into a cell array of strings, one string per row.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-10-12


    fileId = fopen(path, 'r');    
    
    % ASSERTION (Gives better error message than textscan failure.)
    if fileId == -1
        error('BICAS:read_text_file:PathNotFound', 'Can not find file "%s".', path)
    end    
    
    temp = textscan(fileId, '%s', 'Delimiter', '\n', 'Whitespace', '');
    rowStrList = temp{1};
    fclose(fileId);
end

