function rowStrList = read_text_file(path)
% Read text file into a cell array of strings, one string per row.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-10-12

    fileId = fopen(path, 'r');    
    temp = textscan(fileId, '%s', 'Delimiter', '\n', 'Whitespace', '');
    rowStrList = temp{1};
    fclose(fileId);
end

