%
% Non-automatic test of write_CDF_dataobj. Requires manual inspection.
% Far from complete but useful for testing special oddities associated with writing CDF files.
%
%
% NOTE: Writes temporary file to disk.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2018-03-09
%
function write_CDF_dataobj___MTEST()
% NOTE: dataobj requires at least one "Epoch-like" zVar. otherwise it returns an empty dataobj.
%
% TODO: Create CDF with correct record variance.
%
% TESTS:
%   zVars dimensions (incl. #records)
%       different dimensions -- Test permutations 
%       different numbers of dimensions (incl. higher than three per record?)
%       singleton dimensions
%           1 record
%           singleton dimension as first, last, or intermediate dimension (index).
%       zero size dimensions? zero records
%   NaN/Pad values?!
%   Strings, string lists
%   Error on last dimension (per record) size=one.
%
% IMPLEMENTATION NOTE:
%   (1) dataobj contains a bug that prevents it from reading files with exactly one zVar.
%   (2) dataobj requires one zVar with a time-like variable.


%tempCdfFile = [tempname(), 'write_CDF_dataobj___TEST.bicas.cdf'];
tempCdfFile = '/home/erjo/temp/write_CDF_dataobj___TEST.bicas.cdf';    % DEBUG

dop.GlobalAttributes   = struct('A_global_attribute_name', 'A global attribute value');
dop.data               = struct();
dop.Variables          = {};
dop.VariableAttributes = struct();


dop = add_zVar(dop, 'Epoch',           int64([0; 10]),    [1], 'tt2000', -1000);
%dop = add_zVar(dop, 'ZVAR_NAME_1234',  ones(1,2,3,4), 'double', -1000);
%dop = add_zVar(dop, 'ZVAR_NAME_2341',  ones(2,3,4,1), 'double', -4000);   % NOTE: Last dimension vanishes before call to spdfcdfwrite.

%dop = add_zVar_size_test(dop, 'DIM', 0, [2]);
%dop = add_zVar_size_test(dop, 'DIM', 0, [2 3]);

%dop = add_zVar_size_test(dop, 'DIM', 0, []);   % Interpretation: One scalar per record (1x1x...)
%dop = add_zVar_size_test(dop, 'DIM', 1, []);   % Interpretation: One scalar per record (1x1x...)

%dop = add_zVar_size_test(dop, 'DIM', 0, [1]);
%dop = add_zVar_size_test(dop, 'DIM', 1, [1]);
%dop = add_zVar_size_test(dop, 'DIM', 2, [1]);
%dop = add_zVar_size_test(dop, 'DIM', 10, [1]);
%dop = add_zVar_size_test(dop, 'DIM', 2, [1 1]);
%dop = add_zVar_size_test(dop, 'DIM', 2, [1 1 1]);
%dop = add_zVar_size_test(dop, 'DIM', 2, [1 1 1 1]);
%dop = add_zVar_size_test(dop, 'DIM', 2, [1 1 1 2]);


%dop = add_zVar(dop, 'STR_0', [''],             [3, 1], 'char', '   ');
%dop = add_zVar(dop, 'STR_1a', permute(['abc'], [2,1,3]),     [3, 1], 'char', '   ');           % 1 record, 1 string per record
%dop = add_zVar(dop, 'STR_1b', permute(['1abc'; '2abc'], [2,1,3]),    [3, 1], 'char', '   ');   % 2 records, 1 string per record
%dop = add_zVar(dop, 'STR_1c', permute(['1abc'; '2abc'], [2,3,1]),    [3, 1], 'char', '   ');   % 1 record,  2 strings per record
dop = add_zVar(dop, 'STR_1a', ['abc'],             [1, 1], 'char', '   ');   % 1 record, 1 string per record
dop = add_zVar(dop, 'STR_1b', ['1abc'; '2abc'],    [2, 1], 'char', '   ');   % 1 records, 2 string per record
dop = add_zVar(dop, 'STR_1c', ['1abc'; '2abc'],    [1, 1], 'char', '   ');   % 2 records, 1 strings per record
% zValue = '';
% zValue(:,:,1) = ['1___';'2___';'3___'];
% zValue(:,:,2) = ['4___';'5___';'6___'];
% dop = add_zVar(dop, 'STR_3', zValue, [3, 1], 'char', '   ')


%dop = add_zVar(dop, 'STR_1a',  {},                  [3, 1], 'char', '   ')
%dop = add_zVar(dop, 'STR_1b',  {'abc'},             [3, 1], 'char', '   ')
%dop = add_zVar(dop, 'STR_1c',  {'abc','def'},       [3, 1], 'char', '   ')
%dop = add_zVar(dop, 'STR_2',   {'abc'; 'def'},      [3, 1], 'char', '   ')
%dop = add_zVar(dop, 'STR_3',   {'1abc', '2def'; '3ghi', '4jik'; '5lmn', '6opq'},  [3, 1], 'char', '   ')
%dop = add_zVar(dop, 'STR_4',   {'1abc', '2def'; '3ghi', '4jik'; '5lmn', '6opq'}', [3, 1], 'char', '   ')
%dop = add_zVar(dop, 'STR_5',   {'11', '12', '13'; '21', '22', '23'; '31', '32', '33'}, [3, 1], 'char', '   ')
%dop = add_zVar(dop, 'STR_6',   permute({'11', '12', '13', '14'; '21', '22', '23', '24'; '31', '32', '33', '34'}, [2,1,3]), [3, 1], 'char', '   ')

delete(tempCdfFile);    % Only warning (not error) if file does not exist.
options = {};
bicas.utils.write_CDF_dataobj(tempCdfFile, ...
    dop.GlobalAttributes, ...
    dop.data, ...
    dop.VariableAttributes, ...
    dop.Variables, options{:});


Do = dataobj(tempCdfFile);

%Do.data.STR_0
%Do.data.STR_0.data
Do.data.STR_1a
Do.data.STR_1a.data
Do.data.STR_1b
Do.data.STR_1b.data
Do.data.STR_1c
Do.data.STR_1c.data
%Do.data.STR_2
%Do.data.STR_2.data
% Do.data.STR_3
% Do.data.STR_3.data
% Do.data.STR_4
% Do.data.STR_4.data
% Do.data.STR_5
% Do.data.STR_5.data
% Do.data.STR_6
% Do.data.STR_6.data

%delete(tempCdfFile);
end





% Add zVar for the purpose of testing zVariables with different dimensions (both #records and non-record dimensions).
% NOTE: Only data type "double".
function dop = add_zVar_size_test(dop, namePrefix, nRecords, statedSizePerRecord)

% ASSERTION
if numel(nRecords) ~= 1
    error('Illegal argument nRecords.')
end


valueSize = [nRecords, statedSizePerRecord, 1];   % NOTE: Always at least size two.
nElements = prod(valueSize);
if nElements == 0
    value = zeros(valueSize);   % Requires valueSize length>=2 (which it has).
    %value = [];
else
    value = reshape(1:nElements, valueSize);
    %value = ones(valueSize);
end

name = [namePrefix, sprintf('_%i_', nRecords), sprintf('%i', statedSizePerRecord)];    % NOTE: Second sprintf string will be repeated since sizePerRecord is a vector.

dop = add_zVar(dop, name, value, statedSizePerRecord, 'double', -1000);
end



function dop = add_zVar(dop, name, value, statedSizePerRecord, dataType, padValue)
% NOTE: "dataType" can (maybe) in principle be read from "value".

    %temp          = size(value);
    %sizePerRecord = temp(2:end);

    dop.data.(name).data   = value;
    dop.Variables(end+1,:) = {name, statedSizePerRecord, [], dataType, [], [], [], [], padValue, [], [], []};
end



% Shift a 1D vector to be a 1D vector in dimension k.
% function v = vec_1D(v, k)
%     v = v(:);    % Force column vector.
%     
%     % NOTE: permute requires permutations with length>=2 ==> 1:(k+2)
%     % NOTE: circshift requires column vector.
%     v = permute(v, circshift([1:(k+2)]', k));    % NOTE: k can be zero.
% end
