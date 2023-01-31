%
% Wrapper around solo.shk.read_RSH_file() to read many files at
% once.
%
%
% ARGUMENTS
% =========
% xmlFilePathsCa
%       Cell array of paths to SHK files.
% namesCa
%       See "read_RSH_file".
%
%
% RETURN VALUES
% =============
% Struct with data covering multiple files.
% See "read_RSH_file".
%
%
% NOTE: Test code with solo.shk.read_RSH_file___UTEST.
% NOTE: Does not sort timestamps.
%
%
% Author: Erik P G Johansson, IRF Uppsala, Sweden
% First created 2021-09-07.
%
function D = read_RSH_file_many(xmlFilePathsCa, namesCa)

    % NOTE: Hackish/ugly: Constant must be the same as in
    % solo.shk.read_RSH_file().
    XML_ELEMENT_TAG_CA = {'Name', 'EngineeringValue', 'TimeStampAsciiA'};

    D = struct();
    for fnCa = XML_ELEMENT_TAG_CA
        D.(fnCa{1}) = cell(0, 1);
    end

    for i = 1:numel(xmlFilePathsCa)
        %tic
        xmlFilePath = xmlFilePathsCa{i};

        fprintf('Reading "%s"\n', xmlFilePath)
        D1 = solo.shk.read_RSH_file(xmlFilePath, namesCa);

        for iFn = 1:numel(XML_ELEMENT_TAG_CA)
            fn = XML_ELEMENT_TAG_CA{iFn};
            D.(fn) = [D.(fn); D1.(fn)];
        end
        %toc
    end

end
