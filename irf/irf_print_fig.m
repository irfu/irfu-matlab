function [out] = irf_print_fig(x1,x2,x3)
%IRF_PRINT_FIG Exports figure to file.
%   IRF_PRINT_FIG(fileName) exports current figure to a .eps with name
%   fileName or file type specified in name, e.g: 'lorem.png'.
%   IRF_PRINT_FIG(f,fileName) exports figure with handle f to a .eps named
%   fileName. If f is an integer, the figure with that number is printed.
%   IRF_PRINT_FIG(fileName,fileType) exports current figure to a file named
%   fileName with format fileType.
%   IRF_PRINT_FIG(f,fileName,fileType) exports figure with handle f
%   to a file named fileName with format fileType.
%
%   File types:
%       'eps' - vector graphics
%       'pdf' - vector graphics, white space will be included
%       'png' - lossless compression
%
%   TODO: 
%       - Remove white space for pdf


%% Input
figNum = 0;
handleInput = 0;
posFileTypes = {'eps','pdf','png'};

if(nargin == 1)         % only name
  [filePath, fileName, fileType] = fileparts(x1);
  if ~isempty(fileType), fileType=fileType(1:end-1); else, fileType='eps'; end
  fileName = [filePath, filesep, fileName];
elseif(nargin == 2)
    if(ischar(x1))      % name and type
        fileName = x1;
        fileType = x2;
    elseif(isnumeric(x1) && floor(x1)==x1) % number and name
        figNum = x1;
        fileName = x2;
    elseif(~isnumeric(x1) && isvalid(x1)) % handle and name
        handleInput = 1;
        fileName = x2;
    else
        irf.log('critical','unknown input type')
        return
    end
elseif(nargin == 3)     % number/handle, name and type 
    if(isnumeric(x1) && floor(x1)==x1) %number
        figNum = x1;
    elseif(~isnumeric(x1) && isvalid(x1)) % handle
        handleInput = 1;
    else
        irf.log('critical','unknown input type')
        return
    end
    fileName = x2;
    fileType = x3;
end

if(strcmp(fileName,'')) % Do not print if name is empty
    irf.log('w','No file name given, will not print')
    return;
end
if ~ismember(fileType,posFileTypes)
    irf.log('w',['Not supported file type: ',fileType])
    return;
end


%% File flags 
switch fileType
    case 'eps'
        fileFlag = '-depsc';
    case 'pdf'
        fileFlag = '-dpdf';
    case 'png'
        fileFlag = '-dpng';
    otherwise
        irf.log('critical','unknown file type')
        return
end

%% Set figure handle
if(handleInput)
    switch x1.Type
        case 'figure'
            irf.log('warning','figure input')
            f = x1;
        case 'axes'
            irf.log('warning','axes input, printing parent')
            f = x1(1).Parent;
        otherwise
            irf.log('critical','unknown input type')
            return
    end
elseif(figNum == 0) %never assigned
    irf.log('warning','printing current figure')
    f = gcf;
    figNum = f.Number;
else
    f = figure(figNum);
end
irf.log('warning',['printing figure ',num2str(figNum)])
irf.log('warning',['saving as: ',fileName,'.', fileType])

%% Exporting figure
fileStr = [fileName,'.',fileType];

if(nargout == 1)
    out = f;
end

% Exporting the figure
print(f,fileFlag,'-loose','-painters',fileStr) 

end