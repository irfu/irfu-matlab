function varargout = psp_var(vartxt)
%% PSP_VAR - obtain information about PSP data variables
%
% PSP_VAR * - display all variables
% PSP_VAR xx - display variables matching xx
% X=PSP_VAR('varName') - get variable 'varName' information 
%
% See also: PSP

if nargin == 0
  help psp_var;
  return;
end

doFindFile = any(strfind(vartxt,'file='));

%% find and read psp_variables.txt
tmp = which('psp_var');
pspFile = [tmp(1:end-5) 'variables.txt'];
fid = fopen(pspFile,'rt');
C = textscan(fid, '%s', 'Delimiter',''); C = C{1};
fclose(fid);

%% split psp_variables.txt into structure array of variable 
% line numnbers of start/end of each structure
startIdx = find(ismember(C, '%%%%'))+1;
endIdx = find(ismember(C, '----'))-1;
startIdxInfo = find(ismember(C, 'info#start'))+1;
endIdxInfo = find(ismember(C, 'info#end'))-1;

% define array of strucutres
N = numel(startIdx);
arr = struct('varName','','fileName','','hourtag',{''},'directory','','shortVar','','realted',{''},'info',"");
arr = repmat(arr,[N 1]);


% parse and store each variable in the structure array 
for i=1:N
  % parse key/value of struct
  s = C(startIdx(i):startIdxInfo(i));
  s = regexp(s,'(\w+)\s*[:=]\s*([^%$]*)(?:%[^$]*)?','tokens','once');
  s = vertcat(s{:});
  
  %  try to parse as numbers
  v = str2double(s(:,2));
  s(~isnan(v),2) = num2cell(v(~isnan(v)));
  
  % store: struct.key = value
  for j=1:size(s,1)
    arr(i).(s{j,1}) = s{j,2};
  end
  
  % fix hourtag
  if isempty(arr(i).hourtag)
    arr(i).hourtag={''};
  elseif strcmp(arr(i).hourtag,'6h')
    arr(i).hourtag={'00';'06';'12';'18'};
  end
  % read info
  arr(i).info = string(C(startIdxInfo(i):endIdxInfo(i)));
end

%% match 
doMatchNameExact = false(1,N);
doMatchName      = false(1,N);

if strcmp(vartxt,'*') % show all variables
  doMatchName(:) = true;
elseif doFindFile
  nameToMatch = vartxt(6:end); % remove "file=" from the beginning
  for i=1:N
    doMatchName(i) = strcmpi(arr(i).fileName, nameToMatch);
  end
else
  for i=1:N
    if strcmpi(arr(i).varName, vartxt) ...
        || any(regexpi(vartxt,['^' arr(i).varName '$'])) ...
        || any(regexpi(vartxt,['' arr(i).varNameShort '$']))
      doMatchNameExact(i) = true;
      iMatchExact = i; 
    end
    if any(regexpi(arr(i).varName,vartxt)) ||...
        any(regexpi(arr(i).varNameShort,vartxt))
      doMatchName(i) = true;
    end
  end
end

%% Output
if nargout == 0
  % print all matching variables
  if any(doMatchNameExact)
    fprintf('%s\n',C{startIdx(iMatchExact)-1:endIdx(iMatchExact)+1});
  else
    if numel(find(doMatchName))> 1
      for i = find(doMatchName)
        fprintf('%s\n',['<a href="matlab: psp_var ' arr(i).varName '">' arr(i).varName '</a>']);
      end
      else
      for i = find(doMatchName)
        fprintf('%s\n',C{startIdx(i)-1:endIdx(i)+1});
      end
    end
  end
elseif nargout == 4 && any(doMatchNameExact)
  %[fileName,varName,hourTag,shortVar] 
  if any(doMatchNameExact)
    out = arr(iMatchExact);
    if any(regexpi(vartxt,['^' arr(iMatchExact).varName '$'])) % varName matched regexp
      out.varName = vartxt;
    end
    if isempty(out.varNameShort)
      out.varNameShort = out.varName; 
    end
  end
  varargout= {out.fileName, out.varName,out.hourtag,out.varNameShort};
else
  if any(doMatchNameExact)
    out = arr(iMatchExact);
    if any(regexpi(vartxt,['^' arr(iMatchExact).varName '$'])) % varName matched regexp
      out.varName = vartxt;
    end
  else
    out = cell(1,sum(doMatchName == true));
    iMatch = find(doMatchName);
    if numel(iMatch) == 1
      out = arr(iMatch);
    else
      for i = 1:numel(iMatch)
        out{i} = arr(iMatch(i));
      end
    end
  end
  varargout(1)={out};
end
