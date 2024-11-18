function varargout = psp_var(vartxt)
%% PSP_VAR - obtain information about PSP data variables
%
% PSP_VAR * - display all variables
% PSP_VAR xx - display variables matching xx
% X=PSP_VAR('varName') - get variable 'varName' information into structure
%
% See also: PSP
persistent C arr startIdx endIdx N

%% find and read psp_variables.txt
if isempty(C)
  C = importdata(...
    [fileparts(which('psp_var')) filesep 'psp_variables.txt']...
    );
  parse_psp_variables_file;
end

%%
if nargin == 0
  help psp_var;
  return;
end

%%
doFindFile = any(strfind(vartxt,'file='));

%% match
doMatchNameExact  = false(1,N);
doMatchNameExpand = false(1,N);
doMatchNameShort  = false(1,N);
doMatchName       = false(1,N);

if strcmp(vartxt,'*') % show all variables
  doMatchName(:) = true;
elseif doFindFile
  nameToMatch = vartxt(6:end); % remove "file=" from the beginning
  for i=1:N
    doMatchName(i) = any(regexpi(nameToMatch,['^' arr(i).fileName '$']));
  end
else
  for i=1:N
    if strcmpi(arr(i).varName, vartxt)
      doMatchNameExact(i) = true;
      iMatchExact = i;
    elseif any(regexpi(vartxt,['^' arr(i).varName '$']))
      doMatchNameExact(i) = true;
      doMatchNameExpand(i) = true;
      iMatchExact = i;
    elseif any(regexpi(vartxt,['' arr(i).varNameShort '$']))
      doMatchNameExact(i) = true;
      doMatchNameExpand(i) = true;
      doMatchNameShort(i) = true;
      iMatchExact = i;
    elseif any(regexpi(arr(i).varName,vartxt)) ||...
        any(regexpi(arr(i).varNameShort,vartxt))
      doMatchName(i) = true;
    end
  end
end

%% Output
if nargout == 0  % print all matching variables
  display_variable_description(doMatchNameExact);
  if ~any(doMatchNameExact)
    display_variable_description(doMatchName);
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
    if doMatchNameExpand(iMatchExact)
      if doMatchNameShort(iMatchExact)
        varShortList = list_variables(arr(iMatchExact).varNameShort);
        varList      = list_variables(arr(iMatchExact).varName);
        fileList     = list_variables(arr(iMatchExact).fileName);
        out.varNameShort = vartxt;
        out.varName = varList{strcmp(vartxt,varShortList)};
        if numel(fileList) > 1
          out.fileName = fileList{strcmp(vartxt,varShortList)};
        end

      else
        out.varName = vartxt;
        varShortList = list_variables(arr(iMatchExact).varNameShort);
        varList      = list_variables(arr(iMatchExact).varName);
        fileList     = list_variables(arr(iMatchExact).fileName);
        out.varNameShort = varShortList{strcmp(vartxt,varList)};
        if numel(fileList) > 1
          out.fileName = fileList{strcmp(vartxt,varList)};
        end
      end
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

%% NESTED FUNCTIONS
  function display_variable_link(iVarList)
    if islogical(iVarList)
      iVarList = find(iVarList);
    end
    for iVar = iVarList
      fprintf('%s\n',['<a href="matlab: psp_var ' arr(iVar).varName '">' arr(iVar).varName '</a>']);
    end
  end

  function display_variable_description(iVarList)
    % if more than 1 variable then print out list
    if islogical(iVarList)
      iVarList = find(iVarList);
    end
    if isempty(iVarList),return; end
    if numel(iVarList) > 1
      display_variable_link(iVarList);
      return;
    end
    for iVar = iVarList
      fprintf('%s\n',C{startIdx(iVar)-1:endIdx(iVar)+1});
    end
  end

  function varList = list_variables(varFilter)
    % LIST_VARIABLE expands regexp variable filter to a list of all possible variable names
    %  Returns cell array with all possible variable names
    % Ex: p(1|2|3) becomes {'p1','p2','p3'}
    tokens  = regexp(varFilter,'\(([\w\|]*)\)','tokens');
    outbase = regexp(varFilter,'\(([\w\|]*)\)','split');
    nTok    = numel(tokens);
    valArr  = cell(size(tokens));
    nVal    = ones(nTok,1);

    % create template
    for iTok = 1 : nTok
      valArr{iTok} = regexp(tokens{iTok}{:},'\|','split');
      nVal(iTok) = numel(valArr{iTok});
    end

    nOutput = prod(nVal);

    varList = cell(nOutput,1);
    varPerm = cell(nOutput,nTok);
    indSort = zeros(nOutput,1);

    for iTok = 1:nTok
      for iVal = 1:nVal(iTok)
        indPerm = iVal:nVal(iTok):nOutput;
        varPerm(indPerm,iTok)=valArr{iTok}(iVal);
        indSort((1:nOutput/nVal(iTok))+nOutput/nVal(iTok)*(iVal-1)) = indPerm;
      end
      varPerm=varPerm(indSort,:);
    end

    for iOut = 1 : nOutput
      for iBase = 1:nTok
        varList{iOut} = [varList{iOut} outbase{iBase} varPerm{iOut,iBase}];
      end
      varList{iOut} = [varList{iOut} outbase{end}];
    end

  end

  function parse_psp_variables_file
    %% split psp_variables.txt into structure array of variable
    % line numnbers of start/end of each structure
    startIdx = find(ismember(C, '%%%%'))+1;
    endIdx   = find(ismember(C, '----'))-1;
    startIdxInfo = find(ismember(C, 'info#start'))+1;
    endIdxInfo   = find(ismember(C, 'info#end'))-1;

    % define array of strucutres
    N = numel(startIdx);
    arr = struct('varName','','fileName','','hourtag',{''},'directory','','varNameShort','','related',{''},'info',"");
    arr = repmat(arr,[N 1]);


    % parse and store each variable in the structure array
    for i=1:N
      % parse key/value of struct
      s = C(startIdx(i):startIdxInfo(i)-2);
      s = regexp(s,'(\w+)\s*[:=]\s*([^%$]*)(?:%[^$]*)?','tokens','once');

      % store: struct.key = value
      for j=1:numel(s)
        keytoken = s{j};
        arr(i).(keytoken{1}) = keytoken{2};
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

  end
end
