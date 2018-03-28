function dirsave(file, varargin)
%DIRSAVE Like save, but each var in its own file
%
% dirsave filename var1 var2 var3...
if nargin < 1 || isempty(file); file = 'matlab';  end
[tfStruct,loc] = ismember({'-struct'}, varargin);
args = varargin;
args(loc(tfStruct)) = [];
if ~all(cellfun(@isvarname, args))
    error('Invalid arguments. Usage: dirsave filename <-struct> var1 var2 var3 ...');
end
if tfStruct
    structVarName = args{1};
    s = evalin('caller', structVarName);
else
    varNames = args;
    if isempty(args)
        w = evalin('caller','whos');
        varNames = { w.name };
    end
    captureExpr = ['struct(' ...
        join(',', cellfun(@(x){sprintf('''%s'',{%s}',x,x)}, varNames)) ')'];
    s = evalin('caller', captureExpr);
end

% Use Java checks to avoid partial path ambiguity
jFile = java.io.File(file);
if ~jFile.exists()
    ok = mkdir(file);
    if ~ok 
        error('failed creating dirsave dir %s', file);
    end
elseif ~jFile.isDirectory()
    error('Cannot save: destination exists but is not a dir: %s', file);
end
names = fieldnames(s);
for i = 1:numel(names)
    varFile = fullfile(file, [names{i} '.mat']);
    varStruct = struct(names{i}, {s.(names{i})}); %#ok<NASGU>
    save(varFile, '-struct', 'varStruct');
end

function out = join(Glue, Strings)
Strings = cellstr(Strings);
if isempty( Strings )
    out = '';
elseif length( Strings ) == 1
    out = Strings{1};
else
    Glue = sprintf( Glue ); % Support escape sequences
    out = strcat( Strings(1:end-1), { Glue } );
    out = [ out{:} Strings{end} ];
end
