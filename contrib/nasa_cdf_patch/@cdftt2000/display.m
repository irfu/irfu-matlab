function display(obj)
%DISPLAY   DISPLAY for CDFTT2000 object.

if strcmp(get(0, 'FormatSpacing'), 'loose')
    loose = 1;
else
    loose = 0;
end;

%
% Name or ans
%
if loose ~= 0
    disp(' ');
end;

if isempty(inputname(1))
    disp('ans =');
else
    disp([inputname(1) ' =']);
end;

if loose ~= 0
    disp(' ');
end;

tt2000 = todatenum(obj);
if isempty(tt2000)
    disp('     Empty cdftt2000 object');
    return;
elseif isequal(size(obj), [1 1])
    disp('     cdftt2000 object:');
end
disptt2000(tt2000);
