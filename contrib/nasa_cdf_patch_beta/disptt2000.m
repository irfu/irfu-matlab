function disptt2000(obj)
%DISPTT2000   DISP for CDFTT2000 object.

% If obj is not scalar, then just display the size
s = size(obj);
if ~isequal(s,[1 1])
    disp(sprintf(['     [%dx%d cdftt2000]'], s(1), s(2)));
else
  if (isa(obj, 'double'))
    disp(sprintf('     %s',[datestr(obj,'dd-mmm-yyyy HH:MM:SS.FFF')]));
  else
    tt2000 = spdfencodett2000(obj);
    disp(sprintf('     %s',tt2000{:}));
  end
end
