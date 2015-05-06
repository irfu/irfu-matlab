function disptt2000(obj)
%DISPTT2000   DISP for CDFTT2000 object.

% If obj is not scalar, then just display the size
s = size(obj);
if ~isequal(s,[1 1])
    disp(sprintf(['     [%dx%d cdftt2000]'], s(1), s(2)));
else
    disp( [spdfencodett2000(obj)]);
end
