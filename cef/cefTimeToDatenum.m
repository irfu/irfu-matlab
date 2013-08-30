% 
% Conversion from CEF time format to standard matlab format.
% Ex:
%     cefTimeToDatenum('2001-07-09T06:19:28.521Z')
% 
%
%
function ret=cefTimeToDatenum(d)

ret=cefMjsToDatenum(cefTimeToMjs(d));
 