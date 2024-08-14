%
% Create UTC datetime object.
%
% Mostly useful for creating automated tests, hence the short name.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function Dt = UTC(varargin)
% PROPOSAL: Automatic test code.
Dt = datetime(varargin{:}, 'TimeZone', 'UTCLeapSeconds');
end
