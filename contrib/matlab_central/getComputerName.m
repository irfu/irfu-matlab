% GETCOMPUTERNAME returns the name of the computer (hostname)
% name = getComputerName()
%
% WARN: output string is converted to lower case
%
% See also SYSTEM, GETENV, ISPC, ISUNIX
%
% m j m a r i n j (AT) y a h o o (DOT) e s
% (c) MJMJ/2007
% MOD: MJMJ/2013

% Author:  Manuel Marin
% Version: 12/12/2013
% Source: Matlab Central File Exchange
% https://www.mathworks.com/matlabcentral/fileexchange/16450-get-computer-name-hostname/

% SPDX-License-Identifier: BSD-2-Clause
% LICENSE
% Copyright (c) 2013, Manuel Marin
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% * Redistributions of source code must retain the above copyright
%   notice, this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright
%   notice, this list of conditions and the following disclaimer in
%   the documentation and/or other materials provided with the distribution
%
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%   POSSIBILITY OF SUCH DAMAGE.
% END OF LICENSE

function name = getComputerName
[err, name] = system('hostname');
if err
  if ispc
    % Windows system
    name = getenv('COMPUTERNAME');
  else
    % Unix system
    name = getenv('HOSTNAME');
  end
end
% Trim and change to lowercase.
name = strtrim(lower(name));
end
