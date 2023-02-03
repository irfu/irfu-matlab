% DATASTORE is a facility for persisting data across matlab
% sessions. As such, it is useful for matlab applications or
% scripts that need to store permanent settings. The
% following calls demonstrate usual usage.
%
% datastore(app, key, value)
%
%   Write "value" against "key" into the database associated
%   with "app".
%
% (1) out = datastore
%
%   Retrieve all data in the store, in a two-level
%   structure.
%
% (2) out = datastore(app)
%
%   Retrieve all the data associated with "app", in a
%   one-level structure.
%
% (3) out = datastore(app, key)
%
%   Retrieve just one item from the data associated with
%   "app", specified by "key".
%
% The above calls all imply the operation from the argument
% lists. For some other operations, the operation must be
% specified explicitly, as follows.
%
% datatore('@delete', app)
%
%   Delete all data associated with "app".
%
% datatore('@delete', app, key)
%
%   Delete any data stored against "key" in the database
%   associated with "app".
%
% NB: If a value stored against an app/key combination is
% absent, then the read call (3) returns the empty
% matrix, []. Therefore, call (3) returns [] both if the
% value is absent, or if the value is equal to []. Thus, a
% call to @delete only has an impact on what is returned by
% read calls (1) and (2), which only return key/value pairs
% that are actually present in the database.

% Author: Ben Mitch
% Version: 05/12/2011
% Source: Matlab Central File Exchange
% https://www.mathworks.com/matlabcentral/fileexchange/34081-datastore/

% SPDX-License-Identifier: BSD-2-Clause
% LICENSE
% Copyright (c) 2011, Ben Mitch
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

function out = datastore(app, key, value)



%% PREAMBLE

persistent data filename

% get filename
if isempty(filename)
  
  % report
  disp('loading datastore...');
  
  % find path
  if ispc
    path = getenv('APPDATA');
  elseif isunix
    path = '~';
  else
    error('please let me know how to find a suitable path on your system!');
  end
  
  % Get computer name
  hostName = getComputerName;
  % get filename
  filename = [path, '/.matlab_datastore_', hostName];
  
  % load
  if exist(filename, 'file')
    load(filename, '-mat');
  else
    data = struct();
  end
  
end



%% HANDLE OPERATIONS

if nargin >= 1 && ischar(app) && ~isempty(app) && app(1) == '@'
  
  % shift args
  op = app;
  if nargin >= 2
    app = key;
  end
  if nargin >= 3
    key = value;
  end
  
  % switch on operation
  switch op
    
    case '@delete'
      switch nargin
        case 2
          if isfield(data, app)
            data = rmfield(data, app);
            save(filename, 'data');
          end
        case 3
          if isfield(data, app)
            if isfield(data.(app), key)
              data.(app) = rmfield(data.(app), key);
              save(filename, 'data');
            end
          end
        otherwise
          error('bad argument count for @delete');
      end
      
    otherwise
      error(['unrecognised operation "' op '"']);
      
  end
  
  % ok
  return
  
end



%% HANDLE SUMMARY

if nargout == 0 && nargin == 0
  
  % space
  disp(' ')
  
  % for each app
  apps = fieldnames(data);
  for a = 1:length(apps)
    
    % extract
    app = apps{a};
    appdata = data.(app);
    
    % display
    disp(['app "' app '":'])
    disp(appdata)
    
  end
  
  % no apps case
  if isempty(apps)
    disp('No data stored.');
    disp(' ');
  end
  
  % ok
  return
  
end



%% HANDLE READ

if nargout == 1 && nargin <= 2
  
  % get data
  out = data;
  
  % pare by app
  if nargin >= 1
    if isfield(data, app)
      out = data.(app);
    else
      out = struct();
    end
  end
  
  % pare by key
  if nargin >= 2
    if isfield(out, key)
      out = out.(key);
    else
      out = [];
    end
  end
  
  % ok
  return
  
end



%% HANDLE WRITE

if nargout == 0 && nargin == 3
  
  % write data
  data.(app).(key) = value;
  
  % save data
  save(filename, 'data');
  
  % ok
  return
  
end



%% OTHERWISE

error('bad usage of datastore - see "help datastore"');



