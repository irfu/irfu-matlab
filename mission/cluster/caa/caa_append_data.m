function data = caa_append_data(data,app)
%CAA_APPEND_DATA  concatenate two datasets
%
% res = caa_append_data(data,app)
%       append APP to the end of DATA and fill the gap in between
%
% See also CAA_FILL_GAPS
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(2,2)

if size(data,2) ~= size(app,2), error('data has a different dimension'), end

if isstruct(data)
  if ~isstruct(app), error('APP must be also a struct, as DATA'), end
  if ~isfield(data,'t'), error('struct DATA must have field ''t'''), end
  t = data.t;
  tapp = app.t;
else
  t = data(:,1);
  tapp = app(:,1);
end

if tapp(1) <= t(end)
  if tapp(1) == t(end)
    %error('Last point in DATA has same timestamp as first point in APP')
    if length(tapp) == 1, return, end
    app=app(2:end,:);
  else
    app = irf_tlim(app,tapp(1),t(end),1);
    if isempty(app), return, end
  end
end

if isstruct(data)
  fn = fieldnames(data);
  ndata = length(t);
  for fi=1:length(fn)
    if iscell(data.(fn{fi}))
      for i = 1:length(data.(fn{fi}))
        if size(data.(fn{fi}){i},1) == ndata
          data.(fn{fi}){i} = [data.(fn{fi}){i}; app.(fn{fi}){i}];
        end
      end
    else
      if size(data.(fn{fi}),1) == ndata
        data.(fn{fi}) = [data.(fn{fi}); app.(fn{fi})];
      end
    end
  end
else
  % Fill the gap between the datasets
  data = caa_fill_gaps(data,app(1,1));
  
  data = [data; app];
end
