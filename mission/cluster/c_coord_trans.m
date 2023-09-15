function y = c_coord_trans(from,to,x,varargin)
%C_COORD_TRANS  coordinate transform between GSE/GSM/DSC/DSI/ISR2 for Cluster
%
% OUT = C_COORD_TRANS(FROM_CS,TO_CS,INP,[ARGS])
%
% Transform INP from coordinate system FROM_CS into coordinate system TO_CS.
%
% Input:
%     FROM_CS,TO_CS - one of GSE/GSM/DSC/SR2/DSI/ISR2
%     INP           - vector data with or without the time column
%
% ARGS can be parameter/value pairs with following parameters:
%     'SAX'   - spin axis vector in GSE
%     'CL_ID' - Cluster id 1..4
%     'T'     - Time corresponding to INP (ISDAT epoch or ISO time string)
%
% Note on reference frames:
%       DSI=ISR2 - despinned inversed reference frame
%       DSC=SR2  - despinned reference frame
%
% Examples:
%     B1 = C_COORD_TRANS('DSI','GSE',diB1,'cl_id',1);
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

persistent lat long cl_id_saved

CL_SP_AUX = [];

narginchk(3,8)
if isempty(x) && isa(x,'TSeries'), y = TSeries([]); return;
elseif isempty(x), y = []; return; end % if empty input, empty output

allowed_coord_sys={'GSE','DSC','SR2','DSI','ISR2','GSM'};

if ~any(strcmpi(from,allowed_coord_sys)) || ~any(strcmpi(to,allowed_coord_sys))
  error('FROM/TO must be one of GSE,GSM,DSC,SR2,DSI,ISR2')
end
if strcmpi(from,to) % if both reference frames are the same
  y=x;
  return;
end
if strcmpi(from,'ISR2'), from = 'DSI'; end
if strcmpi(from,'SR2'), from = 'DSC'; end
if strcmpi(to,'ISR2'), to = 'DSI'; end
if strcmpi(to,'SR2'), to = 'DSC'; end
if strcmpi(to,'GSM')
  xgse=c_coord_trans(from,'GSE',x,varargin{:});
  y=irf_gse2gsm(xgse);
  return
end
if strcmpi(from,'GSM')
  xgse=irf_gse2gsm(x,-1);
  y=c_coord_trans('GSE',to,xgse,varargin{:});
  return
end

lx=size(x,2);
if lx > 3
  inp = x(:,[2 3 4]); % assuming first column is time
  t = x(:,1);
elseif lx == 3
  inp = x;
  t = [];
else
  error('c_coord_trans: too few components of vector!')
end

cl_id = [];
sax = [];
have_options = 0;
args = varargin;
if size(args,2) > 1, have_options = 1; end

while have_options
  if ~ischar(args{1}), error('expecting ''sax'',''cl_id'' or ''t'''),end

  switch(lower(args{1}))
    case 'sax'
      if length(args)>1
        if isnumeric(args{2})
          sax = args{2};
          if ~all(size(sax)==[1 3]), error('size(sax) ~= [1 3]'), end
          l = 2;
        else, error('SAX value must be numeric')
        end
      else, error('SAX value is missing')
      end
    case 'cl_id'
      if length(args)>1
        if isnumeric(args{2})
          cl_id = args{2};
          if ~any(cl_id == [1 2 3 4]), error('CL_ID value must be 1..4'), end
          l = 2;
        else, error('CL_ID value must be numeric')
        end
      else, error('CL_ID value is missing')
      end
    case 't'
      if length(args)>1
        if isnumeric(args{2})
          t = args{2};
          l = 2;
        elseif ischar(args{2})
          try
            t = irf_time(args{2},'utc>epoch');
            l = 2;
          catch
            error('T is not a valid ISO time string. ')
          end
        else, error('T value must be an ISDAT epoch or ISO time string')
        end
      else, error('T value is missing')
      end
    otherwise
      error(['unknown parameter : ' args{1}])
  end
  args = args(l+1:end);
  if isempty(args), break, end
end

if strcmpi(from,'GSE') || strcmpi(to,'GSE')

  if isempty(sax) && ( isempty(cl_id) || isempty(t) )
    error('Need both CL_ID and T if SAX is not given and wanting GSE')
  end
  if isempty(sax) % get spin axis from latitude longitude
    [ok,sax] = sax_from_lat_long;
    if ~ok, flagReadLat = 1; else, flagReadLat = 0; end
    if flagReadLat % try to read lat and long from CAA files
      if isempty(CL_SP_AUX)
        caa_load CL_SP_AUX % Load CAA data files
      end
      if ~isempty(CL_SP_AUX)
        lat = []; long = [];
        lat = getmat(CL_SP_AUX,sprintf('sc_at%d_lat__CL_SP_AUX',cl_id));
        long = getmat(CL_SP_AUX,sprintf('sc_at%d_long__CL_SP_AUX',cl_id));
        cl_id_saved=cl_id;
      end
      [ok,sax] = sax_from_lat_long;
      if ~ok, flagReadLat = 1; else, flagReadLat = 0; end
    end
    if flagReadLat % try to stream from CSA
      lat = []; long = [];
      latVarName  = ['sc_at' num2str(cl_id) '_lat__CL_SP_AUX'];
      longVarName = ['sc_at' num2str(cl_id) '_long__CL_SP_AUX'];
      tint = [t(1)-120 t(end)+120];
      out = c_caa_cef_var_get({latVarName,longVarName},'tint',tint,'stream');
      lat = out{1};
      long = out{2};
      cl_id_saved=cl_id;
      [ok,sax] = sax_from_lat_long;
      if ~ok, flagReadLat = 1; else, flagReadLat = 0; end
    end
    if flagReadLat % still no data, put to empty sax
      sax=[];
    end
  end

  if isempty(sax) % could not load anywhere SAX, use default 0 0 1
    disp('!!!!!!!!!!!! ERROR !!!!!!!!!!!!')
    disp('c_coord_trans: could not load SAX variable')
    disp('Using SAX=[0 0 1]; ')
    disp('Coordinate transformations are wrong!')
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    sax=[0 0 1];
  end
  Rx = sax(1);
  Ry = sax(2);
  Rz = sax(3);
  a = 1/sqrt(Ry^2+Rz^2);
  M = [[a*(Ry^2+Rz^2) -a*Rx*Ry -a*Rx*Rz];[0 a*Rz	-a*Ry];[Rx	Ry	Rz]];
end

if strcmpi(from,'GSE') && strcmpi(to,'DSC') % GSE -> DSC
  out = M*inp';
  out = out';
elseif strcmpi(from,'GSE') && strcmpi(to,'DSI') % GSE -> DSI
  out = M*inp';
  out = out';
  out(:,2) = -out(:,2);
  out(:,3) = -out(:,3);
elseif strcmpi(from,'DSC') && strcmpi(to,'GSE')  % DSC -> GSE
  out = M\inp';
  out = out';
elseif strcmpi(from,'DSI') && strcmpi(to,'GSE')  % DSI -> GSE
  inp(:,2) = -inp(:,2);
  inp(:,3) = -inp(:,3);
  out = M\inp';
  out = out';
elseif ( strcmpi(from,'DSI') && strcmpi(to,'DSC') ) ||...  % DSI <-> DSC
    ( strcmpi(from,'DSC') && strcmpi(to,'DSI') )
  inp(:,2) = -inp(:,2);
  inp(:,3) = -inp(:,3);
  out = inp;
else
  error('No coordinate transformation done!')
end

y = x;
if lx > 3, y(:,[2 3 4]) = out; % Assuming first column is time
else, y=out;
end

  function [ok,sax]=sax_from_lat_long
    % check if persistent lat and long variables include required time
    % interval
    ok = false;sax=[];
    if any(lat) % exists saved spin latidude files, check they are ok
      if (cl_id == cl_id_saved) && (t(1) > lat(1,1)-60) && (t(end) < lat(end,1)+60)
        latlong   = irf_resamp([lat long(:,2)],t(1));
        [xspin,yspin,zspin] = sph2cart(latlong(3)*pi/180,latlong(2)*pi/180,1);
        sax = [xspin yspin zspin];
        ok = true;
      end
    end
  end

end
