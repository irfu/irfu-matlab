function y = c_cs_trans(from,to,x,varargin)
%C_CS_TRANS  coordinate transform between GSE/DSC/DSI/ISR2 for Cluster
%
% OUT = C_CS_TRANS(FROM_CS,TO_CS,INP,[ARGS])
%
% Transform INP from FROM_CS to TO_CS.
%
% Input:
%     FROM_CS,TO_CS - one of GSE/DSC/DSI/ISR2
%     INP           - vector data with or without the time column 
% 
% ARGS can be:
%     'SAX'   - plot all x values in separate subplots
%     'CL_ID' - Cluster id 1..4
%     'T'     - Time corresponding to INP (ISDAT epoch or ISO time string)
%
% Examples:
%     B1 = c_cs_trans('DSI','GSE',diB1,'cl_id',1);
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

error(nargchk(4,6,nargin))

if ~ischar(from) || ~ischar(to) || ...
        ~( strcmpi(from,'GSE') || strcmpi(from,'DSC') || ...
        strcmpi(from,'DSI') || strcmpi(from,'ISR2') ) || ...
        ~( strcmpi(to,'GSE') || strcmpi(to,'DSC') || ...
        strcmpi(to,'DSI') || strcmpi(to,'ISR2') )
    error('FROM/TO must be one of GSE,DSC,DSI,ISR2')
end
if strcmpi(from,'ISR2'), from = 'DSI'; end
if strcmpi(to,'ISR2'), to = 'DSI'; end
if strcmpi(from,to)
    error('FROM and TO must be different')
end

lx=size(x,2);
if lx > 3
    inp = x(:,[2 3 4]); % assuming first column is time
    t = x(:,1);
elseif lx == 3
    inp = x;
    t = [];
else
    disp('too few components of vector')
    exit
end

cl_id = [];
sax = [];
have_options = 0;
args = varargin;
if size(args,2) > 1, have_options = 1; end

while have_options
    l = 1;
    if ~ischar(args{1}), error('expecting ''sax'',''cl_id'' or ''t'''),end
    
    switch(lower(args{1}))
        case 'sax'
            if length(args)>1
                if isnumeric(args{2})
                    sax = args{2};
                    if ~all(size(SAX2)==[1 3]), error('size(SAX) ~= [1 3]'), end
                    l = 2;
                else error('SAX value must be numeric')
                end
            else error('SAX value is missing')
            end
            
        case 'cl_id'
            if length(args)>1
                if isnumeric(args{2})
                    cl_id = args{2};
                    if ~any(cl_id == [1 2 3 4]), error('CL_ID value must be 1..4'), end
                    l = 2;
                else error('SCL_IDAX value must be numeric')
                end
            else error('CL_ID value is missing')
            end
            
        case 't'
            if length(args)>1
                if isnumeric(args{2})
                    t = args{2};
                    l = 2;
                elseif ischar(args{2})
                    try
                        t = iso2epoch(args{2});
                        l = 2;
                    catch
                        error('T is not a valid ISO time string')
                    end
                else error('SCL_IDAX value must be an ISDAT epoch or ISO time string')
                end
            else error('CL_ID value is missing')
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
    
    if isempty(sax) % Need to fetch/load spin axis
        caa_load CL_SP_AUX % Load CAA data files
        if exist('CL_SP_AUX','var')
            lat = []; long = [];
            c_eval('lat=getmat(CL_SP_AUX,''sc_at?_lat__CL_SP_AUX'');',cl_id);
            c_eval('long=getmat(CL_SP_AUX,''sc_at?_long__CL_SP_AUX'');',cl_id);
            if (t > lat(1,1)-60) && (t < lat(end,1)+60),
                latlong   = irf_resamp([lat long(:,2)],t);
                [xspin,yspin,zspin] = sph2cart(latlong(3)*pi/180,latlong(2)*pi/180,1);
                sax = [xspin yspin zspin];
                irf_log('dsrc',irf_ssub('Loaded SAX? from CAA file',cl_id))
            end
        end
    end
    
    if isempty(sax)
        [ok,sax] = c_load('SAX?',cl_id); % Load from saved ISDAT files or fetch from ISDAT
        if ok
            irf_log('dsrc',irf_ssub('Loaded SAX? from ISDAT file',cl_id))
            % XXX TODO: check that the SAX is from the right time
            %[iso_t,dt] = caa_read_interval();
        else
            tempv = getData(ClusterDB(c_ctl(0,'isdat_db'),c_ctl(0,'data_path')),...
                t(1),120,cl_id,'sax','nosave');
            if isempty(tempv)
                irf_log('dsrc',irf_ssub('Cannot load SAX?',cl_id))
                sax = [];
            else
                sax = tempv{2};
                irf_log('dsrc',irf_ssub('Loaded SAX? from ISDAT',cl_id))
            end
            clear tempv
        end
    end
    if isempty(sax), error('Cannot get SAX'), end
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
else y=out;
end


