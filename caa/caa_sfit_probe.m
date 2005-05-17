function p_res = caa_sfit_probe(cl_id,probe)
%CAA_SFIT_PROBE  set/get probe pair for spin resolution data
%
% caa_sfit_probe(cl_id,probe-p)
%    set probe pair for CL_ID (saves to mInfo.mat)
%    PROBE-P is one of 12, 32, 34
%
% [p_res] = caa_sfit_probe(cl_id)
%    get/display probe pair for CL_ID
%
% $Id$

% Copyright 2005 Yuri Khotyaintsev

error(nargchk(1,2,nargin))
if cl_id<1 | cl_id>4, error('CL_ID must be 1..4'), end

if nargin>1
	if probe~=12 & probe~=32 & probe~=34
		error('PROBE must be 12, 32 or 34')
	end
	c_eval('sfit_probe?=probe;',cl_id)
	if exist('./mInfo.mat','file'), c_eval('save mInfo sfit_probe? -append',cl_id)
	else, c_eval('save mInfo sfit_probe? -append',cl_id)
	end
	irf_log('save',irf_ssub('sfit_probe? -> mInfo.mat',cl_id))
else
	if exist('./mInfo.mat','file')
		warning off
		c_eval('load mInfo sfit_probe?;',cl_id)
		
	end
	if exist(irf_ssub('sfit_probe?',cl_id),'var')
		c_eval('pp=sfit_probe?;',cl_id)
	else
		pp = c_ctl(cl_id,'probe_p');
	end
	
	if nargout>0
		p_res = pp;
	else
		disp(sprintf('C%d sfit probe-p : %d',cl_id,pp))
	end
end
