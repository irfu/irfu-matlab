%script remove_problems
%
% Input: signal,probe,cl_id,problems
% Output: res
%
% $Id$

% Copyright 2005-2007 Yuri Khotyaintsev

if ~exist('CAA_MODE','var'), CAA_MODE = c_ctl(0,'caa_mode'); end

res = signal;
if probe>10
	switch probe
	case 12
		p_list = [1,2];
	case 32
		p_list = [3,2];
	case 34
		p_list = [3,4];
	otherwise
		error('Unknown probe')
    end
elseif probe>0 && probe <=4, p_list = probe;
else
	error('Unknown probe')
end

param = tokenize(problems,'|');
for i=1:length(param)

	switch lower(param{i})
		case 'reset' 
			% Remove bad bias around EFW reset
			[ok,bbias] = c_load('BADBIASRESET?',cl_id);
			if ok
				if ~isempty(bbias)
					irf_log('proc','blanking bad bias due to EFW reset')
					res = caa_rm_blankt(res,bbias);
				end
			else
				if CAA_MODE, error(irf_ssub('Cannot load BADBIASRESET?',cl_id)), end
				irf_log('load',irf_ssub('Cannot load BADBIASRESET?',cl_id))
			end
			clear ok bbias
			
		case 'bbias' 
			% Remove bad bias from bias current indication
			for kk = p_list
				[ok,bbias] = c_load(irf_ssub('BADBIAS?p!',cl_id,kk));
				if ok
					if ~isempty(bbias)
						irf_log('proc',['blanking bad bias on P' num2str(kk)])
						res = caa_rm_blankt(res,bbias);
					end
				else
					if CAA_MODE, error(irf_ssub('Cannot load BADBIAS?p!',cl_id,kk)), end
					irf_log('load',irf_ssub('Cannot load BADBIAS?p!',cl_id,kk))
				end
				clear ok bbias
			end
			
		case 'probesa' 
			% Remove probe saturation
			for kk = p_list
				[ok, sa] = c_load(irf_ssub('PROBESA?p!',cl_id,kk));
				if ok
					if ~isempty(sa)
						irf_log('proc',['blanking saturated P' num2str(kk)])
						res = caa_rm_blankt(res,sa);
					end
				else
					if CAA_MODE, error(irf_ssub('Cannot load PROBESA?p!',cl_id,kk)), end
					irf_log('load',irf_ssub('Cannot load PROBESA?p!',cl_id,kk))
				end
				clear ok sa
			end
			
		case 'probeld' 
			% Remove probe saturation due to low density
			for kk = p_list
				[ok, sa] = c_load(irf_ssub('PROBELD?p!',cl_id,kk));
				if ok
					if ~isempty(sa)
						irf_log('proc',...
							['blanking low density saturation on P' num2str(kk)])
						res = caa_rm_blankt(res,sa);
					end
				else
					if CAA_MODE, error(irf_ssub('Cannot load PROBELD?p!',cl_id,kk)), end
					irf_log('load',irf_ssub('Cannot load PROBELD?p!',cl_id,kk))
				end
				clear ok sa
			end
			
		case 'whip' 
			% Remove whisper pulses
			[ok,whip] = c_load('WHIP?',cl_id);
			if ok
				if ~isempty(whip)
					irf_log('proc','blanking Whisper pulses')
					res = caa_rm_blankt(res,whip);
				end
			else
				if CAA_MODE, error(irf_ssub('Cannot load WHIP?',cl_id)), end
				irf_log('load',	irf_ssub('Cannot load WHIP?',cl_id))
			end
			clear ok whip
			
		case 'sweep' 
			% Remove sweeps
			[ok,sweep] = c_load('SWEEP?',cl_id);
			if ok
				if ~isempty(sweep)
					irf_log('proc','blanking sweeps')
					res = caa_rm_blankt(res,sweep);
					clear sweep
				end
			else
				if CAA_MODE, error(irf_ssub('Cannot load SWEEP?',cl_id)), end
				irf_log('load', irf_ssub(['Cannot load SWEEP?'],cl_id))
			end
			clear ok sweep
			
		case 'bdump' 
			% Remove burst dumps
			[ok,bdump] = c_load('BDUMP?',cl_id);
			if ok
				if ~isempty(bdump)
					irf_log('proc','blanking burst dumps')
					res = caa_rm_blankt(res,bdump);
					clear bdump
				end
			else
				if CAA_MODE, error(irf_ssub('Cannot load BDUMP?',cl_id)), end
				irf_log('load', irf_ssub(['Cannot load BDUMP?'],cl_id))
			end
			clear ok bdump
			
		otherwise
			error('Unknown parameter')
    end
end
