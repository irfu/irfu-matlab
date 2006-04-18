function res = caa_str2errid(id)
%CAA_STR2ERRID convert error ID string to number
%
% res = caa_str2errid(id_s)
%
% $Id$

switch id
case 'no_tm', res = 0;
case 'bad_tm', res = 1;
case 'bad_data', res = 2;
case 'no_p1', res = 11;
case 'no_p2', res = 12;
case 'no_p3', res = 13;
case 'no_p4', res = 14;
case 'no_p12', res = 112;
case 'no_p23', res = 123;
case 'no_p34', res = 134;
case 'no_10Hz_filt', res = 255;
case 'no_spin_fits', res = 260;
case 'spec_bias', res = 261;
case 'TDB', res = 999;
otherwise
	error('Unknown problem ID')
end
