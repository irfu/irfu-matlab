function phase_out = c_phase(t,phase_2)
%C_PHASE spacecraft phase for give time vector
%
% phase_out = c_phase(t,phase_2)
% phase_out = c_phase(t,ic)
%
% Find spacecraft phase for give time vector t
%
% t - column vector with time in isdat epoch
% phase_2 - column vector [time phase_2]
%   if phase_2 is one number assume that it is sc_number and 
%   download phase from isdat 'disco:10'
% phase - column vector [t phase]
%
% $Id$

if nargin < 2, help c_phase;return;end

t=t(:); % t should be column vector

if prod(size(phase_2))==1, % phase_2 is sc number, download phase_2 from isdat 
  ic=phase_2;tint=[min(t) max(t)];
  [t_isdat,ph2]=caa_is_get('disco:10',fromepoch(tint(1)),tint(2)-tint(1),ic,'ephemeris', 'phase_2');
  phase_2=[t_isdat ph2];
end  % 

phase_out=[t interp1q(phase_2(:,1),phase_2(:,2),t)];

