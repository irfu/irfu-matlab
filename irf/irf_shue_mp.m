function r_mp = irf_shue_mp(pos_Re_gsm, bz_nT, swp_nPa)
%IRF_SHUE_MP  estimate distance to model(Shue) magnetopause
%
%  r_mp_Re = irf_shue_mp(pos_Re_gsm, bz_nT, swp_nPa)
%
% Input:
%		pos_Re_gsm - GSM position in Re (3 or 4 components)
%		bz_nT      - IMF Bz in nT
%		swp_nPa    - Solar wind dynamic pressure in nPa
%
% References:
%		Shue et. al., A new functional form to study the solar
% 		wind control of the magnetopause size ans shape,
%		JGR, 102, p.9497, 1997.
%
% See also IRF_GSE2GSM
%

% Copyright 2006 Yuri Khotyaintsev

disp(' ')
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp('irf_shue_mp will be replaced by MODEL.MAGNETOPAUSE_NORMAL')
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp(' ')

if size(pos_Re_gsm,2)>3, pos_Re_gsm = pos_Re_gsm(:,2:4); end

% Shue et. al., Eq. 13
alpha = ( 0.58 -0.01*bz_nT )*( 1.0 +0.01*swp_nPa );

% Shue et. al., Eq. 12
if bz_nT>=0, r0 = ( 11.4 +0.013*bz_nT )*swp_nPa^( -1.0/6.6 );
else,        r0 = ( 11.4 +0.140*bz_nT )*swp_nPa^( -1.0/6.6 );
end

r = irf_abs(pos_Re_gsm,1);
cosTheta = pos_Re_gsm(:,1)./r;
% Shue et. al., Eq. 1
r_mp = r0 *( 2.0./( 1.0 +cosTheta )).^alpha - r;

