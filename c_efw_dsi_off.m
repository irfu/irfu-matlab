function [Ddsi,Damp] = c_efw_dsi_off(t,cl_id)
%C_EFW_DSI_OFF  get EFW offsets 
%
% [Ddsi,Damp] = c_efw_dsi_off(t,[cl_id])
%
% Ddsi is complex: Dx = real(Ddsi), Dy = imag(Ddsi)
%
% See also CAA_COROF_DSI
%
% $Id$

% Copyright 2007 Yuri Khotyaintsev (yuri@irfu.se)

error(nargchk(1,2,nargin))

Damp = 1.1;

if     t>=toepoch([2005 07 01 00 0 0]), Ddsi = [ .31      .60 .52  .64 ];
elseif t>=toepoch([2005 03 01 00 0 0]), Ddsi = [ .35+0.2i .78 .51  .62 ];
elseif t>=toepoch([2004 11 01 00 0 0]), Ddsi = [ .35      .78 .51  .62 ];
elseif t>=toepoch([2004 05 01 00 0 0]), Ddsi = [ .35      .78 .41  .72 ];
elseif t>=toepoch([2004 01 01 00 0 0]), Ddsi = [ .23      .72 .50  .92 ];
elseif t>=toepoch([2003 01 01 00 0 0]), Ddsi = [ .15      .53 .47  .71 ];
elseif t>=toepoch([2002 01 01 00 0 0]), Ddsi = [ .33      .69 .73  .92 ];
elseif t>=toepoch([2001 07 01 00 0 0]), Ddsi = [ .47      .82 .89 1.04 ];
elseif t>=toepoch([2001 06 01 00 0 0]), Ddsi = [ .31      .60 .52  .64 ];
elseif t>=toepoch([2001 05 25 00 0 0]), Ddsi = [ .31     1.35 .52 1.55 ];
elseif t>=toepoch([2001 04 25 00 0 0]), Ddsi = [ .31      .60 .52  .64 ];
elseif t>=toepoch([2001 02 01 00 0 0]), Ddsi = [ .48      .77 .44 1.11 ];
else
	Ddsi = [ 0 0 0 0];
end

if nargin > 1
	Ddsi = Ddsi(cl_id);
else
	Damp = Damp*ones(1,4);
end
	
	
				
				
				
