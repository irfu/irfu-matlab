function caa_pl_summary(sp)
%CAA_PL_SUMMARY  make summary plots of the CAA data
%
% caa_pl_summary([sp])
%   sp - storage path
%
% $Id$

if nargin<1, sp=pwd; end

old_pwd = pwd;
cd(sp);

for j=1:4
	figure(j)
	set(gcf,'position',[601   167   999   960])
end

c_pl_summary
c_export_ps
c_export_png

cd(old_pwd)
