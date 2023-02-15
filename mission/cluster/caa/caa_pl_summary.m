function caa_pl_summary(sp,op)
%CAA_PL_SUMMARY  make summary plots of the CAA data
%
% caa_pl_summary([sp,options])
%   sp - storage path
%   options :
%     'noexport' do not export plot imaged (PS,PDF,PNG)
%

do_export = 1;

if nargin<1, sp = pwd;
elseif nargin==1
  if strcmp(sp,'noexport'), sp = pwd; do_export = 0; end
elseif nargin==2
  if strcmp(op,'noexport'), do_export = 0;
  else, irf_log('fcal',['unknown option ' op])
  end
end

old_pwd = pwd;
cd(sp);

% Save the screen size
sc_s = get(0,'ScreenSize');
if sc_s(3)==1600 && sc_s(4)==1200, scrn_size = 2;
else, scrn_size = 1;
end

for j=1:4
  figure(j)
  if scrn_size==1 ,set(gcf,'position',[91  40 909 640])
  else, set(gcf,'position',[691   159   909   916])
  end
end

c_pl_summary

if do_export
  [st_s,dt1] = caa_read_interval;
  t1 = iso2epoch(st_s);
  c_export_ps(t1,'pdf')
  c_export_png(t1)
end

cd(old_pwd)
