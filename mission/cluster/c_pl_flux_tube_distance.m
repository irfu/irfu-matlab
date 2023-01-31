function [d1,d2,d3,d4]=c_pl_flux_tube_distance(r1,r2,r3,r4,b1,b2,b3,b4,n_skip) %#ok<INUSL>
%C_PL_FLUX_TUBE_DISTANCE  estimate the distance between flux tubes
%
%   [d1,d2,d3,d4]= C_PL_FLUX_TUBE_DISTANCE(r1,r2,r3,r4,b1,b2,b3,b4,n_skip);
%     Estimate the distance between flux tubes on which Cluster spacecraft
%     are located
%     d1..d4,r1..r4.b1..b4 - column vectors with first column time
%     The time resolution and time interval is given by b1 sampling rate
%     b1..b4 - resampled at b1 sampling
%     n_skip - plot only every n_skip point
%   [d1,d2,d3,d4]= C_PL_FLUX_TUBE_DISTANCE('R?','B?',n_skip);
%     load R1 .. R4 from mR and B1 .. B4 from the environment
%   [d1,d2,d3,d4]= C_PL_FLUX_TUBE_DISTANCE(n_skip);
%     if only n_skip is given load R1 .. R4 from mR and B1 .. B4 from mB
%

flag_distance_along_B=1; % to plot also distance along B

if nargin == 1
  n_skip = r1;
  load mR R1 R2 R3 R4; r1=R1;r2=R2;r3=R3;r4=R4; %#ok<NASGU>
  load mB B1 B2 B3 B4; b1=B1;b2=B2;b3=B3;b4=B4; %#ok<NASGU>
elseif nargin == 3
  r_s = r1;
  b_s = r2;
  n_skip = r3;
  for cl_id=1:4
    ttt = evalin('caller',irf_ssub(r_s,cl_id),'[]'); %#ok<NASGU>
    eval(irf_ssub('r? =ttt;',cl_id)); clear ttt
  end
  for cl_id=1:4
    ttt = evalin('caller',irf_ssub(b_s,cl_id),'[]');  %#ok<NASGU>
    eval(irf_ssub('b? =ttt;',cl_id)); clear ttt
  end
  clear r_s b_s
end

if nargin<1, disp('See usage:');help c_pl_flux_tube_distance; return; end
if nargin==9 || nargin == 1
  ind=1:n_skip:size(b1,1);
  b1=b1(ind,:);
end


for ic=1:4,eval(irf_ssub('bn?=irf_norm(irf_resamp(b?,b1,''linear''));',ic)),end
for ic1=1:4
  for ic2=1:4
    eval(irf_ssub('r!r?=irf_add(1,irf_resamp(r!,b1,''linear''),-1,irf_resamp(r?,b1,''linear''));',ic1,ic2));
    eval(irf_ssub('d!d?=irf_abs(irf_cross(r!r?,bn?),1);',ic1,ic2));
    if flag_distance_along_B==1
      eval(irf_ssub('zd!d?=abs(irf_dot(r!r?,bn?,1));',ic1,ic2));
    end
  end
end

figure;clf;
if flag_distance_along_B==1
  subplot(2,1,1);
end

ccol=['k','r','g','b'];
h=plot(b1(:,1),d2d1,'-ko','MarkerFaceColor','r','MarkerSize',2);
title('Distance between flux tubes on which sc are located');
ylabel('R_a(dots)-R_B(line) [km]');
irf_pl_info(['c\_pl\_flux\_tube\_distance() ' char(datetime("now","Format","dd-MMM-uuuu HH:mm:ss"))]); % add information to the plot

grid on
hold on
for ic1=1:4
  for ic2=1:4
    eval(irf_ssub('plot(b1(:,1),d!d?,[''-'',ccol(?),''o''],''MarkerFaceColor'',ccol(!),''MarkerSize'',2);',ic1,ic2));
  end
end
irf_timeaxis(gca,'date');

if flag_distance_along_B==1
  subplot(2,1,2);
  h=plot(b1(:,1),zd2d1,'-ko','MarkerFaceColor','r','MarkerSize',2);
  title('Distance between spacecraft along B');
  ylabel('R_a(dots)-R_B(line) [km]');
  
  grid on
  hold on
  for ic1=1:4
    for ic2=1:4
      eval(irf_ssub('plot(b1(:,1),zd!d?,[''-'',ccol(?),''o''],''MarkerFaceColor'',ccol(!),''MarkerSize'',2);',ic1,ic2));
    end
  end
  irf_timeaxis(gca,'date');
  
end


