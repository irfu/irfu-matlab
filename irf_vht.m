function [vht,eht,dvht,p,cc]=irf_vht(e,b,flag)
% [vht,eht,dvht,p,cc]=irf_vht(e,b,flag)
%  estimate velocity of the deHoffman-Teller frame 
%  from the velocity estimate the electric field eht=-vhtxb
% 
% if flag==2 use version where E=(Ex,Ey,<Not used>) is assumed
% otherwise assumeE=E(Ex,Ey,Ez)
% Assumed units: e [mV/m] b[nT] vht [km/s]
%
% Output:
%   vht - Hofmann Teller frame velocity [km/s]
%   eht - Calculated -vht x b [mV/m]
%  dvht - error of Hofmann Teller frame 
%
% See also IRF_VHT_PLOT
%
% $Id$

if nargin==0, help irf_vht;return;end
if nargin<3, flag=1;end

if size(e,1) ~= size(b,1), b=irf_resamp(b,e);end
%p=av_t_xx(b);
p(1)=sum(b(:,2).*b(:,2))/size(b,1); % Bx*Bx
p(2)=sum(b(:,2).*b(:,3))/size(b,1); % Bx*By
p(3)=sum(b(:,2).*b(:,4))/size(b,1); % Bx*Bz
p(4)=sum(b(:,3).*b(:,3))/size(b,1); % By*By
p(5)=sum(b(:,3).*b(:,4))/size(b,1); % By*Bz
p(6)=sum(b(:,4).*b(:,4))/size(b,1); % Bz*Bz


if (nargin > 2) && (flag == 2), % assume only Ex and Ey
 e(:,4)=0; % put z component to 0 when using only Ex and Ey
 K=[[p(6) 0 -p(3)];[0 p(6) -p(5)];[-p(3) -p(5) p(1)+p(4)]];
 comm= 'Hofmann-Teller frame is calculated using 2 components of E=(Ex,Ey,0)';
else
 K=[[p(4)+p(6) -p(2) -p(3)];[-p(2) p(1)+p(6) -p(5)];[-p(3) -p(5) p(1)+p(4)]];
 comm= 'Hofmann-Teller frame is calculated using all 3 components of E=(Ex,Ey,Ez)';
end

xxExB=irf_cross(e,b);
ind_number=find(~isnan(xxExB(:,2))); % remove NaN from calculation
ExB=sum(xxExB(ind_number,2:4),1)/size(xxExB,1);
VHT=K\ ExB'.*1e3; % 9.12 in ISSI book 
vht=VHT';
strvht=['V_{HT}=' num2str(irf_abs(vht,1),3) ' [ ' num2str(irf_norm(vht),' %5.2f') '] =[' num2str(vht,' %5.2f') '] km/s'];
disp(comm)
disp(strvht);


%
% Calculate the goodness of the Hofmann Teller frame 
% 
eht=irf_e_vxb([0 vht],b); 

if flag == 2,
  ep=[e(ind_number,2);e(ind_number,3)];
  ehtp=[eht(ind_number,2);eht(ind_number,3)];
  delta_e=[e(:,1) e(:,2)-eht(:,2) e(:,3)-eht(:,3) e(:,1)*0];
else
  ep=[e(ind_number,2);e(ind_number,3);e(ind_number,4)];
  ehtp=[eht(ind_number,2);eht(ind_number,3);eht(ind_number,4)];
  delta_e=[e(:,1) e(:,2)-eht(:,2) e(:,3)-eht(:,3) e(:,4)-eht(:,4)];
end
[p,s]=polyfit( ehtp,ep,1);
cc=corrcoef(ep,ehtp);

disp(['slope=' num2str(p(1),3) '  offs=' num2str(p(2),2)]);
disp(['cc=' num2str(cc(1,2),3)]);

%
% Calculate error in velocity estimate
%
delta_xxExB=abs(irf_cross(delta_e,b));
delta_ExB=sum(delta_xxExB(ind_number,2:4),1)/length(delta_xxExB(ind_number,1));
delta_VHT=K\ delta_ExB'.*1e3; % 9.12 in ISSI book 
delta_vht=delta_VHT';dvht=delta_vht;
delta_strvht=['\delta V_{HT}=' num2str(irf_abs(delta_vht,1),3) ' [ ' num2str(irf_norm(delta_vht),' %5.2f') '] =[' num2str(delta_vht,' %5.2f') '] km/s'];
disp(comm)
disp(delta_strvht);

