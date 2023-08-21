function [vht,eht,dvht,p,cc]=irf_vht(e,b,flag)
% [vht,eht,dvht,p,cc]=irf_vht(e,b,flag)
%  estimate velocity of the De Hoffmann-Teller frame
%  from the velocity estimate the electric field eht=-vhtxb
%
% if flag==2 use version where E=(Ex,Ey,<Not used>) is assumed
% otherwise assumeE=E(Ex,Ey,Ez)
% Assumed units: e [mV/m] b[nT] vht [km/s]
%
% if flag==1, use version where E=(<not used>,Ey,<not used>)
% is assumed. Otherwise assume E=E(Ex,Ey,Ez)
% Assumed units: e [mV/m] b[nT] vht [km/s]
%
% Output:
%   vht - De Hoffmann Teller frame velocity [km/s]
%   eht = - vht x b [mV/m]
%  dvht - error of De Hoffmann Teller frame
%
% See also IRF_VHT_PLOT

if nargin==0, help irf_vht;return;end
if nargin<3, flag=1;end

%% Version with TSeries input
if isa(e,'TSeries') && isa(b,'TSeries')
  nSamples = e.length;
  if nSamples ~= b.length, b=b.resample(e.time);end
  bx = b.x.data;
  by = b.y.data;
  bz = b.z.data;
  ex = e.x.data;
  ey = e.y.data;
  ez = e.z.data;

  p(1)=sum(bx .* bx)/nSamples; % Bx*Bx
  p(2)=sum(bx .* by)/nSamples; % Bx*By
  p(3)=sum(bx .* bz)/nSamples; % Bx*Bz
  p(4)=sum(by .* by)/nSamples; % By*By
  p(5)=sum(by .* bz)/nSamples; % By*Bz
  p(6)=sum(bz .* bz)/nSamples; % Bz*Bz


  if (nargin > 2) && (flag == 2) % assume only Ex and Ey
    z=0; % put z component to 0 when using only Ex and Ey
    K=[[p(6) 0 -p(3)];[0 p(6) -p(5)];[-p(3) -p(5) p(1)+p(4)]];
    comm= 'De Hoffmann-Teller frame is calculated using 2 components of E=(Ex,Ey,0)';
  elseif (nargin > 2) && (flag == 1)
    x=0; z=0;

    K=[p(6),-p(3);-p(3),p(1)];
    comm= {'De Hoffmann-Teller frame is calculated using 1 component of E=(0,Ey,0)';'Output velocities [vx, vz], as vy cannot be calculated (assumed = 0)'};
  else
    K=[[p(4)+p(6) -p(2) -p(3)];[-p(2) p(1)+p(6) -p(5)];[-p(3) -p(5) p(1)+p(4)]];
    comm= 'De Hoffmann-Teller frame is calculated using all 3 components of E=(Ex,Ey,Ez)';
  end
  ExB=cross(e,b);
  indData=find(~isnan(ExB.x.data)); % exclude NaN from calculation
  %   revised by Wenya LI; 2015-11-21, wyli @ irfu
  tmp1 = ExB(indData);
  averExB=sum(tmp1.data,1)/nSamples;
  if (nargin > 2) && (flag == 1)
    averExB(:,2)=[];
  end

  % averExB=sum(ExB(indData).data,1)/nSamples;
  %   end revise.
  VHT=K\ averExB'.*1e3; % 9.12 in ISSI book
  vht=VHT';
  strvht=['V_{HT}=' num2str(irf_abs(vht,1),3) ' [ ' num2str(irf_norm(vht),' %5.2f') '] =[' num2str(vht,' %5.2f') '] km/s'];
  disp(comm)
  disp(strvht);


  %
  % Calculate the goodness of the Hofmann Teller frame
  %
  if flag==1
    vht=[vht(1),0,vht(2)];
  end
  eht=irf_e_vxb(vht,b);

  if flag == 2
    ep=e(indData);
    ehtp=eht(indData);
    ep.data(:,3) = 0;
    ehtp.data(:,3) = 0;
    deltaE=ep.data -ehtp.data;
    [p,s]=polyfit( ehtp.data,ep.data,1);
    cc=corrcoef(ep.data,ehtp.data);
  elseif flag == 1
    ep=e(indData);
    ehtp=eht(indData);
    ep.data(:,1)=0;
    ep.data(:,3)=0;
    ehtp.data(:,1)=0;
    ehtp.data(:,3)=0;

    deltaE=ep.data -ehtp.data;
    [p,s]=polyfit( ehtp.y.data,ep.y.data,1);
    cc=corrcoef(ep.y.data,ehtp.y.data);

  else
    ep=e(indData);
    ehtp=eht(indData);

    deltaE=ep.data -ehtp.data;
    [p,s]=polyfit( ehtp.data,ep.data,1);
    cc=corrcoef(ep.data,ehtp.data);
  end

  disp(['slope=' num2str(p(1),3) '  offs=' num2str(p(2),2)]);
  disp(['cc=' num2str(cc(1,2),3)]);

  %
  % Calculate error in velocity estimate
  %
  % 9.16 in ISSI book
  DVHT=sum(irf_abs(deltaE,1).^2)/length(indData);
  lambda=eig(K);
  S=(DVHT/(2*length(indData)-3))*inv(K);
  if flag==1
    dvht(1)=sqrt([1 0]*S*[1;0])*1e3;
    dvht(2)=sqrt([0 1]*S*[0;1])*1e3;
    strdvht=['\delta V_{HT}=' num2str(irf_abs(dvht,1),3) ' [ ' num2str(irf_norm(dvht),' %5.2f') '] =[' num2str(dvht,' %5.2f') '] km/s'];

  else
    dvht(1)=sqrt([1 0 0]*S*[1;0;0])*1e3;
    dvht(2)=sqrt([0 1 0]*S*[0;1;0])*1e3;
    dvht(3)=sqrt([0 0 1]*S*[0;0;1])*1e3;
    % delta_xxExB=abs(irf_cross(delta_e,b));
    % delta_ExB=sum(delta_xxExB(ind_number,2:4),1)/length(delta_xxExB(ind_number,1));
    % delta_VHT=K\ delta_ExB'.*1e3; % 9.12 in ISSI book
    % delta_vht=delta_VHT';dvht=delta_vht;
    strdvht=['\delta V_{HT}=' num2str(irf_abs(dvht,1),3) ' [ ' num2str(irf_norm(dvht),' %5.2f') '] =[' num2str(dvht,' %5.2f') '] km/s'];
  end
  disp(comm)
  disp(strdvht);
  return
end

%% Old version with numerical input
if size(e,1) ~= size(b,1), b=irf_resamp(b,e);end
p(1)=sum(b(:,2).*b(:,2))/size(b,1); % Bx*Bx
p(2)=sum(b(:,2).*b(:,3))/size(b,1); % Bx*By
p(3)=sum(b(:,2).*b(:,4))/size(b,1); % Bx*Bz
p(4)=sum(b(:,3).*b(:,3))/size(b,1); % By*By
p(5)=sum(b(:,3).*b(:,4))/size(b,1); % By*Bz
p(6)=sum(b(:,4).*b(:,4))/size(b,1); % Bz*Bz


if (nargin > 2) && (flag == 2) % assume only Ex and Ey
  e(:,4)=0; % put z component to 0 when using only Ex and Ey
  K=[[p(6) 0 -p(3)];[0 p(6) -p(5)];[-p(3) -p(5) p(1)+p(4)]];
  comm= 'Hofmann-Teller frame is calculated using 2 components of E=(Ex,Ey,0)';
else
  K=[[p(4)+p(6) -p(2) -p(3)];[-p(2) p(1)+p(6) -p(5)];[-p(3) -p(5) p(1)+p(4)]];
  comm= 'Hofmann-Teller frame is calculated using all 3 components of E=(Ex,Ey,Ez)';
end

xxExB=irf_cross(e,b);
indData=find(~isnan(xxExB(:,2))); % exclude NaN from calculation
xxExB=sum(xxExB(indData,2:4),1)/size(xxExB,1);
VHT=K\ xxExB'.*1e3; % 9.12 in ISSI book
vht=VHT';
strvht=['V_{HT}=' num2str(irf_abs(vht,1),3) ' [ ' num2str(irf_norm(vht),' %5.2f') '] =[' num2str(vht,' %5.2f') '] km/s'];
disp(comm)
disp(strvht);


%
% Calculate the goodness of the Hofmann Teller frame
%
eht=irf_e_vxb([0 vht],b);

if flag == 2
  ep=[e(indData,2);e(indData,3)];
  ehtp=[eht(indData,2);eht(indData,3)];
  deltaE=[e(:,1) e(:,2)-eht(:,2) e(:,3)-eht(:,3) e(:,1)*0];
else
  ep=[e(indData,2);e(indData,3);e(indData,4)];
  ehtp=[eht(indData,2);eht(indData,3);eht(indData,4)];
  deltaE=[e(:,1) e(:,2)-eht(:,2) e(:,3)-eht(:,3) e(:,4)-eht(:,4)];
end
[p,s]=polyfit( ehtp,ep,1);
cc=corrcoef(ep,ehtp);

disp(['slope=' num2str(p(1),3) '  offs=' num2str(p(2),2)]);
disp(['cc=' num2str(cc(1,2),3)]);

%
% Calculate error in velocity estimate
%
% 9.16 in ISSI book
DVHT=sum(irf_abs(deltaE(indData,:),1).^2)/length(indData);
lambda=eig(K);
S=(DVHT/(2*length(indData)-3))/K;
dvht(1)=sqrt([1 0 0]*S*[1;0;0])*1e3;
dvht(2)=sqrt([0 1 0]*S*[0;1;0])*1e3;
dvht(3)=sqrt([0 0 1]*S*[0;0;1])*1e3;
% delta_xxExB=abs(irf_cross(delta_e,b));
% delta_ExB=sum(delta_xxExB(ind_number,2:4),1)/length(delta_xxExB(ind_number,1));
% delta_VHT=K\ delta_ExB'.*1e3; % 9.12 in ISSI book
% delta_vht=delta_VHT';dvht=delta_vht;
strdvht=['\delta V_{HT}=' num2str(irf_abs(dvht,1),3) ' [ ' num2str(irf_norm(dvht),' %5.2f') '] =[' num2str(dvht,' %5.2f') '] km/s'];
disp(comm)
disp(strdvht);

