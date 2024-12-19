function e = R_c_despin(es,phase,coef,db)
%
% the original function is c_despin
% modified by Robert in the summer of -03
% the change is made in the choice of databas: from disco:10 to unix:99
% disp on 57,96 and 98 has been supprest
% open and close db has been replaced by an open db in function call
% ceil has been added to row 58 and floor has been added the row 59
%
% function e = c_despin(es,phase,coef)
% function e = c_despin(es,phase,spacecraft_number)
% function e = c_despin(es,phase,flag)
% function e = c_despin(es,spacecraft_number,flag)
%                   get time from es and use corresponding calibration coef
% e =[t Ex_DSC Ey_DSC Ez_DSC] - electric field in despinned satellite reference frame;
% es=[t p12 p34] - electric field in satellite reference frame;
% es=[t WEC_X WEC_Y WEC_Z] - vector in WEC coordinates
% es=[t SAT_X SAT_Y SAT_Z] - vector in SAT coordinates if flag=='sat'
% phase = [t phase]
% t - time
% coef - calibration coefficients [[A_12 E_offs_12_s E_offs_12_xy];[A_34 E_offs_34_s E_offs_34_xy]]
%        A_12 = Real (relative amplitude error)
%        E_offs_12_s = Real (p12 boom offset)
%        E_offs_12_xy = Complex (real part tells E offset in DSC_X and imaginary in DSC_Y)
%
% !Phase is calculated as a linear fit to the data, to make it fast and simple
% In some case with many and large data gaps this can fail.
%
% despinning requires callibration coeffcients
% the despin algorithm is
%  A_12*(((p12+E_offs_12_s)-> rotate into DSC )+E_offs_12_xy)
%  A_34*(((p34+E_offs_34_s)-> rotate into DSC )+E_offs_34_xy)
%  ---------------------- add both
%  = total field (complex, real along DSC_X and imaginary along DSC_Y)
%
t=es(:,1);
if nargin == 2
  coef=[[1 0 0];[1 0 0]];
end

if nargin == 4
  if isnumeric(coef)
    ref_frame='wec';
    if size(coef,1) == 1
      ic=coef;
      [c1,c2,c3,c4]=c_efw_calib(es(1,1)); %#ok<ASGLU>
      clear coef;
      eval(av_ssub('coef=c?;',ic));
    end
  elseif strcmp(coef,'sat')
    ref_frame='sat';
    coef=[[1 0 0];[1 0 0]];
  end
end

if numel(phase)==1 % load phase from isdat database
  ic=phase;phase=[]; %disp(['load phase for sc' num2str(ic)]);
  start_time=ceil(fromepoch(es(1,1))); % time of the first point
  Dt=floor(es(end,1)-es(1,1)+1);
  %  db = Mat_DbOpen('unix:99');
  [phase_t,phase_data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'phase', ' ', ' ', ' ');
  phase=[double(phase_t) double(phase_data)]; clear phase_t phase_data;
  %  Mat_DbClose(db);
end

switch ref_frame
  case 'wec'
    phi_12=3*pi/4;phi_34=pi/4; % angles when phase =0
    if size(es,2)==3
      p12=es(:,2);p34=es(:,3);
    elseif size(es,2)==4
      p12=es(:,4);p34=es(:,3);
    end
  case 'sat'
    phi_12=pi/2;phi_34=0; % angles when phase =0
    if size(es,2)==3
      p12=es(:,2);p34=es(:,3);
    elseif size(es,2)==4
      p12=es(:,3);p34=es(:,2);
    end
end

%contPhase=unwrap(double(real(phaseVal))/180*pi);
ph=phase;ph(:,1)=ph(:,1)-phase(1,1);
phc=unwrap(ph(1:2,2)/180*pi);
phc_coef=polyfit(ph(1:2,1),phc,1);
for j=1:floor(log10(length(ph(:,1))))
  ii=10^j;
  dp=ph(ii,2)-mod(polyval(phc_coef,ph(ii,1))*180/pi,360);
  dpm=[dp dp-360 dp+360];in=find(abs(dpm)<180);
  dph=dp(in);
  phc_coef(1)=phc_coef(1)+dph*pi/180/ph(ii,1);
end
dphc=exp(1i*ph(:,2)/180*pi)-exp(1i*polyval(phc_coef,ph(:,1)));
dd=dphc.*conj(dphc);
err=sum(dd)/length(dd);
if err>.001
  % disp(strcat('Be careful, can be problems with despin, err=',num2str(err)))
end
%disp(strcat('rotation period=',num2str(2*pi/phc_coef(1)),' s'));

phase=polyval(phc_coef,t-phase(1,1));
% take away offsets in the satellite ref frame
p12=p12-coef(1,2);
p34=p34-coef(2,2);

% take away dc offsets of the despinned ref frame
p12=p12-abs(coef(1,3))*cos(angle(coef(1,3))-phase-phi_12);
p34=p34-abs(coef(2,3))*cos(angle(coef(2,3))-phase-phi_34);

% despin each component separately
dp12=p12.*exp(1i*phase)*exp(1i*phi_12);
dp34=p34.*exp(1i*phase)*exp(1i*phi_34);

% add necessary scaling
dp12=coef(1,1)*dp12;
dp34=coef(2,1)*dp34;

% create the final field
e=[t real(dp12+dp34) imag(dp12+dp34)];

if size(es,2)==4
  switch ref_frame
    case 'wec'
      e(:,4)=es(:,2);
    case 'sat'
      e(:,4)=es(:,4);
  end
end
