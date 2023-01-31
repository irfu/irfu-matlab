function e = c_efw_despin(es,phase,coef,options)
%C_EFW_DESPIN despin EFW data
%
% function e = c_efw_despin(es,phase)
% function e = c_efw_despin(es,phase,coef)
% function e = c_efw_despin(es,phase,spacecraft_number)
% function e = c_efw_despin(es,phase,spacecraft_number,options)
% function e = c_efw_despin(es,phase,options)
% function e = c_efw_despin(es,spacecraft_number)
% function e = c_efw_despin(es,spacecraft_number,options)
%                   get time from es and use corresponding calibration coef
% e =[t Ex_DSC Ey_DSC Ez_DSC] - electric field in despinned satellite reference frame;
% es=[t p12 p34] - electric field in satellite reference frame;
% es=[t WEC_X WEC_Y WEC_Z] - vector in WEC coordinates
% es=[t SAT_X SAT_Y SAT_Z] - vector in SAT coordinates if options=='sat'
% phase = [t phase]
% t - time
% coef - calibration coefficients [[A_12 E_offs_12_s E_offs_12_xy];[A_34 E_offs_34_s E_offs_34_xy]]
%        A_12 = Real (relative amplitude error)
%        E_offs_12_s = Real (p12 boom offset)
%        E_offs_12_xy = Complex (real part tells E offset in DSC_X and imaginary in DSC_Y)
% options - 'efw'      despin from WEC ref frame + use the closest callibration [good for very short or handcalibrated time intervals]
%        'efw_a'    despin from WEC, subtract mean value of probe signals to get rid of offsets
%        'efw_b'    same as 'efw_a' + correct with nearest sunward offset
%        'staff' or 'wec' despin from WEC
%        'sat'      despin from SR
%        'asym'     despin asymmetric probe configuration
%        'interf'   despin p14/p24 data
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
if numel(phase)==1, ic=phase;end  % if only one number then it is sc number
if nargin == 2
  coef=[[1 0 0];[1 0 0]];
  ref_frame='wec';
end

if size(es,2)==3 % if input is [t p12 p34] convert to [t 0 p34 p12]
  es=es(:,[1 3 3 2]);es(:,2)=0;
end
use_asym = 0;
use_interf = 0;
if nargin == 4
  if strcmp(options,'asym') || strcmp(options,'asym32'), use_asym = 32;
  elseif strcmp(options,'asym42'), use_asym = 42;
  elseif strcmp(options,'interf')
    use_interf = 1;
    ref_frame='wec';
  end
end

if nargin >= 3
  if isnumeric(coef)
    ref_frame='wec';
    if size(coef,1) == 1
      ic=coef;
      [c1,c2,c3,c4]=c_efw_calib(es(1,1)); %#ok<NASGU,ASGLU>
      clear coef;
      eval(irf_ssub('coef=c?;',ic));
      if nargin == 4
        if strcmp(options,'efw_b')
          coef(1,2)=mean(es(:,4));
          coef(2,2)=mean(es(:,3));
        elseif strcmp(options,'efw_a')
          coef=[[1 0 0];[1 0 0]];
          coef(1,2)=mean(es(:,4));
          coef(2,2)=mean(es(:,3));
        end
      end
    end
  elseif strcmp(coef,'sat')
    ref_frame='sat';
    coef=[[1 0 0];[1 0 0]];
  elseif strcmp(coef,'wec')
    ref_frame='wec';
    coef=[[1 0 0];[1 0 0]];
  elseif strcmp(coef,'efw')
    ref_frame='wec';
    [c1,c2,c3,c4]=c_efw_calib(es(1,1)); %#ok<NASGU,ASGLU>
    eval(irf_ssub('coef=c?;',ic));
  elseif strcmp(coef,'efw_b')
    ref_frame='wec';
    [c1,c2,c3,c4]=c_efw_calib(es(1,1)); %#ok<NASGU,ASGLU>
    eval(irf_ssub('coef=c?;',ic));
    coef(1,2)=mean(es(:,4));
    coef(2,2)=mean(es(:,3));
  elseif strcmp(coef,'efw_a')
    ref_frame='wec';
    coef=[[1 0 0];[1 0 0]];
    coef(1,2)=mean(es(:,4));
    coef(2,2)=mean(es(:,3));
  elseif strcmp(coef,'staff')
    ref_frame='wec';
    coef=[[1 0 0];[1 0 0]];
  elseif strcmp(coef,'interf')
    use_interf = 1;
    ref_frame='wec';
    coef=[[1 0 0];[1 0 0]];
  end
end

if numel(phase)==1 % load phase from isdat database
  ic=phase;disp(['load phase for sc' num2str(ic)]);
  start_time=fromepoch(es(1,1)); % time of the first point
  Dt=es(end,1)-es(1,1)+1;
  DB_S = c_ctl(0,'isdat_db');
  db = Mat_DbOpen(DB_S);
  [phase_t,phase_data] = isGetDataLite( db, start_time, Dt,'Cluster', num2str(ic), 'ephemeris', 'phase', ' ', ' ', ' ');
  phase=[double(phase_t) double(phase_data)]; clear phase_t phase_data;
  Mat_DbClose(db);
end

if numel(phase) < 3
  tf=t(1,1);
  if isinf(t(1,1)) || isnan(t(1,1))
    tf=0;
  end
  irf_log('proc',['c_efw_despin failed: no phase information for interval starting ' epoch2iso(tf)])
  e=[t zeros(length(t),1)+NaN zeros(length(t),1)+NaN];
  return
end

switch ref_frame
  case 'wec'
    if use_interf
      phi_12 = pi/2;   % p14
      phi_34 = 0; % p24
    else
      phi_12=3*pi/4;
      phi_32=pi/2;
      phi_42=pi;
      phi_34=pi/4; % angles when phase =0
    end
    p12=es(:,4);p34=es(:,3);
  case 'sat'
    phi_12=pi/2;phi_34=0; % angles when phase =0
    p12=es(:,3);p34=es(:,2);
end

%contPhase=unwrap(double(real(phaseVal))/180*pi);
ph=phase;tref=phase(1,1);ph(:,1)=ph(:,1)-tref;
phc=unwrap(ph(1:2,2)/180*pi);
phc_coef=polyfit(ph(1:2,1),phc,1);
for j=1:floor(log10(length(ph(:,1))))
  ii=10^j;
  dp=ph(ii,2)-mod(polyval(phc_coef,ph(ii,1))*180/pi,360);
  dpm=[dp dp-360 dp+360];
  dph=dpm(abs(dpm)<180);
  if isempty(dph)
    continue;
  end
  phc_coef(1)=phc_coef(1)+dph*pi/180/ph(ii,1);
end
%dphc=exp(1i*ph(:,2)/180*pi)-exp(1i*polyval(phc_coef,ph(:,1)));
%dd=dphc.*conj(dphc);
%err=sum(dd)/length(dd);
diffangle=mod(ph(:,2)-polyval(phc_coef,ph(:,1))*180/pi,360);
diffangle=abs(diffangle);
diffangle=min([diffangle';360-diffangle']);
err_angle_mean=mean(diffangle);
err_angle=std(diffangle);
if err_angle>1 || err_angle_mean>1
  irf_log('proc',['Using standard despinning! Polyn. fit despinning errors would be >1deg. err=' num2str(err_angle) 'deg.']);
  unwrap_phase=phase;
  if max(diff(phase(:,1)))>4
    disp('There are data gaps in the phase data.Despinned data can have problems near data gaps.');
    for j=2:length(phase(:,1))
      ddphase=unwrap_phase(j,:)-unwrap_phase(j-1,:);
      if ddphase(2)<0 || ddphase(1)>3.9
        unwrap_phase(j:end,2)=unwrap_phase(j:end,2)+360*round((ddphase(1)*360/4-ddphase(2))/360);
      end
    end
  else % no data gaps in phase
    for j=2:length(phase(:,1))
      ddphase=unwrap_phase(j,:)-unwrap_phase(j-1,:);
      if ddphase(2)<0
        unwrap_phase(j:end,2)=unwrap_phase(j:end,2)+360;
      end
    end
  end
  phase=interp1q(unwrap_phase(:,1),unwrap_phase(:,2),t)*pi/180;
else % use polynomial fit
  % disp(['Despinning using polynomial fit. Average deviation ' num2str(err_angle_mean) 'deg. Standard deviation ' num2str(err_angle) ' deg.']);
  phase=polyval(phc_coef,t-tref);
  % if err>.001
  %  disp(strcat('Be careful, can be problems with despin, err=',num2str(err)))
  % end
end

irf_log('proc',strcat('rotation period=',num2str(2*pi/phc_coef(1)),' s'));

% take away offsets in the satellite ref frame
p12=p12-coef(1,2);
p34=p34-coef(2,2);

% take away dc offsets of the despinned ref frame
p34=p34-abs(coef(2,3))*cos(angle(coef(2,3))-phase-phi_34);
if use_asym==32
  p12 = p12-abs(coef(1,3))*cos(angle(coef(1,3))-phase-phi_32);
  p12 = 1.4142*p12 - p34;
  irf_log('proc','Using ASYMMETRIC probe 32 conf')
elseif use_asym==42
  p12 = p12-abs(coef(1,3))*cos(angle(coef(1,3))-phase-phi_42);
  p12 = 1.4142*p12 - p34;
  irf_log('proc','Using ASYMMETRIC probe 42 conf')
else, p12 = p12-abs(coef(1,3))*cos(angle(coef(1,3))-phase-phi_12);
end

% despin each component separately
dp12=p12.*exp(1i*phase)*exp(1i*phi_12);
dp34=p34.*exp(1i*phase)*exp(1i*phi_34);

% add necessary scaling
dp12=coef(1,1)*dp12;
dp34=coef(2,1)*dp34;

% create the final field
e=[t real(dp12+dp34) imag(dp12+dp34)];

switch ref_frame
  case 'wec'
    e(:,4)=es(:,2);
  case 'sat'
    e(:,4)=es(:,4);
end
