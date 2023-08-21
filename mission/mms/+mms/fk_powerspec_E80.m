function out = fk_powerspec_E80(varargin)
%
% out = fk_powerspec_E80(SCpot,T1)
%
% This function calculates the power spectral density in
% frequency-wavenumber space of the electric field in 2 orthogonal
% directions in the spin plane. Electric field from probes 2-3 and that
% from probes 4-1 are timed to get wavevector in the direction 23-->41 (x
% direction in this coordinate system) same thing is done for E field from
% probes 4-2 and 1-3 to get wavevector in the direction 42-->13 (y
% direction in this coordinate system).
%
% This function is based on mms.fk_powerspectrum written by Daniel B.
% Graham.
%
% This function is written by Ahmad Lalti.
%
% Inputs:
%
% SCpot - TSeries containing probe potentials (Note that time stamp
% correction are applied to SCpot in this function, so do not apply them
% before running the function.
%
% T1 - Time interval of the waveburst of interest
%
% Output:
%
% out - a cell array of 2 structures containing: pow_x/y the power, k_x/y
% wavenumber and f the frequency for both directions respectively.
%
% Options:
% boom_shortening - default is 0 where the boom shortening effect (Pederson
% et. al. 1998) is not taken into consideration. 1 take into consideration
% the boom shortening effect.
%
% w0 - Wavelet width of the wavelet transform, the default is w0 =5.36 the
% Morlet width.
%
% numf - number of elements in the frequency vector
%
% numk - number of elements in the wavenumber vector
%
% f - frequency range for the wavelet transform
%
% return_fields - flag to return electric fields in E80 and E120 coordinate
% system along with all 6 probe potentials
%
% correct_timeshifts - default is 1 where the time lag between probes is
% accounted for. 0 do not account for time lag.
%
% If no output is given the function plots the results
%
% Examples:
% out = fk_powerspec_E80(SCpot,T1)
% out = fk_powerspec_E80(SCpot,T1,'boom_shortening',1,'w0',4*5.36,'numf',200,'numk',200,'f',[100 4000])



if numel(varargin)<2
  help mms.fk_powerspec_E80
  return
end


SCpot = varargin{1};
T1 = varargin{2};

boom_shortening = 0;
w0 = 5.36;
numf = 100;
numk = 100;
f_range_flag = 0;
return_fields = 0;
correct_timeshifts = 1;

%% defining different options
if numel(varargin)>2

  in = varargin(3:end);
  flag = ~isempty(in);
  if mod(length(in),2) ~= 0
    disp('Check that each flag has a corresponding value')
    help mms.fk_powerspec_E80
    out=[];
    return
  end
else
  flag = 0;
end
while flag

  switch in{1}
    case 'boom_shortening'

      boom_shortening = in{2};
      in(1:2) = [];
      flag = ~isempty(in);

    case 'w0'

      w0 = in{2};
      in(1:2) = [];
      flag = ~isempty(in);

    case 'numf'

      numf = in{2};
      in(1:2) = [];
      flag = ~isempty(in);

    case 'numk'

      numk = in{2};
      in(1:2) = [];
      flag = ~isempty(in);
    case 'f'
      f_range_flag = 1;
      frange = in{2};
      in(1:2) = [];
      flag = ~isempty(in);
    case 'return_fields'
      return_fields = in{2};
      in(1:2) = [];
      flag = ~isempty(in);
    case 'correct_timeshifts'
      correct_timeshifts = in{2};
      in(1:2) = [];
      flag = ~isempty(in);
    otherwise

      irf.log('warning',['Unknown flag: ' in{1}]);
      return

  end


end

if correct_timeshifts
  %% Correct for timing in spacecraft potential data.
  E12 = TSeries(SCpot.time,(SCpot.data(:,1)-SCpot.data(:,2))/0.120);
  E34 = TSeries(SCpot.time,(SCpot.data(:,3)-SCpot.data(:,4))/0.120);
  E56 = TSeries(SCpot.time,(SCpot.data(:,5)-SCpot.data(:,6))/0.02815);
  V1 = TSeries(SCpot.time,SCpot.data(:,1));
  V3 = TSeries(SCpot.time+ 7.629e-6,SCpot.data(:,3));
  V5 = TSeries(SCpot.time+15.259e-6,SCpot.data(:,5));
  E12.time = E12.time + 26.703e-6;
  E34.time = E34.time + 30.518e-6;
  E56.time = E56.time + 34.332e-6;

  test_resamp=1;
  if test_resamp
    V3 = mms.dft_timeshift(V3,-7.629e-6);
    V5 = mms.dft_timeshift(V5,-15.259e-6);
    E12 = mms.dft_timeshift(E12,-26.703e-6);
    E34 = mms.dft_timeshift(E34,-30.518e-6);
    E56 = mms.dft_timeshift(E56,-34.332e-6);
  else
    V3 = V3.resample(V1.time); %#ok<UNRCH>
    V5 = V5.resample(V1.time);
    E12 = E12.resample(V1.time); %These resamples need to be changed.
    E34 = E34.resample(V1.time);
    E56 = E56.resample(V1.time);

  end

  V2 = V1 - E12 * 0.120;
  V4 = V3 - E34 * 0.120;
  V6 = V5 - E56 * 0.0292;
  % Make new SCpot with corrections
  SCpot = irf.ts_scalar(V1.time,[V1.data V2.data V3.data V4.data V5.data V6.data]);
  clear V1 V2 V3 V4 V5 V6 E12 E34 E56
end

%% defining Tseries and focusing on the waveburst

V=SCpot.data;Time=SCpot.time;
V1=irf.ts_scalar(Time,V(:,1));V2=irf.ts_scalar(Time,V(:,2));
V3=irf.ts_scalar(Time,V(:,3));V4=irf.ts_scalar(Time,V(:,4));
V5=irf.ts_scalar(Time,V(:,5));V6=irf.ts_scalar(Time,V(:,6));
V1=V1.tlim(T1);V2=V2.tlim(T1);V3=V3.tlim(T1);V4=V4.tlim(T1);V5=V5.tlim(T1);V6=V6.tlim(T1);
V1.name='V1';V2.name='V2';V3.name='V3';V4.name='V4';V5.name='V5';V6.name='V6';

%% define baseline

if boom_shortening
  dl56 = 14.5;dl12 = 96;dl13 = sqrt(2)*dl12/2;
else
  dl56 = 28.15;dl12 = 120;dl13 = sqrt(2)*dl12/2;
end



%% calculate Efields
E56=-(V5-V6)/dl56;

E13=-(V3-V1)/dl13;E14=-(V4-V1)/dl13;
E80=[E13.data,E14.data,E56.data]*1000;E80=irf.ts_vec_xyz(V1.time,E80);
E80.coordinateSystem = 'E80';

E23=-(V3-V2)/dl13;E24=-(V4-V2)/dl13;
E80_2=[-E24.data,-E23.data,E56.data]*1000;E80_2=irf.ts_vec_xyz(V1.time,E80_2);%%the negative so the timing is between the same vectors (13 and 42, 14 and 32)
E80_2.coordinateSystem = 'E80_2';

%%%calculate Efield with probes 1 2 3 and 4
E12=-(V2-V1)/dl12;E34=-(V4-V3)/dl12;
E120=[-E12.data,-E34.data,E56.data]*1000;E120=irf.ts_vec_xyz(V1.time,E120);%%the negative sign is to have the positive x in the direction of probe 1
E120.coordinateSystem = 'E120';

%% calculating wavelet transforms

fmin=E80.time(2)-E80.time(1);fmin=1/fmin;fmin=fmin/length(E80.data);fmin=ceil(fmin);
if ~ f_range_flag
  frange = [fmin 4000];
end

w80=irf_wavelet(E80,'nf',numf,'returnpower',0,'f',frange,'wavelet_width',w0);
w80_2=irf_wavelet(E80_2,'nf',numf,'returnpower',0,'f',frange,'wavelet_width',w0);
f=w80.f;

%% calculating powerspectrum in the f,k_x space

p80x=w80.p{2};p80x_2=w80_2.p{2};
[k_80_x,pow_80_x]=disprel(numf,numk,p80x,p80x_2,dl13);

%% calculating powerspectrum in the f,k_y space

p80y=w80.p{1};p80y_2=w80_2.p{1};
[k_80_y,pow_80_y]=disprel(numf,numk,p80y,p80y_2,dl13);

out{1} = struct('pow_x',pow_80_x,'k_x',k_80_x,'f',f);
out{2} = struct('pow_y',pow_80_y,'k_y',k_80_y,'f',f);
TS = struct('E80',E80,'E80_2',E80_2,'E120',E120,'V1',V1,'V2',V2,'V3',V3,'V4',V4,'V5',V5,'V6',V6);
if return_fields
  out{length(out)+1} = TS;
end

if nargout == 0

  figure;
  h80x=subplot('Position',[0.07 0.25 0.4 0.4]);
  h80y=subplot('Position',[0.55 0.25 0.4 0.4]);
  surf(h80x,k_80_x,f,log10(pow_80_x)')
  shading(h80x,'interp')
  view(h80x,0,90)
  axis(h80x,'tight')
  xlabel(h80x,'k 1/m')
  ylabel(h80x,'f Hz')
  title(h80x,'E80x')
  colorbar(h80x)
  c=caxis(h80x);

  surf(h80y,k_80_y,f,log10(pow_80_y)')
  shading(h80y,'interp')
  view(h80y,0,90)
  axis(h80y,'tight')
  xlabel(h80y,'k 1/m')
  ylabel(h80y,'')
  yticklabels(h80y,'')
  title(h80y,'E80y')
  colorbar(h80y)
  caxis(h80y,c)
  colormap jet

end
end


function [kvec,pow]=disprel(numf,numk,p1,p2,dl)

%%p1 should be the power in the positive direction of the coordinate system


phi12=p2.*conj(p1);
phi12=atan2(imag(phi12),real(phi12));
k=phi12/(dl);
Powerav=(p1.*conj(p1)+p2.*conj(p2))/2;


N=size(Powerav,1);
pow = zeros(numk+1,numf);
mink=min(min(k));maxk=max(max(k));

dk = (maxk - mink)/numk;
kvec = mink + [0:1:numk]*dk;

for m = 1:N
  for q = 1:numf

    knumber = floor((k(m,q)-mink)/dk)+1;

    if ~isnan(knumber)
      pow(knumber,q) = pow(knumber,q) + Powerav(m,q);
    end


  end
end
end