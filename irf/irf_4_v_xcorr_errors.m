function [resspec] = irf_4_v_xcorr_errors(varargin)
% IRF_4_V_XCORR_ERRORS Four-spacecraft timing analysis with error estimates
%
%   SOMEWHAT EXPERIMENTAL
%
%   Based on:
%   Vogt, J., Haaland, S., and Paschmann, G.:
%   Accuracy of multi-point boundary crossing time analysis,
%   Ann. Geophys., 29, 2239-2252,
%   https://doi.org/10.5194/angeo-29-2239-2011, 2011.
%
%   resspec = IRF_4_V_XCORR_ERRORS(tint,S1,S2,S3,S4,R1,R2,R3,R4) calculates
%   time differences between the signals S1,...S4, and discontinuity
%   velocity given spacecraft position R1,...R4. The analysis is done in
%   the time interval tint. S1,...S4 are scalar TSeries objects.
%   R1,...R4 are TSeries objects. Returns a structure with the results of
%   the timing analysis.
%
%   IRF_4_V_XCORR_ERRORS(tint,'S?','R?') reads signals and
%   position from workspace.
%
%   IRF_4_V_XCORR_ERRORS(...,'Opt1',OptVal1,...) set options in function
%
%
%   resspec contains:
%       tau     -   time difference between spacecraft signal and pattern,
%                   1x4 array, in [s]
%       dtau2   -   <dta dtb> matrix, see eq (53), in [s^2], to get an
%                   estimate of average time uncertainties: try
%                   dt = diag(sqrt(dtau2))
%       V       -   discontinuity speed in [km/s]
%       n       -   normal unit vector vector
%       m       -   slowness vector, n/V, in [s/km]
%       dV      -   root-mean-square error of V in [km/s]
%       dTh     -   errors in propagation angle (semi-major/minor axis of
%                   error ellipse) in [radians]
%       eTh     -   corresponding vectors for min and max dTh, eTh(:,1)
%                   corresponds to min error and eTh(:,2) corresponds to
%                   max error
%       dm2     -   error matrix for m, dm2 = <dm*dm'> in [s^2/km^2]
%                   -> dm = diag(sqrt(dm2))
%
%
%   Options:
%   'p'     -   manually input TSeries pattern which to fit signals
%   'tau'   -   starting guess for time delays, 1x4 array, in [s]. For
%               example, run irf_4_v_gui first and guess the results from
%               there
%   'calcerrors'    -  boolean value whether to calculate errors, default
%                   is 1, set to 0 for faster calculation
%   'showfigure'    -  boolean value whether to show plots of the
%                   calculation, good for debugging, default is 0.
%                   'calcerrors' must be 1 otherwise it will crash.
%
%   Examples:
%       % interplanetary shock example
%       tint = irf.tint('2018-01-08T06:41:10/2018-01-08T06:41:11.5');
%       % shorter time interval
%       tint2 = irf.tint('2018-01-08T06:41:10.5/2018-01-08T06:41:11');
%       c_eval('B? = mms.get_data(''B_gse_fgm_brst_l2'',tint,?);')
%       R = mms.get_data('R_gse',tint+[-60,60]);
%       c_eval('R? = irf.ts_vec_xyz(R.time,R.gseR?);')
%       tspec = irf_4_v_xcorr_errors(tint2,'B?.z','R?','showfigure',1)
%
%   Some notes:
%       - The function finds local minima of the mean square deviation to
%       get the time delays. Therefore a good starting guess might be
%       necessary.
%       - If no pattern is in input, the function constructs a pattern. It
%       does this by using S1 as a pattern and fitting S2-S4 to it. Then
%       it constructs the new pattern as the average signal S1-S4
%       time-shifted to S1. Then it fits all signals to the new pattern.
%       - Make sure there's data on both sides of of the time interval used
%       for the timing. This is because the function shifts signals to
%       calculate errors.
%       - Function assumes spacecraft formation does not change during the
%       crossing time, which means that the crossing time should be short,
%       which it always is for MMS. It also means that the speed is
%       calculated in the spacecraft frame and not the Earth frame.
%       - The function assumes that the relative errors in the spacecraft
%       positions are small, which should hold for both MMS and Cluster.
%       - Negative or close-to-zero values in the off-diagonal elements of
%       dtau2 is generally a good sign. It means that the residuals are not
%       correlated.
%       - To get errors in speed and propagation angle, the function
%       assumes a scalar time uncertainty dt = trace(sqrt(dtau2))/4. This
%       could possibly be improved.
%       - The function is typically rather slow. Performance is improved in
%       case of shorter time intervals and lower resolution data.
%
%   See also: IRF_4_V_GUI, C_4_V_XCORR

%   Written by: Andreas Johlander, andreasj@irfu.se


%% Input
args = varargin;

% tint, S?, and R? must always be inputted

% tint must be first input
tint = args{1};

% always remove input that has been handled
args = args(2:end);

% get signals
if ischar(args{1}) && ~isempty(strfind(args{1},'?')) % like 'S?'
  varstr = args{1};
  c_eval('S?=evalin(''base'',irf_ssub(varstr,?));');
  args = args(2:end);
elseif isa(args{1},'TSeries') % All Ss are Tseries objects
  c_eval('S? = args{?}');
  args = args(5:end);
else
  error('Unkown signal input')
end

% get sc positions
if ischar(args{1}) && ~isempty(strfind(args{1},'?')) % like 'R?'
  varstr = args{1};
  c_eval('R?=evalin(''base'',irf_ssub(varstr,?));');
  args = args(2:end);
elseif isa(args{1},'TSeries') % All Rs are Tseries objects
  c_eval('R? = args{?}');
  args = args(5:end);
else
  error('Unkown position input')
end

% Handle the rest of inputs
nargs = length(args);

patternInput = 0;
startTau = -0.01*[1,1,1,1]; % starting guess for time differences
calcErrors = 1; % if error estimates should be calculated
showFigure = 0;

have_options = nargs > 1;
while have_options
  switch(lower(args{1}))
    case 'p' % pattern
      P = args{2};
      patternInput = 1;
    case 'tau' % times
      % minus sign because I don't know
      startTau = -args{2};
    case 'calcerrors' % times
      calcErrors = args{2};
    case 'showfigure' % times
      showFigure = args{2};
  end
  args = args(3:end);
  if isempty(args), break, end
end


%% Prepare analysis
disp('Performing Timing analysis:')
disp('---------------------------')

% get sample time (assume constant)
Tsamp = median(diff(S1.time.epochUnix));

%% If no pattern in input -> fit to S1 and then construct average pattern
if ~patternInput
  % Get time differences with first sc as reference
  % The goal is to minimize eq (45) for each sc pair (mean square deviation)
  % get best time differences (for some reason start guess must be negative)
  disp('Getting time differences with first sc as reference...')
  tau12 = fminsearch(@(x) msdTS(x,S2,S1,tint),startTau(2));
  tau13 = fminsearch(@(x) msdTS(x,S3,S1,tint),startTau(3));
  tau14 = fminsearch(@(x) msdTS(x,S4,S1,tint),startTau(4));

  % Construct an average signal to use as pattern
  disp('Constructing an average signal to use as new pattern...')
  c_eval('S?p = S?; S?p.time = S?p.time+tau1?;',2:4)
  % new pattern
  P = 1/4*(S1+S2p.resample(S1)+S3p.resample(S1)+S4p.resample(S1));
  % only in tint
  P = P.tlim(tint);

  %redefine startTau, negative starting guess for good measure
  startTau = [-1e-9,tau12,tau13,tau14];

end

%% Get time differences using the pattern
P = P.tlim(tint); % P defines the time interval

disp('Getting time differences using pattern...')
tau1 = fminsearch(@(x) msdTS(x,S1,P,tint),startTau(1));
tau2 = fminsearch(@(x) msdTS(x,S2,P,tint),startTau(2));
tau3 = fminsearch(@(x) msdTS(x,S3,P,tint),startTau(3));
tau4 = fminsearch(@(x) msdTS(x,S4,P,tint),startTau(4));

% the minus is empirical
tau = -[tau1,tau2,tau3,tau4];

%% Get time difference errors
% The goal is to calculate eq (52)
if calcErrors
  % the time interval used for averaging
  Tw = diff(P.time([1,end]).epochUnix)/2;

  % Define a time diff vector Delta
  Delta = P.time.epochUnix-mean(P.time.epochUnix);
  % number of points
  N = length(Delta);

  % get G (49)
  disp('Calculating G...')
  G = zeros(1,N);
  for ii = 1:N; G(ii) = corrG(Delta(ii),Tw,Tsamp,P); end

  % Get H (50)
  c_eval('H?! = zeros(1,N);',1:4,1:4); % initiate arrays
  % 16 times! wow!
  disp('Calculating H...')
  fprintf(1,'%d%%',0) % display progress
  for ii = 1:N; [H11(ii),h1] = corrH(Delta(ii),tau1,tau1,S1,S1,P); end
  for ii = 1:N; [H12(ii),h2] = corrH(Delta(ii),tau1,tau2,S1,S2,P); end
  for ii = 1:N; [H13(ii),h3] = corrH(Delta(ii),tau1,tau3,S1,S3,P); end
  for ii = 1:N; [H14(ii),h4] = corrH(Delta(ii),tau1,tau4,S1,S4,P); end
  fprintf(1,'\b\b%d%%',25);
  for ii = 1:N; H21(ii) = corrH(Delta(ii),tau2,tau1,S2,S1,P); end
  for ii = 1:N; H31(ii) = corrH(Delta(ii),tau3,tau1,S3,S1,P); end
  for ii = 1:N; H41(ii) = corrH(Delta(ii),tau4,tau1,S4,S1,P); end
  for ii = 1:N; H22(ii) = corrH(Delta(ii),tau2,tau2,S2,S2,P); end
  fprintf(1,'\b\b\b%d%%',50);
  for ii = 1:N; H23(ii) = corrH(Delta(ii),tau2,tau3,S2,S3,P); end
  for ii = 1:N; H24(ii) = corrH(Delta(ii),tau2,tau4,S2,S4,P); end
  for ii = 1:N; H32(ii) = corrH(Delta(ii),tau3,tau2,S3,S2,P); end
  for ii = 1:N; H42(ii) = corrH(Delta(ii),tau4,tau2,S4,S2,P); end
  fprintf(1,'\b\b\b%d%%',75);
  for ii = 1:N; H33(ii) = corrH(Delta(ii),tau3,tau3,S3,S3,P); end
  for ii = 1:N; H34(ii) = corrH(Delta(ii),tau3,tau4,S3,S4,P); end
  for ii = 1:N; H43(ii) = corrH(Delta(ii),tau4,tau3,S4,S3,P); end
  for ii = 1:N; H44(ii) = corrH(Delta(ii),tau4,tau4,S4,S4,P); end
  fprintf(1,'\b\b\b%d%%\n',100);



  Imin1 = msdTS(tau1,S1,P,tint);
  Imin2 = msdTS(tau2,S2,P,tint);
  Imin3 = msdTS(tau3,S3,P,tint);
  Imin4 = msdTS(tau4,S4,P,tint);
  dtau2 = zeros(4,4);

  % derivative of P
  Pp = P;
  Pp.data = [0;diff(P.data)/Tsamp]; % zero fill value

  % error through eq (52)
  c_eval('dtau2(?,!) = sqrt(Imin?*Imin!)/(N*mean(Pp.data.^2))*sum(G.*H?!);',1:4,1:4)
else
  dtau2 = 0;
end

%% Get propagation velocity from time lags
disp('Getting discontinuity velocity with corresponding errors...')

% time in center of interval
meanT = tint(1)+diff(tint.epochUnix)/2;

% call the function to get velocity with corresponding errors
[V,n,m,rmsV,dTh,eTh,dm2] = get_disc_vel(meanT,tau,dtau2,R1,R2,R3,R4);

%% Set output
resspec = [];

resspec.tau = tau;
resspec.dtau2 = dtau2;
resspec.V = V;
resspec.n = n;
resspec.m = m;
resspec.dV = rmsV;
resspec.dTh = dTh;
resspec.eTh = eTh;
resspec.dm2 = dm2;

%% Plot the results if requested
if showFigure
  disp('Plotting...')
  % set figure
  h = irf_plot(3,'newfigure');
  fn = h(1).Parent;
  set(fn,'color','white'); % white background for figures (default is grey)
  set(fn,'PaperUnits','centimeters')
  xSize = 10; ySize = 16;
  xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
  set(fn,'PaperPosition',[xLeft yTop xSize ySize])
  set(fn,'Position',[10 10 xSize*50 ySize*50])
  set(fn,'paperpositionmode','auto') % to get the same printing as on screen

  delete(h); % make new ones
  nPlots = 4+~patternInput;
  unitStr = S1.units;
  if isempty(unitStr); unitStr = 'a.u.'; end
  h = gobjects(1,nPlots);
  tintPlot = [S1.time(1);S1.time(end)];

  % original signals and possibly pattern
  hca = subplot(nPlots,1,1);
  irf_pl_tx(hca,'S?')
  hold(hca,'on');
  if patternInput; irf_plot(hca,P); end
  ylabel(hca,['$S$ [',unitStr,']'],'interpreter','latex')
  title(hca,'Original signal, and pattern if prescribed')
  irf_legend(hca,{'SC1','SC2','SC3','SC4'},[.02,.98],'color','cluster')
  irf_zoom(hca,'x',tintPlot)
  irf_pl_mark(hca,tint)
  xlabel(hca,'')
  h(1) = hca;

  % first time-shifted signals and new pattern
  if ~patternInput
    hca = subplot(nPlots,1,2);
    S1p = S1; % copy signal
    irf_pl_tx(hca,'S?p')
    hold(hca,'on')
    irf_plot(hca,P);
    ylabel(hca,['$S$ [',unitStr,']'],'interpreter','latex')
    title(hca,'Signals time-shifted to SC1 and constructed pattern')
    irf_legend(hca,{'SC1','SC2','SC3','SC4'},[.02,.98],'color','cluster')
    irf_zoom(hca,'x',tintPlot)
    irf_pl_mark(hca,tint)
    xlabel(hca,'')
    h(2) = hca;
  end

  % final time-shifted signals
  hca = subplot(nPlots,1,2+~patternInput);
  c_eval('S?pp = S?;'); c_eval('S?pp.time = S?.time+-tau(?);')
  irf_pl_tx(hca,'S?pp')
  hold(hca,'on')
  irf_plot(hca,P);
  ylabel(hca,['$S$ [',unitStr,']'],'interpreter','latex')
  title(hca,'Final time-shifted signals and pattern')
  irf_legend(hca,{'SC1','SC2','SC3','SC4'},[.02,.98],'color','cluster')
  irf_zoom(hca,'x',tintPlot)
  irf_pl_mark(hca,tint)
  xlabel(hca,'')
  h(2+~patternInput) = hca;

  % residuals
  hca = subplot(nPlots,1,3+~patternInput);
  % Plot like Fig 5a. in the paper
  plot(hca,Delta,h1,'Color','k')
  hold(hca,'on')
  plot(hca,Delta,h2,'Color','r')
  plot(hca,Delta,h3,'Color',[0,0.5,0])
  plot(hca,Delta,h4,'Color','b')
  title(hca,'Residuals')
  ylabel(hca,['$h_{\alpha}$ [',unitStr,']'],'interpreter','latex')
  xlabel(hca,'$\Delta$ [s]','interpreter','latex')
  irf_legend(hca,{'\alpha = 1','\alpha = 2','\alpha = 3','\alpha = 4'},[.02,.98],'color','cluster')
  hca.XLim = Delta([1,end]);
  h(3+~patternInput) = hca;

  % Product of correlation functions
  hca = subplot(nPlots,1,4+~patternInput);
  % Plot like Fig 5b. in the paper
  hold(hca,'on')
  title(hca,'Example of correlation functions')
  plot(hca,Delta,G.*H11,'Color','k')
  plot(hca,Delta,G.*H12,'Color','r')
  plot(hca,Delta,G.*H13,'Color',[0,0.5,0])
  plot(hca,Delta,G.*H14,'Color','b')
  ylabel(hca,'$G(\Delta)*H_{1,\beta}(\Delta)$','interpreter','latex')
  irf_legend(hca,{'\beta = 1','\beta = 2','\beta = 3','\beta = 4'},[.02,.98],'color','cluster')
  xlabel(hca,'$\Delta$ [s]','interpreter','latex')

  hca.XLim = Delta([1,end]);
  h(4+~patternInput) = hca;

  % axes loop
  for jj = 1:length(h)
    hca = h(jj);
    hca.Box = 'on';
    hca.LineWidth = 1.3;
    grid(hca,'on')
  end

end
disp('---------------------------')
end


% mean square deviation function
function [I,h] = msdTS(tau,S,P,tint)
% S is signal and P is pattern, here first sc is considered pattern

% time shift TS object
S.time = S.time+tau;

% start and stop times
tstart = tint(1);
tstop = tint(2);

S = S.resample(P);
% important to throw away data outside time interval
S.data(S.time<tstart | S.time>tstop) = 0;

% get residual h (48)
h = (S.data-P.data);

% mean square deviation I (45)
I = mean(abs(h).^2);
end


function G = corrG(Delta,Tw,Tsamp,P)
% get G in eq (49)
% p-prime is the time derivative! See appendix B

% start and stop times
% tstart = tint(1).epochUnix;
% tstop = tint(2).epochUnix;

% get time derivative Pprime(t)
Pp = P;
Pp.data = [0;diff(P.data)/Tsamp]; % zero fill value

% get G in (49)
Ppd = Pp; Ppd.time = P.time+Delta; % P(t+Delta)
Ppd = Ppd.resample(P);
% important to throw away data outside time interval
Ppd.data(Ppd.time<P.time(1)+Delta) = 0;
Ppd.data(Ppd.time>P.time(end)+Delta) = 0;

% average only over non-zero values
G = (1-abs(Delta)/Tw)*mean((Pp.data(Ppd.data~=0).*Ppd.data(Ppd.data~=0)))/(mean(Pp.data.^2));
end


function [Hab,hb] = corrH(Delta,tauStarA,tauStarB,Sa,Sb,P)
% get H in (50)

% shifted signals
Ssa = Sa;
Ssb = Sb;

Ssa.time = Sa.time+tauStarA;
Ssb.time = Ssb.time+tauStarB;

Ssa = Ssa.resample(P);
Ssb = Ssb.resample(P);

% residuals aafo time, eq (48)
ha = Ssa.data-P.data;
hb = Ssb.data-P.data;

hbTsDelta = irf.ts_scalar(P.time+Delta,hb);
hbTsDelta = hbTsDelta.resample(P);

% set to zero when curves are not overlapping (kind off a hack but should
% not matter)
hbTsDelta.data(hbTsDelta.time<P.time(1)+Delta) = 0;
hbTsDelta.data(hbTsDelta.time>P.time(end)+Delta) = 0;

hbDelta = hbTsDelta.data;

% average over non-zero values
Hab = mean(ha(hbDelta~=0).*hbDelta(hbDelta~=0))/(sqrt(mean(ha.^2))*sqrt(mean(hb.^2)));
end


%% The timing function
function [V,n,m,rmsV,dTh,eTh,dm2] = get_disc_vel(T,dTp,dT,R1,R2,R3,R4)
%GET_DISC_VEL Calculate velocity with errors with timing method.
%
% [V,n,m,rmsV,dTh,eTh,dm2] = GET_DISC_VEL(T,dTp,dT,R1,R2,R3,R4) Performs timing
%   analysis with:
%   Output:
%       V       - discontinuity speed
%       n       - normal vector
%       m       - slowness vector (n/V)
%       rmsV    - root-mean-square of V
%       dTh     - errors in propagation angle (semi-major/minor axis of
%               error ellipse)
%       eTh     - corresponding vectors for min and max dTh
%       dm2     - "squared" error for m, scalar or tensor, dm2 = <dm*dm'>
%               -> dm = diag(sqrt(dm2))
%   Input:
%       T       - epochTT object of crossing time
%       dTp     - time difference of discontinuity between spacecraft with
%               arbitrary time offset
%       dT      - error in time differences can be a scalar value or a
%               matrix of "squared" values <dt_a dt_b>
%       R1-R4   -  TSeries ocjects of spacecraft position.
%
% For rmsV, dTh, and eTh, a scalar time error is assumed even if inputted
% tensor. The operation for the scalar is:
% >> dT = trace(sqrt(dT_tensor))/4; % should be a real number
%
% This function follows the procedure and notation of section 3 & 4 in:
% Vogt, J., et al.,% Ann. Geophys., 29, 2239-2252, 2011.
%

% Rs are timeseries objects
% T is an epochTT object
% Time differences are doubles with unit seconds


% first check if time error is a scalar or tensor
if size(dT,2) == 4 % is a tensor
  dT_tensor = dT; % tensor value
  dT = trace(sqrt(dT_tensor))/4; % scalar value
elseif size(dT,1) == 1 % is a scalar
  dT_tensor = dT; % "tensor" value
else
  error('unkown format of dTp')
end

% tetrahedron center
rstar = 1/4*(R1.resample(T).data+R2.resample(T).data+...
  R3.resample(T).data+R4.resample(T).data)';

% position vectors with origin in tetrahedron center
c_eval('rstar? = R?.resample(T).data''-rstar;')

% get position tensor Rstar (1)
Rstar = rstar1*rstar1'+rstar2*rstar2'+rstar3*rstar3'+rstar4*rstar4';

% generalized reciprocal vectors q (3), (q=k)
c_eval('q? = inv(Rstar)*rstar?;')

% reciprocal tensor
K = q1*q1'+q2*q2'+q3*q3'+q4*q4';

% get slowness vector (17), which should hold for arbitrary time offsets
m = q1*dTp(1)+q2*dTp(2)+q3*dTp(3)+q4*dTp(4);

% normal vector from (9)
n = m/norm(m);

% discontinuity velocity from (9)
V = 1/norm(m);

% Get eigenvetors to R
% [eigV,eigD] = eig(Rstar);

% just get some normal vector perp to n
e1 = cross([0;1,;0],n)/norm(cross([0;1;0],n));
% and another
e2 = cross(n,e1);

% get errors in the timing
% Get absolute rms error in velocity (33) (deltaV/V = sqrt(<dV^2>)/V)
rmsV = V^2*dT*sqrt(n'*K*n);

% get error in angle (32)
% not an elegant solution, should be able to predict max and min errors and
% corresponding directions
% get many linear combinations of e1 and e2 (all are perp to n and norm=1)
eArr = e1*cos(linspace(0,2*pi,1e4))+e2*sin(linspace(0,2*pi,1e4));
% loop through all linear combinations and find max and min errors
dThArr = zeros(1,length(eArr));
for ii = 1:length(eArr)
  dThArr(ii) = V*dT*sqrt(eArr(:,ii)'*K*eArr(:,ii));
end
dTh = [min(dThArr),max(dThArr)];
eTh = [eArr(:,find(dThArr==min(dThArr),1)),eArr(:,find(dThArr==max(dThArr),1))];

% Get error in m
if size(dT_tensor,2) == 4 % tensor
  % from eq (21) and assumes dk (= dq) = 0
  dm2 = zeros(3);
  c_eval('dm2 = dm2+dT_tensor(?,!)*q?*q!'';',1:4,1:4)
else
  % from eq (27) and assumes dr = 0
  dm2 = dT^2*K;
end

end

