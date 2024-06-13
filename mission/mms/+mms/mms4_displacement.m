function [Rpar,Rperp,thetaBR] = mms4_displacement(varargin)
%MMS4_DISPLACEMENT Calculation displacements of MMS parallel and
%perpendicular to B
%
% Input:
%     'B?' - Magnetic fields for the four spacecraft.
%     'R?' - Positions of the four spacecraft.
%     'Tint' - Time interval (TSeries) to use if 'B?' and 'R?' are not
%     given.
%
% Output:
%     Rpar - TSeries of parallel displacements between spacecraft.
%     Rperp - TSeries of perpendicular displacements between spacecraft.
%     thetaBR - TSeries of angles between dispacements and B.
%
% Options:
%     plot - set to 1 to plot overview figure (otherwise no figure)
%
% Notes:
%     Rpar, Rperp, thetaBR orders are 12, 13, 14, 23, 24, 34.
%     If Tint is longer than 10 minutes survey data is used, otherwise
%     burst mode B is used if it spans Tint.
%
% Example:
%     [Rpar,Rperp,thetaBR] = mms.mms4_displacement('B?','R?','plot',1);
%     [Rpar,Rperp,thetaBR] = mms.mms4_displacement(Tint,'plot',0);

% if no arguments are given, display help
if (nargin < 1)
  help mms.mms4_displacement;
  Rpar = NaN;
  Rperp = NaN;
  thetaBR = NaN;
  return;
end

% setting default parameters
ic = 1:4; % MMS s/c numbers
getBR = false;
plotfig = false;
% varargin is cell array so index w/ {}
if isa(varargin{1},'EpochTT') % if Tint provided as 1st arg
  if length(varargin{1}) == 2 % interval = start and end time
    Tint = varargin{1};
    getBR = true; % Tint provided so need to find B and R
    argsstart = 2; % start args after Tint when called later
    irf.log('notice','Tint passed.') % logging that Tint has been passed
  end
else % if B and R provided, evaluate in base workspace
  c_eval('B?=evalin(''base'',irf_ssub(varargin{1},?));',ic); % ssub = string substitution
  c_eval('R?=evalin(''base'',irf_ssub(varargin{2},?));',ic);
  irf.log('notice','B and R are passed.')
  argsstart = 3;
end

if getBR % only true if Tint passed
  Tintl = Tint+[-60 60]; % assuming this adds 1min either side of interval
  R = mms.get_data('R_gse',Tintl); %#ok<NASGU> % get R_gse for all s/c
  c_eval('R? = irf.ts_vec_xyz(R.time,R.gseR?);',ic); % create TSeries obj
  if Tint(2)-Tint(1) > 600 % if time interval > 10 mins, use survey mode
    c_eval('B? = mms.get_data(''B_gse_fgm_srvy_l2'',Tint,?);',ic);
    irf.log('notice','Survey mode B is used.');
  else % if time interval <= 10 mins, use burst mode, unless Bint<Tint/1.2
    c_eval('B? = mms.get_data(''B_gse_fgm_brst_l2'',Tint,?);',ic);
    if B1.time.stop-B1.time.start < (Tint(2)-Tint(1))/1.2
      c_eval('B? = mms.get_data(''B_gse_fgm_srvy_l2'',Tint,?);',ic);
      irf.log('notice','Survey mode B is used.')
    else
      irf.log('notice','Burst mode B is used.')
    end
  end
end

% check for options args, i.e. 'plot' case
args=varargin(argsstart:end);
if numel(args)>0
  haveoptions=1;
else
  haveoptions=0;
end

% if plotting, check for 0 (no plot) or 1 (plot figure)
while haveoptions
  l = 2;
  switch(lower(args{1}))
    case 'plot'
      if numel(args)>1 && isnumeric(args{2})
        if args{2} > 0
          plotfig = true;
        end
      end
    otherwise
      irf.log('warning',['Unknown flag:' args{1}])
      l=1;
      break
  end
  args = args(l+1:end);
  if isempty(args), haveoptions=0; end
end
TSeries
% resample the B and R data at the B1 data sampling rate
c_eval('B? = B?.resample(B1);',ic);
c_eval('R? = R?.resample(B1);',ic);

% relative displacements between s/c
R_12=R1-R2; %#ok<NASGU>
R_13=R1-R3; %#ok<NASGU>
R_14=R1-R4; %#ok<NASGU>
R_23=R2-R3; %#ok<NASGU>
R_24=R2-R4; %#ok<NASGU>
R_34=R3-R4; %#ok<NASGU>

% relative B-field between s/c
B_12=(B1+B2)/2; %#ok<NASGU>
B_13=(B1+B3)/2; %#ok<NASGU>
B_14=(B1+B4)/2; %#ok<NASGU>
B_23=(B2+B3)/2; %#ok<NASGU>
B_24=(B2+B4)/2; %#ok<NASGU>
B_34=(B3+B4)/2; %#ok<NASGU>

c_eval('absR? = median(R_?.abs.data);',[12 13 14 23 24 34]);
c_eval('Rnorm_? = R_?/R_?.abs;',[12 13 14 23 24 34]);
c_eval('Bnorm_? = B_?/B_?.abs;',[12 13 14 23 24 34]);
c_eval('theta_BR_? = dot(Bnorm_?,Rnorm_?);',[12 13 14 23 24 34]);
c_eval('theta_BR_?.data = acosd(theta_BR_?.data);',[12 13 14 23 24 34]);
thetaBR = irf.ts_scalar(theta_BR_12.time,[theta_BR_12.data theta_BR_13.data ...
  theta_BR_14.data theta_BR_23.data theta_BR_24.data theta_BR_34.data]);

c_eval('[R_par_?,R_perp_?]=irf_dec_parperp(B_?,R_?);',[12 13 14 23 24 34]);
c_eval('R_perp_? = R_perp_?.abs;',[12 13 14 23 24 34]);

Rpar = irf.ts_scalar(R_par_12.time,[R_par_12.data R_par_13.data R_par_14.data ...
  R_par_23.data R_par_24.data R_par_34.data]);
Rperp = irf.ts_scalar(R_perp_12.time,[R_perp_12.data R_perp_13.data R_perp_14.data ...
  R_perp_23.data R_perp_24.data R_perp_34.data]);

if plotfig
  h = irf_plot(3,'newfigure');
  xSize=750; ySize=700;
  set(gcf,'Position',[10 10 xSize ySize]);

  xwidth = 0.88;
  ywidth = 0.30;
  set(h(1),'position',[0.10 0.95-ywidth xwidth ywidth]);
  set(h(2),'position',[0.10 0.95-2*ywidth xwidth ywidth]);
  set(h(3),'position',[0.10 0.95-3*ywidth xwidth ywidth]);

  h(1)=irf_panel('PosBpar');
  irf_plot(h(1),Rpar);
  irf_legend(h(1),{'12','13','14','23','24','34'},[1.0 1.0])
  ylabel(h(1),{'R_{||} (km)'},'Interpreter','tex');
  title(h(1),['|R_{12}| = ' num2str(round(absR12,1,'decimals')) ...
    ', |R_{13}| = ' num2str(round(absR13,1,'decimals')) ...
    ', |R_{14}| = ' num2str(round(absR14,1,'decimals')) ...
    ', |R_{23}| = ' num2str(round(absR23,1,'decimals')) ...
    ', |R_{24}| = ' num2str(round(absR24,1,'decimals')) ...
    ', |R_{34}| = ' num2str(round(absR34,1,'decimals')) ' km']);

  h(2)=irf_panel('PosBperp');
  irf_plot(h(2),Rperp);
  ylabel(h(2),{'|R_{\perp}| (km)'},'Interpreter','tex');

  h(3)=irf_panel('ThetaBR');
  irf_plot(h(3),thetaBR);
  ylabel(h(3),{'\theta_{BR} (^{o})'},'Interpreter','tex');

  Tint = irf.tint(B1.time.start.utc,B1.time.stop.utc);
  irf_plot_axis_align(1,h(1:3))
  irf_zoom(h(1:3),'x',Tint);
  xtickangle(h,0)
end

end

