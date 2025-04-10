function [outp] = irf_shock_rh(inp,varargin)
%IRF_SHOCK_RH Solves the Rankine-Hugoniot relations for fast mode shocks.
%
% outParam = IRF_SHOCK_RH(upParam) Returns structure with mostly
% downstream parameters. upParams is a structure with fields specifying
% input.
%
% IRF_SHOCK_RH(flag1,value1,flag2,value2,...) Set input with flags and
% values instead of a structure.
%
% IRF_SHOCK_RH('gui') Creates a gui for Rankine-Hugoniot testing
%
% Dimensonless input (fields or flags):
%   Mf    - Fast Mach number in the NI frame (Mf>=1)
%   Ma    - Alfven Mach number in the NI frame
%   Ms    - Sonic Mach number in the NI frame
%             NB: only one Mach number in input is required, Mf
%             overwrites Ma which overwrites Ms
%   th    - Shock normal angle in degrees
%   beta  - Plasma beta
%  optional:
%   gamma - Heat capacity ratio (5/3 if omitted)
%   rr    - Array of compression ratios to test (only for debugging)
%   out   - Cell of strings or single string specifying outputs,
%             e.g 'nd' for only compression ratio, default is a large
%             number of outputs
%
% Output if dimensionless input:
%   Bd    - Downstream magnetic field in units of norm(Bu)
%   Vd    - Downstream flow velocity in units of upstream Alfven
%             speed
%   Vd_Vu - Downstream flow velocity in units of Vu
%   Vu    - Upstream velocity in units of Alfven speed
%   nd    - Downstream number density in units of nu
%   Pd    - Downstream thermal pressure Pd = nd*kB*Td
%   Pu    - Upstream thermal pressure Pu = nu*kB*Tu
%   thd   - Downstream shock angle in degrees
%   Bdt   - Downstream tangential B (same as Bd(2))
%   betad - Downstream plasma beta
%   Mdf   - Downstream fast Mach number
%   Mds   - Downstream sonic Mach number
%   Mda   - Downstream Alfven Mach number
%   MdI   - Downstream intermediate Mach number
%   Muf   - Upstream fast Mach number
%   Mus   - Upstream sonic Mach number
%   Mua   - Upstream Alfven Mach number
%   MuI   - Upstream intermediate Mach number
%   vda   - Downstream Alfven speed in units of upstream Alfven
%             speed
%   thdVn - Downstream angle between velocity and normal vector
%   nSol  - Number of solutions to the RH shock conditions, there are
%             normally two solutions (trivial and shock) but switch-on
%             shocks have four solutions
%   f4n   - Normalized energy difference |dW|/Wu given rr (for
%             debugging, not in default output)
%
% Input in physical units:
%   B - Magnetic field vector in [nT]
%   V - Plasma flow velocity vector in [km/s]
%   n - Plasma number density in [cm^-3]
%   T - Plasma temperature in [eV]
%   Ti/Te - Separate ion/electron temperatures (not recommended)
%   nvec - Normal vector of shock ([1,0,0] if omitted)
%
% Output in physical units:
%   Bd - Downstream magnetic field vector in [nT]
%   Vd - Downstream plasma flow velocity vector in [km/s]
%   nd - Downstream plasma number density in [cm^-3]
%   Td - Downstream plasma temperature in [eV]
%   Tid/Ted - Ion/electron temperatures in [eV] (assumes adiabatic
%               electron heating)
%   dimLessOut - All dimensionless outputs
%
%
% Notes:
%   - The function uses a semi-analytical approach where the difference in
%   energy flux between up- and downstream is calculated for a given
%   density compression ratio. The solution is then found where the energy
%   fluxes are the same.
%   - If Mf in the input is close to 1, the function might have
%   difficulties separating the shock solution from the trivial no-shock
%   solution
%   - Shock angle th close to 0 can be cause issues
%   - The function can deal with switch-on shocks but care should be taken
%   when this is the case
%
% Examples:
%   % Example 1: simple use
%   inp = []; inp.Mf = 3; inp.th = 60; inp.beta = .5;
%   outp = IRF_SHOCK_RH(inp)
%
%   % Example 2: calculate and plot compression ratio as a function of Ma
%   Ma = linspace(1.01,10,50); nd = zeros(1,50);
%   for ii = 1:50;
%     nd(ii) = IRF_SHOCK_RH('beta',0,'th',90,'Ma',Ma(ii),'out','nd');
%   end
%   figure; plot(Ma,nd); xlabel('M_A'); ylabel('n_d/n_u'); set(gca,'FontSize',18);
%
%   % Example 3: GUI
%   irf_shock_rh gui
%
%   % Example 4: First critical Mach number. Famous Figure 4c in (Edmiston &
%   % Kennell, 1984, 10.1017/S002237780000218X)
%   Mf = linspace(2.7,1.001,128);
%   beta = linspace(0,1,16); th = linspace(0.5,90,16);
%   % solve RH conditions a bunch of times
%   Mds = zeros(length(th),length(beta),length(Mf));
%   for ii = 1:length(th); disp([num2str(ii),'/',num2str(length(th))])
%     for jj = 1:length(beta)
%       for kk = 1:length(Mf)
%         Mds(ii,jj,kk) = irf_shock_rh('beta',beta(jj),'th',th(ii),'Mf',Mf(kk),'Np',1e3,'out','Mds');
%       end
%     end
%   end
%   % get critical Mach number
%   Mc = zeros(length(th),length(beta));
%   for ii = 1:length(th)
%     for jj = 1:length(beta)
%       idTemp = find(squeeze(Mds(ii,jj,:))>1,1,'first'); % find where Mds=1
%       Mc(ii,jj) = Mf(idTemp);
%     end
%   end
%   % plot
%   figure; contour(th,beta,Mc',1:.1:2.7,'k','linewidth',3)
%   xlabel('\theta_B_n'); ylabel('\beta');
%   set(gca,'XLim',[0,90]); set(gca,'YLim',[0,1]); set(gca,'FontSize',18);
%
%   % Example 5: physical units
%   inp = [];
%   inp.n = 3; % cm^-3
%   inp.B = [-1,1,5]; % nT
%   inp.V = [-425,30,0]; % km/s
%   inp.T = 20; % eV
%   inp.nvec = [1,-1,0]./sqrt(2);
%   outp = irf_shock_rh(inp)


% TODO: Probably a lot

% Some upstream parameters are:
%   Bu = [cosd(th),sind(th),0]
%   Vu = [Mf*vms(th),0,0]
%   nu = 1
%   va = 1
% constants:
%   m = 1
%   mu0 = 1
%   kB = 1


% Written by A. Johlander


%% Handle input

% check if input is a string
% if it is, it is either an option ("gui","plot"), or a flag
% if it's a flag, then all inputs are flag, value, flag, value,....
if isa(inp,'char')

  % options
  if strcmp(inp,'gui')
    % create gui
    % error('GUI not implemented yet!')
    ud = [];
    ud = def_outputs(ud);
    ud = init_gui(ud);
    irf_shock_rh('plot');
    return;
  elseif strcmp(inp,'plot')
    % only 'plot' now
    ud = get(gcf,'userdata');
    ud = set_par(ud);
    ud = plot_rankhug(ud);
    set(gcf,'userdata',ud);
    return;
  else % flag
    % redefine inp structure input
    if nargin<2; error('not enough inputs'); end
    % first input is special
    tempinp.(inp) = varargin{1};
    inp = tempinp;
    for ii = 2:2:nargin-1
      inp.(varargin{ii}) = varargin{ii+1};
    end

  end
end


%% input with units (i.e. spacecraft measurement)

% check if B is in input and th is not
if isfield(inp,'B') && ~isfield(inp,'th')
  dimInp = 1;
else
  dimInp = 0;
end


if dimInp
  % define dimensional constants
  u = irf_units;
  mp_Dim = u.mp;
  mu0_Dim = u.mu0;
  qe_Dim = u.e;

  % get paramters from input
  B = inp.B; % [nT]
  n = inp.n; % [cm^-3]
  V = inp.V; % [km/s]
  % if Ti and Te is in input, calculate Mf with yi=3 and ye=1
  if isfield(inp,'Ti') && isfield(inp,'Te')
    useKineticMf = 1;
    Ti = inp.Ti; % [eV]
    Te = inp.Te; % [eV]
  else % otherwise use MHD-like temperature and Mf
    useKineticMf = 0;
    T = inp.T; % [eV]
    % copy pasted from further down
    if isfield(inp,'gamma')
      gamma = inp.gamma;
    else
      % default value
      gamma = 5/3;
    end
  end

  % get coordinate system
  if ~isfield(inp,'nvec')
    inp.nvec = [1,0,0];
    irf.log('w','setting normal vector to x')
  end
  nvec = inp.nvec;

  t2vec = cross(nvec,B)/norm(cross(nvec,B));
  t1vec = cross(t2vec,nvec);

  % get normal incidence frame
  if ~isfield(inp,'Vsh')
    inp.Vsh = 0;
  end
  Vsh = inp.Vsh;

  vNIF = V-(dot(V,nvec-Vsh)*nvec);
  VnData = dot(V,nvec);

  % set dimensionless parameters to input
  bvec = B/norm(B);
  inp.th = acosd(dot(nvec,bvec));
  if inp.th>90; inp.th = 180-inp.th; end

  if useKineticMf
    % this might be slow
    inp.ref_sys = 'nif';
    shp = irf_shock_parameters(inp);
    inp.Mf = shp.Mf; % fast Mach
    inp.beta = shp.bi+shp.be; % beta
  else
    % this is wrong! (units mostly)
    inp.beta = 2*mu0_Dim*n*1e6*T*qe_Dim/norm(B*1e-9)^2;
    va = norm(B*1e-9)/sqrt(mu0_Dim*n*mp_Dim*1e6);
    cs = sqrt(gamma/2*inp.beta)*va; % sound speed
    cms = sqrt(va^2+cs^2);
    % magnetosonic speed
    th = acosd(dot(B/norm(B),nvec));
    vms = sqrt(cms^2/2+sqrt(cms^4/4-va^2*cs^2*cosd(th)^2)); % [m/s]
    inp.Mf = -VnData*1e3/vms; % fast Mach
  end

  % other goodems
  BnDataSign = sign(dot(B,nvec));
end


%% normal input

% internally the progam has the following assumptions
% Vn is positive (unintuitive)
% Bn is always positive
% Bt is positive
% Therefore Vdt is also positive

% Set parameters
% theta and beta should always be in input
th = inp.th;
beta = inp.beta;

% check for mach numbers
MfInput = 0;
MaInput = 0;
MsInput = 0;
if isfield(inp,'Mf')
  Mf = inp.Mf;
  MfInput = 1;
end
if isfield(inp,'Ma')
  Ma = inp.Ma;
  MaInput = 1;
end
if isfield(inp,'Ms')
  Ms = inp.Ms;
  MsInput = 1;
end

% Mf overwrites Ma overwrites Ms
if MfInput && MaInput
  irf.log('w','Both Mf and Ma in input, Mf is used')
  MaInput = 0;
end
if MfInput && MsInput
  irf.log('w','Both Mf and Ms in input, Mf is used')
  MsInput = 0;
end
if MaInput && MsInput
  irf.log('w','Both Ma and Ms in input, Ma is used')
  MsInput = 0;
end
if MfInput && MaInput && MsInput
  irf.log('w','Mf, Ma and Ms in input, Mf is used')
  MaInput = 0;
  MsInput = 0;
end

if isfield(inp,'gamma')
  gamma = inp.gamma;
else
  % default value
  gamma = 5/3;
end
% maximum compression ratio
rg = gamma/(gamma-1);

% default values
if ~isfield(inp,'out'); inp.out = {'full'}; end
if ~iscell(inp.out); inp.out = {inp.out}; end
if ~isfield(inp,'Np'); inp.Np = 1e3; end

% set output strucutre
outp = [];
% special
if strcmpi(inp.out{1},'full')
  inp.out = {'Bd','Vd','Vd_Vu','Vu','nd','Pd','Pu','thd','Bdt','betau','betad',...
    'Mdf','Mds','Mda','MdI','Muf','Mus','Mua','MuI','vda','thdVn','nSol'};
end


%% Set parameters

% Bu and nu are defined as 1
Bu = 1; %
nu = 1;

Bn = Bu*cosd(th);
But = Bu*sind(th);

% all velocities are normalized to the Alfven speed
va = 1; % means that mu0*m = 1. Let's set them both to 1.
mu0 = 1;
m = 1;

cs = sqrt(gamma/2*beta)*va; % sound speed
cms = sqrt(va^2+cs^2);
% magnetosonic speed
vms = sqrt(cms^2/2+sqrt(cms^4/4-va^2*cs^2*cosd(th)^2));

if MfInput
  Vu = Mf*vms;
elseif MaInput
  Vu = Ma*va;
  Mf = Vu/vms;
elseif MsInput
  Vu = Ms*cs;
  Mf = Vu/vms;
end

% return nans if no fast shock solution
if Mf<1
  % return nans
  for k = 1:length(inp.out)
    % could be the wrong size but who cares?
    eval(['outp.',inp.out{k},'=nan;'])
  end
  return;
end

% define a number of upstream parameters for output and/or calculations
% mach numbers
Muf = Mf;
Mua = Vu/va;
Mus = Vu/cs;
MuI = Mua*secd(th); % [Schwartz, 1998] (ISSI book)
% other
Pu = beta*Bu^2/2;
betau = beta;


%% Solve the Rankine-Hugoniot conditions
% The function "guesses" compression ratio (nd) and solves for the
% remaining four unknowns

% analytical solutions as functions of compression ratio
% density
nd = @(r)nu*r;
% normal speed
Vdn = @(r)Vu./r;
% tangential velocity
Vdt = @(r)(Bn.*But.*(Vu./Vdn(r)-1))./(nd(r).*m.*mu0.*Vdn(r)-Bn^2./Vdn(r));
% tangential B field
Bdt = @(r)(Vdt(r).*Bn+Vu*But)./Vdn(r);
% pressure
Pd = @(r)m*nu*Vu^2+Pu+(But^2)/(2*mu0)-m*nd(r).*Vdn(r).^2-(Bdt(r).^2)/(2*mu0);

% energy flux difference between up- and downstream
% Solve f4(r) = 0, numerically to find compression ratio
f4 = @(r)(m*nu*Vu*(Vu^2/2+rg*Pu/(m*nu)+But^2/(mu0*m*nu))...
  -m*nd(r).*Vdn(r).*((Vdn(r).^2+Vdt(r).^2)/2+gamma/(gamma-1)*Pd(r)./...
  (m*nd(r))+Bdt(r).^2./(mu0*m*nd(r)))+(Vdt(r).*Bdt(r))*Bn/mu0);



%% Solve the jump conditions - ACTUAL ENGINE PART OF FUNCTION

% find zero crossing
if isfield(inp,'rr') % user-defined grid
  rr = inp.rr;
else % automatic
  Mmin = 0.99;
  Mmax = 1.01*(gamma+1)/(gamma-1);
  rr = linspace(Mmin,Mmax,inp.Np);
end

f4v = f4(rr);
% find zero crossings
idz = find(f4v(1:end-1).*f4v(2:end)<0);

% if no zero crossings are found, there is no shock since Mf is close to 1
if isempty(idz)
  idz = [1,1];
  nSol = 0;
else
  % number of solutions found, normally it should be 2, but SO shocks
  % have 4 (?)
  nSol = length(idz);
end

% test if increasing
if f4(rr(idz(1)+1))-f4(rr(idz(1)))>0
  idz = idz(1);
else
  % try the next one
  idz = idz(2);
end
% interpolate to get a better estimate, sort of slow
% r0 = interp1([f4(rr(idz+1)),f4(rr(idz))],[rr(idz),rr(idz+1)],0);

% just using mean value is faster and seems to work just as well
r0 = (rr(idz)+rr(idz+1))/2;
r = r0;


%% define output

% output functions
Bd = @(r)[Bn,Bdt(r),0];
Vd = @(r)[-Vdn(r),Vdt(r),0];
Vd_Vu = @(r)[-Vdn(r),Vdt(r),0]/Vu;

thd = @(r)atand(Bdt(r)/Bn);
betad = @(r)Pd(r)*2*mu0/norm(Bd(r))^2;

thdVn = @(r)atand(Vdt(r)/Vdn(r));

vda = @(r)norm(Bd(r))/sqrt(nd(r)*m*mu0);
cds = @(r)sqrt(gamma*betad(r)/2)*vda(r);
vdf = @(r)(vda(r)^2+cds(r)^2)/2+sqrt((vda(r)^2+cds(r)^2)^2/4-vda(r)^2*cds(r)^2*cosd(thd(r)).^2);
vdt = @(r)sqrt(2*Pd(r)/(m*nd(r))); % ion thermal speed
Mdf = @(r)Vdn(r)/vdf(r);
Mda = @(r)Vdn(r)/vda(r);
Mds = @(r)Vdn(r)/cds(r);
MdI = @(r)Mda(r)*secd(thd(r)); % [Schwartz, 1998] (ISSI book)
Mdt = @(r)Vdn(r)/vdt(r);
f4n = @(r)f4(rr)./(m*nu*Vu*(Vu^2/2+rg*Pu/(m*nu)+But^2/(mu0*m*nu))); %normalized


% Outputs are either values or functions, they are called slightly
% diffrently
for k = 1:length(inp.out)
  if isnumeric(eval(inp.out{k}))
    eval(['outp.',inp.out{k},'=',inp.out{k},';'])
  else % assume it's function
    % call function with final r
    eval(['outp.',inp.out{k},'=',inp.out{k},'(r);'])
  end
end

% no struct if only one output
fn = fieldnames(outp);
if isscalar(fn); outp = outp.(fn{1}); end


%% If input is dimensional, output should be too

if dimInp
  % overwrite everything
  outpTemp = outp;
  outp = [];
  % add the dimensionless output to structure
  outp.dimLessOut = outpTemp;

  R = [nvec;t1vec;t2vec];
  Vdpp = -outpTemp.Vd_Vu*VnData;
  Vdp = R\Vdpp';
  outp.Vd = Vdp'+vNIF;
  Bdp = outpTemp.Bd*norm(inp.B);
  Bdp(1) = BnDataSign*Bdp(1);
  outp.Bd = (R\Bdp')';

  outp.nd = outpTemp.nd*inp.n;

  % temperatures
  % first assume adiabatic electron heating
  if useKineticMf
    Ted = Te*(outp.nd/nu)^(gamma-1);
    outp.Ted = Ted;
    % rest of energy goes to ions
    Tid = (outp.dimLessOut.betad*norm(outp.Bd*1e-9)^2/(2*mu0_Dim*outp.nd*1e6))/qe_Dim-Ted;
    outp.Tid = Tid;
  else % normal MHD
    Td = (outp.dimLessOut.betad*norm(outp.Bd*1e-9)^2/(2*mu0_Dim*outp.nd*1e6))/qe_Dim;
    outp.Td = Td;
  end
end

end


% -------------- GUI part of function --------------------
function ud = def_outputs(ud)

% list in pop-up menus
ud.out.names = {'nd/nu',...
  'Bd/Bu',...
  'Btd/Btu',...
  'thd',...
  'thd-th',...
  'Mdf',...
  'MdA',...
  'MdI',...
  'betad',...
  'Vd',...
  'Vdn',...
  'Vdt',...
  'th_dVn',...
  'Td/Tu',...
  'nSol'};

ud.out.labels = {'$n_d/n_u$',...
  '$B_d/B_u$',...
  '$B_{dt}/B_{tu}$',...
  '$\theta_{d,Bn}$',...
  '$\theta_{d,Bn}-\theta_{Bn}$',...
  '$M_{df}$',...
  'MdA',...
  'MdI',...
  '$\beta_d$',...
  '$V_d/V_{uA}$',...
  '$V_{dn}/V_{uA}$',...
  '$V_{dt}/V_{uA}$',...
  '$\theta_{d,Vn}$',...
  '$T_d/T_u$',...
  'Number of solutions'};
end


function ud = init_gui(ud)
% init figure
ud.fig = figure;
ud.ax = axes(ud.fig);
ud.ax.Position = [0.12,0.2,0.52,.7];

%% input panel
% panel
ud.ph1.panel = uipanel('Units', 'normalized',...
  'position',[0.69 0.62 0.29 0.37],...
  'fontsize',14,...
  'Title','Input');

% x text
ud.ph1.tx = uicontrol('style','text',...
  'Units', 'normalized',...
  'Parent',ud.ph1.panel,...
  'position',[0.05 0.8 0.3 0.15],...
  'fontsize',14,...
  'HorizontalAlignment','left',...
  'string','X: ');
% x pop-up menu
ud.ph1.px = uicontrol('style','popupmenu',...
  'Units', 'normalized',...
  'Parent',ud.ph1.panel,...
  'position',[0.3 0.8 0.65 0.15],...
  'fontsize',14,...
  'string',{'th','Mf','beta'},...
  'Value',1);

% y text
ud.ph1.ty = uicontrol('style','text',...
  'Units', 'normalized',...
  'Parent',ud.ph1.panel,...
  'position',[0.05 0.65 0.3 0.15],...
  'fontsize',14,...
  'HorizontalAlignment','left',...
  'string','Y: ');
% y pop-up menu
ud.ph1.py = uicontrol('style','popupmenu',...
  'Units', 'normalized',...
  'Parent',ud.ph1.panel,...
  'position',[0.3 0.65 0.65 0.15],...
  'fontsize',14,...
  'string',{'th','Mf','beta'},...
  'Value',2);

% th text
ud.ph1.tth = uicontrol('style','text',...
  'Units', 'normalized',...
  'Parent',ud.ph1.panel,...
  'position',[0.05 0.45 0.3 0.15],...
  'fontsize',14,...
  'HorizontalAlignment','left',...
  'string','th: ');
% th text box
ud.ph1.bth = uicontrol('style','edit',...
  'Units', 'normalized',...
  'Parent',ud.ph1.panel,...
  'position',[0.3 0.45 0.65 0.15],...
  'fontsize',12,...
  'string','45');

% Mf text
ud.ph1.tmf = uicontrol('style','text',...
  'Units', 'normalized',...
  'Parent',ud.ph1.panel,...
  'position',[0.05 0.25 0.3 0.15],...
  'fontsize',14,...
  'HorizontalAlignment','left',...
  'string','Mf: ');
% Mf text box
ud.ph1.bmf = uicontrol('style','edit',...
  'Units', 'normalized',...
  'Parent',ud.ph1.panel,...
  'position',[0.3 0.25 0.65 0.15],...
  'fontsize',12,...
  'string','3');

% beta text
ud.ph1.tb = uicontrol('style','text',...
  'Units', 'normalized',...
  'Parent',ud.ph1.panel,...
  'position',[0.05 0.05 0.3 0.15],...
  'fontsize',14,...
  'HorizontalAlignment','left',...
  'string','beta: ');
%  beta text box
ud.ph1.bb = uicontrol('style','edit',...
  'Units', 'normalized',...
  'Parent',ud.ph1.panel,...
  'position',[0.3 0.05 0.65 0.15],...
  'fontsize',12,...
  'string','0');

%% interval panel
% panel
ud.ph2.panel = uipanel('Units', 'normalized',...
  'position',[0.69 0.4 0.29 0.22],...
  'fontsize',14,...
  'Title','Intervals');

% th text
ud.ph2.tth = uicontrol('style','text',...
  'Units', 'normalized',...
  'Parent',ud.ph2.panel,...
  'position',[0.05 0.8 0.3 0.2],...
  'fontsize',14,...
  'HorizontalAlignment','left',...
  'string','th: ');
%  th1 text box
ud.ph2.bth1 = uicontrol('style','edit',...
  'Units', 'normalized',...
  'Parent',ud.ph2.panel,...
  'position',[0.3 0.8 0.3 0.2],...
  'fontsize',12,...
  'string','1');
%  th2 text box
ud.ph2.bth2 = uicontrol('style','edit',...
  'Units', 'normalized',...
  'Parent',ud.ph2.panel,...
  'position',[0.65 0.8 0.3 0.2],...
  'fontsize',12,...
  'string','90');


% Mf text
ud.ph2.tmf = uicontrol('style','text',...
  'Units', 'normalized',...
  'Parent',ud.ph2.panel,...
  'position',[0.05 0.45 0.3 0.2],...
  'fontsize',14,...
  'HorizontalAlignment','left',...
  'string','Mf: ');
%  Mf1 text box
ud.ph2.bmf1 = uicontrol('style','edit',...
  'Units', 'normalized',...
  'Parent',ud.ph2.panel,...
  'position',[0.3 0.45 0.3 0.2],...
  'fontsize',12,...
  'string','2');
%  Mf2 text box
ud.ph2.bmf2 = uicontrol('style','edit',...
  'Units', 'normalized',...
  'Parent',ud.ph2.panel,...
  'position',[0.65 0.45 0.3 0.2],...
  'fontsize',12,...
  'string','5');

% Mf text
ud.ph2.tb = uicontrol('style','text',...
  'Units', 'normalized',...
  'Parent',ud.ph2.panel,...
  'position',[0.05 0.1 0.3 0.2],...
  'fontsize',14,...
  'HorizontalAlignment','left',...
  'string','beta: ');
%  Mf1 text box
ud.ph2.bb1 = uicontrol('style','edit',...
  'Units', 'normalized',...
  'Parent',ud.ph2.panel,...
  'position',[0.3 0.1 0.3 0.2],...
  'fontsize',12,...
  'string','0');
%  Mf2 text box
ud.ph2.bb2 = uicontrol('style','edit',...
  'Units', 'normalized',...
  'Parent',ud.ph2.panel,...
  'position',[0.65 0.1 0.3 0.2],...
  'fontsize',12,...
  'string','2');


%% another panel
% panel
ud.ph3.panel = uipanel('Units', 'normalized',...
  'position',[0.69 0.1 0.29 0.3],...
  'fontsize',14,...
  'Title','Output');

% text
ud.ph3.tz = uicontrol('style','text',...
  'Units', 'normalized',...
  'Parent',ud.ph3.panel,...
  'position',[0.05 0.8 0.3 0.2],...
  'fontsize',14,...
  'HorizontalAlignment','left',...
  'string','Z: ');
% pop-up menu
ud.ph3.pz = uicontrol('style','popupmenu',...
  'Units', 'normalized',...
  'Parent',ud.ph3.panel,...
  'position',[0.3 0.8 0.65 0.2],...
  'fontsize',14,...
  'string',ud.out.names,...
  'Value',1);

% text
ud.ph3.tc = uicontrol('style','text',...
  'Units', 'normalized',...
  'Parent',ud.ph3.panel,...
  'position',[0.05 0.5 0.3 0.2],...
  'fontsize',14,...
  'HorizontalAlignment','left',...
  'string','C: ');
% pop-up menu
ud.ph3.pc = uicontrol('style','popupmenu',...
  'Units', 'normalized',...
  'Parent',ud.ph3.panel,...
  'position',[0.3 0.5 0.65 0.2],...
  'fontsize',14,...
  'string',ud.out.names,...
  'Value',2);
% Push button
ud.ph3.tc = uicontrol('style','pushbutton',...
  'Units', 'normalized',...
  'Parent',ud.ph3.panel,...
  'position',[0.15 0.2 0.65 0.2],...
  'fontsize',14,...
  'HorizontalAlignment','left',...
  'string','PLOT',...
  'Callback','irf_shock_rh(''plot'');');

set(gcf,'userdata',ud);
end


function ud = set_par(ud)

ud.N = 15;

% IN
ud.params.xpar = ud.ph1.px.String{ud.ph1.px.Value};
ud.params.ypar = ud.ph1.py.String{ud.ph1.py.Value};


xy = 'xy';

for l = 1:2
  switch ud.params.([xy(l),'par'])
    case 'th'
      ud.params.([xy(l),'val']) = linspace(str2double(ud.ph2.bth1.String),str2double(ud.ph2.bth2.String),ud.N);
    case 'Mf'
      ud.params.([xy(l),'val']) = linspace(str2double(ud.ph2.bmf1.String),str2double(ud.ph2.bmf2.String),ud.N);
    case 'beta'
      ud.params.([xy(l),'val']) = linspace(str2double(ud.ph2.bb1.String),str2double(ud.ph2.bb2.String),ud.N);
  end
end

% Out (redo)
ud.params.zpar = ud.ph3.pz.String(ud.ph3.pz.Value);
ud.params.cpar = ud.ph3.pc.String(ud.ph3.pc.Value);

end

function ud = plot_rankhug(ud)

% calculate first
upst = [];

% set fixed th if not selected
if ~strcmp(ud.params.xpar,'th') && ~strcmp(ud.params.ypar,'th')
  upst.th = str2double(ud.ph1.bth.String);
end
% set fixed Mf if not selected
if ~strcmp(ud.params.xpar,'Mf') && ~strcmp(ud.params.ypar,'Mf')
  upst.Mf = str2double(ud.ph1.bmf.String);
end
% set fixed beta if not selected
if ~strcmp(ud.params.xpar,'beta') && ~strcmp(ud.params.ypar,'beta')
  upst.beta = str2double(ud.ph1.bb.String);
end

z = zeros(ud.N,ud.N);
c = zeros(ud.N,ud.N);

zc = 'zc';

for jj = 1:ud.N % x
  for kk = 1:ud.N % y
    upst.(ud.params.xpar) = ud.params.xval(jj);
    upst.(ud.params.ypar) = ud.params.yval(kk);

    % do the calculation
    dnst = irf_shock_rh(upst);

    for l = 1:2
      switch ud.params.([zc(l),'par']){:}
        case ud.out.names{1} % nd/nu
          eval([zc(l),'(jj,kk)=dnst.nd;'])
        case ud.out.names{2} % Bd/Bu
          eval([zc(l),'(jj,kk)=norm(dnst.Bd);'])
        case ud.out.names{3} % Btd/Btu
          eval([zc(l),'(jj,kk)=dnst.Bd(2)/sind(upst.th);'])
        case ud.out.names{4} % thd
          eval([zc(l),'(jj,kk)=dnst.thd;'])
        case ud.out.names{5} % thd-th
          eval([zc(l),'(jj,kk)=dnst.thd-upst.th;'])
        case ud.out.names{6} % Mdf
          eval([zc(l),'(jj,kk)=dnst.Mdf;'])
        case ud.out.names{7} % MdA
          eval([zc(l),'(jj,kk)=dnst.Mda;'])
        case ud.out.names{8} % MdI
          eval([zc(l),'(jj,kk)=dnst.MdI;'])
        case ud.out.names{9} % betad
          eval([zc(l),'(jj,kk)=dnst.betad;'])
        case ud.out.names{10} % Vd
          eval([zc(l),'(jj,kk)=norm(dnst.Vd);'])
        case ud.out.names{11} % Vdn
          eval([zc(l),'(jj,kk)=dnst.Vd(1);'])
        case ud.out.names{12} % Vdt
          eval([zc(l),'(jj,kk)=dnst.Vd(2);'])
        case ud.out.names{13} % thdVn
          eval([zc(l),'(jj,kk)=dnst.thdVn;'])
        case ud.out.names{14} % Td/Tu
          eval([zc(l),'(jj,kk)=dnst.betad*norm(dnst.Bd)^2/(upst.beta*dnst.nd);'])
        case ud.out.names{15} % number of solutions
          eval([zc(l),'(jj,kk)=dnst.nSol;'])
      end
    end
  end
end

% edges for plotting (assume uniform)
dxval = median(diff(ud.params.xval));
dyval = median(diff(ud.params.yval));
ud.params.xvale = [ud.params.xval(1)-dxval/2,ud.params.xval+dxval/2];
ud.params.yvale = [ud.params.yval(1)-dyval/2,ud.params.yval+dyval/2];

%ud.hsf = surf(ud.ax,ud.params.xval,ud.params.yval,z');
zpl = nan(size(z)+[1,1]);
zpl(1:end-1,1:end-1) = z;
cpl = nan(size(z)+[1,1]);
cpl(1:end-1,1:end-1) = c;

ud.hsf = surf(ud.ax,ud.params.xvale,ud.params.yvale,zpl',cpl');
% ud.hsf.CData = c';

ud.hcb = colorbar(ud.ax,'Location','south');

ud.hcb.Position(2) = 0.1;
ud.hcb.AxisLocation = 'out';

% set x and y labels
xy = 'xy';
for l = 1:2
  switch ud.params.([xy(l),'par'])
    case 'th'
      guiLabel(ud.ax,xy(l),'$\theta_{\mathrm{Bn}}$ [$^\circ$]')
    case 'Mf'
      guiLabel(ud.ax,xy(l),'$M_f$')
    case 'beta'
      guiLabel(ud.ax,xy(l),'$\beta$')
  end
end

% z and c labels
guiLabel(ud.ax,'z',ud.out.labels{ud.ph3.pz.Value})
guiLabel(ud.hcb,'y',ud.out.labels{ud.ph3.pc.Value})


% set title
if ~strcmp(ud.params.xpar,'th') && ~strcmp(ud.params.ypar,'th')
  title(ud.ax,['$\theta_{\mathrm{Bn}} = ',num2str(upst.th),'^\circ$'],'interpreter','latex')
elseif ~strcmp(ud.params.xpar,'Mf') && ~strcmp(ud.params.ypar,'Mf')
  title(ud.ax,['$M_{f} = ',num2str(upst.Mf),'$'],'interpreter','latex')
elseif ~strcmp(ud.params.xpar,'beta') && ~strcmp(ud.params.ypar,'beta')
  title(ud.ax,['$\beta = ',num2str(upst.beta),'$'],'interpreter','latex')
end
end


function [] = guiLabel(ax,dir,lab)

labStr = [upper(dir),'Label'];

ax.(labStr).String = lab;
ax.(labStr).Interpreter = 'latex';
ax.FontSize = 16; % could be somewhere better

end

