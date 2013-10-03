function Output = run(PlasmaModel,InputParameters)
% WHAMP.RUN run WHAMP code 
%
% Output = WHAMP.RUN(PlasmaModel,InputParameters);
%
% PlasmaModel is a structure:
%    PlasmaModel.B        - B field in nT                 /default 10nT/
%    PlasmaModel.Species  - cell array of Species         /default 1 species/
%
% Species is a structure:
%    Species.m      - mass/mass_proton (0 means electron) /default electron/
%    Species.n      - density in cm^-3                    /default 1/
%    Species.t      - temperature in eV                   /default 1/
%    Species.a      - anisotropy, Tperp/Tpar              /default 1/ 
%    Species.d      - loss cone parameter 1               /default ??/
%    Species.b      - loss cone parameter 2               /default ??/
%    Species.vd     - Vpar_drift/Vthermal                 /default 0/
%
% InputParameters is a structure: 
%     InputParameters.fstart - start frequency
%     InputParameters.kperp  - kp (perpendicular wave vector) as scalar or
%                              [kp_start kp_step kp_end]
%     InputParameters.kpar   - kz (parallel wave vector) or 
%                              [kz_start kz_step kz_end]
%     InputParameters.varyKzFirst - 1 for each kp value, step through all kz values
%                                 - 0 for each kz value, step through all kp values
%     InputParameters.useLog - 1 input log10(p) and log10(z),(default)
%                              0 input given as p and z
%
% Output is a structure:
% Output.InputParameters   - copy of InputParameters
% Output.PlasmaModel       - copy of PlasmaModel
% Output.(results)         - scalars, vectors or matrices 
%   (results) = kperp,kpar,f,E,Ex,Ey,Ez,B,Bx,By,Bz,EB,VGP,VGZ,
%               S,Sx,Sy,Sz,u,SGP,SGZ
%
% Examples:
%    Output = whamp.run([],[]); % default run
%
%        Electrons = struct('m',0,'n',1,'t',1);
%           Oxygen = struct('m',16,'n',1,'t',2);
%      PlasmaModel = struct('B',10,'Species',{Electrons,Oxygen});
%  InputParameters = struct('fstart',0.4,'kperp',-1,'kpar',-1);
%           Output = whamp.run(PlasmaModel,InputParameters);
%

%% Check input number
if nargin==0
	help whamp.run;
	return;
elseif nargin ~= 2,
	disp('whamp.run: incorrect number of input parameters, see help');
	return;
end

%% Define defaults
DefaultSpecies=struct('m',0,'n',1,'t',1,'a',1,'b',1,'d',1,'vd',0);
defaultB = 10;
DefaultInputParameters = struct('fstart',0.5,'kperp',-1,'kpar',-1,'varyKzFirst',1,'useLog',1);
Output = [];

%% Define full PlasmaModel 
% use defaults where needed

% define plasma species
if isempty(PlasmaModel)
	PlasmaModel.B=defaultB;
	PlasmaModel.Species = {DefaultSpecies};
else
	if isfield(PlasmaModel,'Species')
		if iscell(PlasmaModel.Species)
			% everything ok
		elseif isstruct(PlasmaModel.Species)
			PlasmaModel.Species = {PlasmaModel.Species}; % make cell array of lenght 1
		else
			disp('whamp.run: input PlasmaModel not properly defined!');
			disp('           using default model.');
			PlasmaModel.B=defaultB;
			PlasmaModel.Species = {DefaultSpecies};
		end
	else
		disp('whamp.run: input PlasmaModel does not define Species.');
		disp('           using default species (electrons, 1eV, 1cc).');
		PlasmaModel.Species = {DefaultSpecies};		
	end
end

% define magnetic field 
if ~isfield(PlasmaModel,'B')
	PlasmaModel.B=defaultB;
end

%% Define full InputParameters
% define plasma species
if isempty(InputParameters)
	InputParameters = DefaultInputParameters;
else
	inputParameterFields = fieldnames(DefaultInputParameters);
	for fieldName = inputParameterFields
		if ~isfield(InputParameters,fieldName)
			InputParameters.(fieldName) = DefaultInputParameters.(fieldName);
		end
	end
end

%% Summarize input
disp('----------------- WHAMP run -------------');
disp(' ');
disp('PlasmaModel');
disp('-----------');
disp(['  B=' num2str(PlasmaModel.B) ' nT']);
disp(['    Plasma consist of ' num2str(numel(PlasmaModel.Species)) ' components']);
for iSpecies = 1:numel(PlasmaModel.Species)
disp(['  ' num2str(iSpecies) '.' ...
	' m[mp]=' num2str(PlasmaModel.Species{iSpecies}.m) ...
	' n[cc]=' num2str(PlasmaModel.Species{iSpecies}.n) ...
	' T[eV]=' num2str(PlasmaModel.Species{iSpecies}.t) ...
	]);
end
disp(' ');
disp('InputParameters');
disp('---------------');
disp(['       fstart = ' num2str(InputParameters.fstart) ]);
disp(['        kperp = ' num2str(InputParameters.kperp) ]);
disp(['         kpar = ' num2str(InputParameters.kpar) ]);
disp(['  varyKzFirst = ' num2str(InputParameters.varyKzFirst) ]);
disp(['       useLog = ' num2str(InputParameters.useLog) ]);
disp('-----------------------------------------');
%% Define WHAMP matrices for mexwhamp
% default values
  nWHAMP = zeros(1,10);
  tWHAMP = zeros(1,10);
  dWHAMP = zeros(1,10)+1;
  aWHAMP = zeros(1,10)+1;
  bWHAMP = zeros(1,10);
assWHAMP = zeros(1,10);
 vdWHAMP = zeros(1,10);

% plasma species matrices
for iSpecies = 1:numel(PlasmaModel.Species)
	  nWHAMP(iSpecies) = PlasmaModel.Species{iSpecies}.n;
	  tWHAMP(iSpecies) = PlasmaModel.Species{iSpecies}.t;
	  dWHAMP(iSpecies) = PlasmaModel.Species{iSpecies}.d;
	  aWHAMP(iSpecies) = PlasmaModel.Species{iSpecies}.a;
	  bWHAMP(iSpecies) = PlasmaModel.Species{iSpecies}.b;
	assWHAMP(iSpecies) = PlasmaModel.Species{iSpecies}.m;
	 vdWHAMP(iSpecies) = PlasmaModel.Species{iSpecies}.vd;
end
fceWHAMP = 0.0279928*PlasmaModel.B; % fce in kHz
pzlWHAMP = InputParameters.varyKzFirst;

% define p
if numel(InputParameters.kperp) == 1,
	pWHAMP = InputParameters.kperp;
elseif numel(InputParameters.kperp) == 1,
	pWHAMP = InputParameters.kperp([1 3 2]);
else
	disp('ERROR: InputParameters.kperp wrong format');
	return;
end

% define z
if numel(InputParameters.kpar) == 1,
	zWHAMP = InputParameters.kpar;
elseif numel(InputParameters.kpar) == 1,
	zWHAMP = InputParameters.kpar([1 3 2]);
else
	disp('ERROR: InputParameters.kpar wrong format');
	return;
end

%% call mexwhamp


%% Define Output






