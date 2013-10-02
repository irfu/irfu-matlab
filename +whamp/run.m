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
%     InputParameters.kp     - kp (perpendicular wave vector) as scalar or
%                              [kp_start kp_step kp_end]
%     InputParameters.kz     - kz (parallel wave vector) or 
%                              [kz_start kz_step kz_end]
%     InputParameters.varyKzFirst - 1 for each kp value, step through all kz values
%                                 - 0 for each kz value, step through all kp values
%     InputParameters.useLog - 1 input log10(p) and log10(z),(default)
%                              0 input given as p and z
%
% Output.InputParameters
% Output.PlasmaModel
% Output.p
% Output.z
% Output.f
% Output.E
% Output.Ex
% Output.Ey
% Output.Ez
% Output.B
% Output.Bx
% Output.By
% Output.Bz
% Output.EB
% Output.VGP
% Output.VGZ
% Output.S
% Output.Sx
% Output.Sy
% Output.Sz
% Output.u
% Output.SGP
% Output.SGZ
%
% Examples:
%    Inp=struct('fstart',0.4,'kp',ad);Out=whamp.run(Plasma,Inp);

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
DeafultInputParameter = struct('fstart',0.5,'kp',-1,'kz',-1,'varyKzFirst',1,'useLog',1);

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
	inputParameterFields = fieldnames(DefaultInputParameter);
	for fieldName = inputParameterFields
		if ~isfield(InputParameters,fieldName)
			InputParameters.(fieldName) = DefaultInputParameters.(fieldName);
		end
	end
end

%% Summarize input
disp('PlasmaModel');
disp('-----------');
disp(['B=' PlasmaModel.B ' nT']);

%% Define WHAMP matrices for mexwhamp

%% call mexwhamp

%% Define Output






