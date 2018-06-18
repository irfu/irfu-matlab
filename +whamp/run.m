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
%     InputParameters.maxIterations - maximum number of iterations when
%                              searching for solution (default 50)
%
% Output is a structure:
% Output.InputParameters   - copy of InputParameters
% Output.PlasmaModel       - copy of PlasmaModel
% Output.(results)         - scalars, vectors or matrices [size(kpar) x size(kperp)]
%   (results) = kperp,kpar,f,E,Ex,Ey,Ez,B,Bx,By,Bz,EB,VGP,VGZ,
%               S,Sx,Sy,Sz,u,SGP,SGZ
%
% Examples:
%    Output = whamp.run([],[]) % default run (equivalent to below)
%
%              Oxygen = struct('m',16,'n',1,'t',10,'a',5,'vd',1);
%           Electrons = struct('m',0,'n',1,'t',100);
%         PlasmaModel = struct('B',100);
%	PlasmaModel.Species = {Oxygen,Electrons};
%     InputParameters = struct('fstart',0.1,'kperp',0,'kpar',0.0022,'useLog',0);
%              Output = whamp.run(PlasmaModel,InputParameters)
%


%% Check input number
if nargin==0
	help whamp.run;
	return;
elseif nargin ~= 2
	disp('whamp.run: incorrect number of input parameters, see help');
	return;
end

%% Define defaults
DefaultSpecies={...
	struct('m',16,'n',1,'t',10,'a',5,'b',0,'d',1,'vd',1),...
	struct('m',0,'n',1,'t',100,'a',1,'b',0,'d',1,'vd',0),...
	};
DefaultSpeciesParameters = struct('m',0,'n',1,'t',1,'a',1,'b',0.0,'d',1,'vd',0);
defaultB = 100;
DefaultInputParameters = struct('fstart',0.1,'kperp',[0.0 10 0.0],'kpar',[0.0022 10 0.0022],'varyKzFirst',1,'useLog',0,'maxIterations',50);
Output = [];

%% Define full PlasmaModel 
% use defaults where needed

% define plasma species
if isempty(PlasmaModel)
	PlasmaModel.B=defaultB;
	PlasmaModel.Species = DefaultSpecies;
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
			PlasmaModel.Species = DefaultSpecies;
		end
	else
		disp('whamp.run: input PlasmaModel does not define Species.');
		disp('           using default species (electrons, 1eV, 1cc).');
		PlasmaModel.Species = DefaultSpecies;		
	end
	speciesParameterFields = fieldnames(DefaultSpeciesParameters);
	for iSpecies = 1:numel(PlasmaModel.Species)
		for iFieldName = 1:numel(speciesParameterFields)
			fieldName = speciesParameterFields{iFieldName};
			if ~isfield(PlasmaModel.Species{iSpecies},fieldName)
				PlasmaModel.Species{iSpecies}.(fieldName) = ...
					DefaultSpeciesParameters.(fieldName);
			end
		end
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
	for iFieldName = 1:numel(inputParameterFields)
		fieldName = inputParameterFields{iFieldName};
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
	[' ' species_symbol(PlasmaModel.Species{iSpecies}.m) ' '] ...
	' n[cc]=' num2str(PlasmaModel.Species{iSpecies}.n) ...
	' T[eV]=' num2str(PlasmaModel.Species{iSpecies}.t) ...
	' Tperp/Tpar=' num2str(PlasmaModel.Species{iSpecies}.a) ...
	' d=' num2str(PlasmaModel.Species{iSpecies}.d) ...
	' b=' num2str(PlasmaModel.Species{iSpecies}.b) ...
	' vdrift/vtpar=' num2str(PlasmaModel.Species{iSpecies}.vd) ...
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
disp(['maxIterations = ' num2str(InputParameters.maxIterations) ]);
disp('-----------------------------------------');

%% Define WHAMP matrices for mexwhamp
% default values
  nWHAMP = zeros(1,10);
  tWHAMP = zeros(1,10);
  dWHAMP = zeros(1,10)+1;
  aaWHAMP = zeros(10,2)+1;
%  bWHAMP = zeros(1,10);
assWHAMP = zeros(1,10);
 vdWHAMP = zeros(1,10);

% plasma species matrices
for iSpecies = 1:numel(PlasmaModel.Species)
	  nWHAMP(iSpecies) = PlasmaModel.Species{iSpecies}.n*1e6;
	  tWHAMP(iSpecies) = PlasmaModel.Species{iSpecies}.t/1e3;
	  dWHAMP(iSpecies) = PlasmaModel.Species{iSpecies}.d;
   aaWHAMP(iSpecies,1) = PlasmaModel.Species{iSpecies}.a;
   aaWHAMP(iSpecies,2) = PlasmaModel.Species{iSpecies}.b;
	assWHAMP(iSpecies) = PlasmaModel.Species{iSpecies}.m;
	 vdWHAMP(iSpecies) = PlasmaModel.Species{iSpecies}.vd;
end
   fceWHAMP = 0.0279928*PlasmaModel.B; % fce in kHz
   pzlWHAMP = InputParameters.useLog;
zfirstWHAMP = InputParameters.varyKzFirst;
fstartWHAMP = InputParameters.fstart;

% define p
if numel(InputParameters.kperp) == 1
	pWHAMP = [InputParameters.kperp*[1 1] 10];
elseif numel(InputParameters.kperp) == 3
	pWHAMP = InputParameters.kperp([1 3 2]);
else
	disp('ERROR: InputParameters.kperp wrong format');
	return;
end

% define z
if numel(InputParameters.kpar) == 1
	zWHAMP = [InputParameters.kpar*[1 1] 10];
elseif numel(InputParameters.kpar) == 3
	zWHAMP = InputParameters.kpar([1 3 2]);
else
	disp('ERROR: InputParameters.kpar wrong format');
	return;
end

% define root finding parameters
maxIterationsWHAMP = int32(InputParameters.maxIterations);


%% call mexwhamp
[...
	kperpOUT,...   %1
 	kparOUT,...    %2
 	fOUT,...       %3
	ExOUT,...      %4
	EyOUT,...      %5
	EzOUT,...      %6
	BxOUT,...      %7
	ByOUT,...      %8
	BzOUT,...      %9
 	SxOUT,...      %10
 	SyOUT,...      %11
 	SzOUT,...      %12
 	EBOUT,...      %13
 	VGPOUT,...     %14
 	VGZOUT,...     %15
 	SGPOUT,...     %16
 	SGZOUT,...     %17
 	uOUT,...       %18
 	flagSolutionFoundOUT,...      %19
 	flagTooHeavilyDampedOUT,...   %20
	flagNoConvergenceOUT,...       %21
	]=whamp.mexwhamp(...
	fceWHAMP,...      %1
	pzlWHAMP,...      %2
	zfirstWHAMP,...   %3
	nWHAMP,...        %4
	tWHAMP,...        %5
	dWHAMP,...        %6
	aaWHAMP,...       %7
	assWHAMP,...      %8
	vdWHAMP,...       %9
	pWHAMP,...        %10
	zWHAMP,...        %11
	fstartWHAMP,...   %12
    maxIterationsWHAMP...%13
	);


%% Define Output

Output.InputParameters = InputParameters;
Output.PlasmaModel     = PlasmaModel;
if InputParameters.useLog == 1
	Output.kperp           = 10.^kperpOUT;
	Output.kpar            = 10.^kparOUT;
else
	Output.kperp           = kperpOUT;
	Output.kpar            = kparOUT;
end
Output.f               = fOUT.';
Output.Ex              = ExOUT.';
Output.Ey              = EyOUT.';
Output.Ez              = EzOUT.';
Output.Bx              = BxOUT.';
Output.By              = ByOUT.';
Output.Bz              = BzOUT.';
Output.Sx              = SxOUT.';
Output.Sy              = SyOUT.';
Output.Sz              = SzOUT.';
Output.EB              = EBOUT.';
Output.VGP             = VGPOUT.';
Output.VGZ             = VGZOUT.';
Output.SGP             = SGPOUT.';
Output.SGZ             = SGZOUT.';
Output.u               = uOUT.';
Output.flagSolutionFound    = flagSolutionFoundOUT.';
Output.flagTooHeavilyDamped = flagTooHeavilyDampedOUT.';
Output.flagNoConvergence    = flagNoConvergenceOUT.';

% set to NaN output
outFields = {'f','Ex','Ey','Ez','Bx','By','Bz',...
	'Sx','Sy','Sz','EB','VGP','VGZ','SGP','SGZ','u'};
indNaN = (Output.flagNoConvergence==1);
if any(indNaN)
	for i=1:numel(outFields)
		Output.(outFields{i})(indNaN) = NaN;
	end
end

function symbol = species_symbol(mass)
if (mass==0)
	symbol = 'e-';
elseif (mass==1)
	symbol = 'H+';
elseif (mass==2)
	symbol = 'He++';
elseif (mass==4)
	symbol = 'He+';
elseif (mass==16)
	symbol = 'O+';
else
	symbol=[];
end

