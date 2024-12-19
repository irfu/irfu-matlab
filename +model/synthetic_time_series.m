function timeSeries=synthetic_time_series(varargin)
% MODEL.SYNTHETIC_TIME_SERIES generate synthetic time series
%
%	timeSeries = MODEL.SYNTHETIC_TIME_SERIES(InputParameters)
%
% 	InputParameters is a structure with fields
%		samplingFrequency         [Hz]        (default 100)
%		timeInterval              [s]         (default 100)
%		peakFrequency             [Hz]        (default 10)
%		peakHalfWidth             [Hz]        (default 1)
%		components                [#]         (default 1)
%       amplitude                 [s.u.^2/Hz] (default 1)
%		gyrotropy                 true/false  (default false)
%
%	Components tell how many components of signal to model.
%	Amplitude unit is singal unit squared per hertz.
%   Amplitude is a complex scalar or complex vector of length equal
%	to the number of componets. If amplitude is scalar, all components have the
%	same amplitude but are random with respect to each other. If amplitude
%	is vector then signal is generated with specified relation between
%	the amplitudes of separate components but a random phase.
%	If signal is set to gyrotropic, then amplitude is randomly rotated
%	in the plane spanned by 1st and 2nd component (x and y).
%
%	timeSeries = MODEL.SYNTHETIC_TIME_SERIES('field1',field1Value,...)
%		InputParameters fields can be specified explicitely.
%		Shortenings can be used: fs=samplingFrequency, t=timeInterval,
%		f=peakFrequency, df=peakHalfWidth, c=components, a=amplitude,
%		g=gyrotropy.
%
%	Example:
%		Inp = struct('samplingFrequency',666,'timeInterval',99.9,'peakFrequency',3.14,'peakHalfWidth',.1);
%		ts = model.synthetic_time_series(Inp);
%			irf_plot(ts); % to visualize
%		ts = model.synthetic_time_series('fs',25,'f',4);     % sampling freq 25Hz, peak freq 4 Hz
%       ts = model.synthetic_time_series('n',2);             % two random components
%		ts = model.synthetic_time_series('n',2,'a',[1 -1i]); % right hand polarized
%		ts = model.synthetic_time_series('n',2,'a',[1 -.3i],'g',true); % gyrotropic elliptically polarized right hand wave

%% Defaults
samplingFrequency = 100;
timeInterval = 100;
peakFrequency = 10;
peakHalfWidth = 1;
components = 1;
amplitude = 1;
randomPhaseForEachComponent = false; %#ok<NASGU>
gyrotropicSignal = false;

%% Input check
if nargin == 0 && nargout == 0
  help model.synthetic_time_series;
  return;
elseif nargin == 1 && isstruct(varargin{1})
  InputParameters = varargin{1};
  inputParameterFields = fieldnames(InputParameters);
  for j=1:numel(inputParameterFields)
    fieldname = inputParameterFields{j};
    eval([fieldname ' = InputParameters.' fieldname ';']);
  end
elseif nargin > 1
  args = varargin;
  while numel(args)>=2
    switch lower(args{1})
      case {'samplingfrequency','fs'}
        samplingFrequency = args{2};
      case {'timeinterval','t'}
        timeInterval = args{2};
      case {'peakfrequency','f'}
        peakFrequency = args{2};
      case {'peakhalfwidth','df'}
        peakHalfWidth = args{2};
      case {'components','n'}
        components = args{2};
      case {'amplitude','a'}
        amplitude = args{2};
      case {'gyrotropy','g'}
        gyrotropicSignal = args{2};
      otherwise
        irf.log('critical','unrecognized input');
        return;
    end
    args(1:2)=[];
  end
end

%% Initialization
ts = 0:1/samplingFrequency:timeInterval;ts=ts(:);
fs = 0:1/timeInterval:samplingFrequency;fs=fs(:);
nPoints = numel(ts);

%% Scale amplitude to get right units
% Sso that amplitude corresponds to spectral density
% in units - [signal units ^2 / Hz]
amplitude = amplitude/sqrt(norm(amplitude))*sqrt(1/timeInterval)*sqrt(2);

%% Initialize fftSignal
if numel(amplitude) ~= components
  if isscalar(amplitude) % use the same amplitude for all components
    fftSignal = amplitude(ones(nPoints,components));
    randomPhaseForEachComponent = true;
  else
    irf.log('critical','error in amplitude input');
    return;
  end
else
  fftSignal = ones(nPoints,1)*amplitude;
  randomPhaseForEachComponent = false;
end

%% Apply gyrotropization if needed
if gyrotropicSignal
  if components == 1
    irf.log('critical','should be at least 2 components for gyrotropic signal!');
    error('should be at least 2 components for gyrotropic signal!');
  end
  randEllipseAngle = exp(1i*rand(nPoints,1)*2*pi);
  realP = real(randEllipseAngle);
  imagP = imag(randEllipseAngle);
  newX = sum(fftSignal(:,1:2).*[realP -imagP],2);
  newY = sum(fftSignal(:,1:2).*[imagP realP],2);
  fftSignal(:,1:2)=[newX newY];
end

%% Apply exponential envelope
fftSignal = fftSignal .* ...
  (	sqrt(...
  (1/sqrt(pi)) ...
  * exp(-abs(fs-peakFrequency).^2/peakHalfWidth^2)...
  ) ...
  * ones(1,components)...
  );

%% Apply random phase
if randomPhaseForEachComponent
  randomPhase = exp(1i*rand(nPoints,components)*2*pi);
else
  randomPhase = exp(1i*rand(nPoints,1)*2*pi)*ones(1,components);
end
fftSignal = fftSignal.*randomPhase;

%% Construct time series
timeSeries = nPoints * real(ifft(fftSignal,[],1));
timeSeries = [ts timeSeries]; % add time axis as first column

