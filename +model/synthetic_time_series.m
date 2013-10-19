function timeSeries=synthetic_time_series(varargin)
% MODEL.SYNTHETIC_TIME_SERIES generate synthetic time series
%   
%	timeSeries = MODEL.SYNTHETIC_TIME_SERIES(InputParameters)
%
% 	InputParameters is a structure with fields 
%		samplingFrequency         [Hz] (default 10)
%		timeInterval              [s]  (default 3600)
%		peakFrequency             [Hz] (default 1)
%		peakHalfWidth             [Hz] (default 0.5)
%		components                [#]  (deafult 1)
%       amplitude                 complex vector of length # (default 1)
%
%	Components tell how many components of signal to model.
%	If amplitude is single scalar, all components are random with given
%	amplitude. If amplitude is vector of the length equal to the number of 
%	components then signal is generated with specified relation between 
%	the amplitudes of separate signals but random phase. 
%
%	timeSeries = MODEL.SYNTHETIC_TIME_SERIES('field1',field1Value,...)
%		InputParameters fields can be specified explicitely.
%		Shortenings can be used: fs=samplingFrequency, t=timeInterval, 
%		f=peakFrequency, df=peakHalfWidth, c=components and a=amplitude.
%
%	Example:
%		Inp = struct('samplingFrequency',10,'timeInterval',100,'peakFrequency',2,'peakHalfWidth',1);
%		ts = model.synthetic_time_series(Inp);
%		ts = model.synthetic_time_series('fs',100,'f',10);   % sampling freq 100Hz, peak freq 10 Hz
%       ts = model.synthetic_time_series('n',2);             % two random components
%		ts = model.synthetic_time_series('n',2,'a',[1 -1i]); % right hand polarized

%% Defaults
samplingFrequency = 10;
timeInterval = 3600;
peakFrequency = 1;
peakHalfWidth = 0.1;
components = 1;
amplitude = 1;
randomPhaseForEachComponent = false;

if nargin == 0 && nargout == 0,
	help model.synthetic_time_series;
	return;
elseif nargin == 1 && isstruct(varargin{1}),
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
			otherwise
				irf.log('critical','unrecognized input');
				return;
		end
		args(1:2)=[];
	end
end

ts = 0:1/samplingFrequency:timeInterval;ts=ts(:);
fs = 0:1/timeInterval:samplingFrequency;fs=fs(:);
timeSeries = zeros(numel(ts),components);
timeSeries(:,1) = ts;

if numel(amplitude) ~= components,
	if numel(amplitude) == 1, % use the same amplitude for all components
		amplitude = amplitude(ones(1,components));
	else
		irf.log('critical','error in amplitude input');
		return;
	end
else
	randomPhaseForEachComponent = false; 
end

if randomPhaseForEachComponent
	randPhase = exp(1i*rand(numel(fs),components)*2*pi);
else
	randPhase = exp(1i*rand(numel(fs),1)*2*pi)*ones(1,components);
end

for iComponent = 1:components
	amplitudeComplex = amplitude(iComponent)*sqrt(2);
	fftSignal = amplitudeComplex * ...
		exp(-abs(fs-peakFrequency)/peakHalfWidth);
	
	fftSignalComplex = fftSignal.*randPhase(:,iComponent);
	
	timeSeries(:,iComponent+1) = real(ifft(fftSignalComplex))';
end

