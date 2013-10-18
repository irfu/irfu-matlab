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
%	timeSeries = MODEL.SYNTHETIC_TIME_SERIES('field1',field1Value,...)
%		InputParameters fields can be specified explicitely.
%		Shortenings can be used: fs=samplingFrequency, t=timeInterval, 
%		f=peakFrequency, df=peakHalfWidth, c=components and a=amplitude.
%
%	Example:
%		Inp = struct('samplingFrequency',10,'timeInterval',100,'peakFrequency',2,'peakHalfWidth',1);
%		ts = model.synthetic_time_series(Inp);
%		ts = model.synthetic_time_series('fs',100,'f',10);
%

%% Defaults
samplingFrequency = 10;
timeInterval = 3600;
peakFrequency = 1;
peakHalfWidth = 0.1;
components = 1;
amplitude = 1;

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

ts = 0:1/samplingFrequency:timeInterval;
fs = 0:1/timeInterval:samplingFrequency;
timeSeries = zeros(numel(ts),components);
timeSeries(:,1) = ts';

if numel(amplitude) ~= components,
	if numel(amplitude) == 1, % use the same amplitude for all components
		amplitude = amplitude(ones(1,components));
	else
		irf.log('critical','error in amplitude input');
		return;
	end
end

for iComponent = 1:components
	amplitudeComplex = amplitude(iComponent)*sqrt(2);
	fftSignal = amplitudeComplex * ...
		exp(-abs(fs-peakFrequency)/peakHalfWidth);
	randPhase = exp(1i*rand(1,numel(fftSignal))*2*pi);
	
	fftSignalComplex = fftSignal.*randPhase;
	
	timeSeries(:,iComponent+1) = real(ifft(fftSignalComplex))';
end

