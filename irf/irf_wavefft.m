

% NAME
%   irf_wavefft: Short-Time Fourier Transform. 
% SYNOPSIS
%   [S, F, T] = irf_wavefft(x, [window], [frame_overlap], [frame_length],
%               [fs])
% DESCRIPTION
%   Short-Time Fourier Transform
% INPUTS
%   x      (vector) of size Nx1. --- onedimension array
%   window (string) the window function such as rectwin, hamming.  
%           if not specified, equivalent to hamming
%   frame_overlap 
%          (scalar) the length of each frame overlaps in micro second. 
%           The default is frame_length / 2. 
%   frame_length
%          (scalar) the length of each frame in micro second. 
%           The default is 20ms.  
%   fs     (scalar) the sampling rate in Hz. The default assumes 16kHz
% OUTPUTS
%   S      (matrix) of size MxK where M is the size of each frame(fft))
%           and K is the number of frames. S is spectrogram (complex number)
%   F      (vector) of size Mx1 that is a vector of frequencies (y-axis)
%           in Hz. 
%   T      (vector) of size 1xK whose value corresponds to the center of 
%           each frame (x-axis) in sec.

%----writen by Huishan Fu at IRFU (2012-08-20)----

 
function [S, F, T] = irf_wavefft(x, window, frame_overlap, frame_length, fs, show)
 %% Initialization
 N = length(x);
 if ~exist('window', 'var') || isempty(window)
     window = 'hamming';
 end
 if ~exist('fs', 'var') || isempty(fs)
     fs = 16000;
 end
 if ~exist('frame_length', 'var') || isempty(frame_length)
     frame_length = 20;
 end
 if ~exist('frame_overlap', 'var') || isempty(frame_overlap)
     frame_overlap = frame_length / 2;
 end
 if ~exist('show', 'var') || isempty(show)
     show = 0;
 end
 nsample = round(frame_length  * fs / 1000); % convert ms to points
 noverlap = round(frame_overlap * fs / 1000); % convert ms to points
 window   = eval(sprintf('%s(nsample)', window)); % e.g., hamming(nsample)

 %% spectrogram
 [S, F, T] = gen_spectrogram(x, window, noverlap, nsample, fs); % below



function [S, F, T] = gen_spectrogram(x, window, noverlap, nsample, fs)
 N = length(x);
 S = [];
 pos = 1;
 while (pos+nsample <= N)
     frame = x(pos:pos+nsample-1);
     pos = pos + (nsample - noverlap);
     Y = fft(frame .* window, nsample);
     % see also goertzel. signal/spectrogram.m is using it. 
     S = [S Y(1:round(nsample/2), 1)]; % half is enough, another half is just mirror
 end
 [M, K] = size(S);
 F = (0:round(nsample/2)-1)' / nsample * fs; % [0, fs/2) Hz 
 % F = psdfreqvec(nsample, fs, 'half');
 T = (round(nsample/2):(nsample-noverlap):N-1-round(nsample/2))/fs;

