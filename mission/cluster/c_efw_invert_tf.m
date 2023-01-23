function eout = c_efw_invert_tf(einp,filt,tm,method,edge)
%C_EFW_INVERT_TF  Invert EFW analogue filters response
%
% EOUT = C_EFW_INVERT_TF(EINP,FILT,[TM],[METHOD],[EDGE])
%
% Input:
%    EINP - data in the spinning frame, e.g. wE1p34
%    FILT - analogue filter used : 'L' (10Hz), 'M' (180Hz), 'H' (4kHz),
%           'BP' (50 Hz - 8 kHz bandpass), 'U' (unfiltered)
%      TM - optional, 'HX' or 'IB' (default).
%           For HX we uncorrect the filter group delay which was put
%           there by ISDAT. Note TM is mandatory if method is to be
%           chosen.
%  METHOD - optional, 'frequency' (default) or 'time'.
%           For time we use time domain filtering and block convolution
%           with "overlap save method".
%           If method other than default is to be used then parameter
%           TM must also be defined!
%    EDGE - optional, for 'time' method only. Edge treatment can either
%           be 'edge_zero'  (fill with zero outside of edges)
%           or 'edge_wrap'  (wrap the edge around outside of edges)
%           or 'edge_trunc' (truncate the edge outside of edges)
%           or 'edge_non'   (do nothing about the edges).
%           Default edge treatment is 'edge_wrap'.
%           * Truncate is usually considered a bad treatment and it is
%           therefor NOT recommended.
%           * No edge treatment is also considered a bad treatment and
%           it is also NOT recommended.

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

% Updated input check
narginchk(2, 5)

% Default tm to 'IB' and method to 'frequency' if not set by user.
if nargin == 2, tm = 'IB'; method = 'frequency'; end

if nargin == 3
  if strcmpi(tm,'HX')||strcmpi(tm,'IB')
    % tm correct but no method, Default method to 'old'.
    method='frequency';
  else
    % Optional tm is mandatory and must be 'HX' or 'IB' if selecting method.
    error('If you want to define method, you must also define the TM parameter correct first.');
  end
end

if nargin == 4
  if strcmpi(method,'time')
    % Default to wrap edges
    edge='edge_wrap';
  elseif strcmpi(method,'frequency')
    %(old method, do nothing with edge).
  else
    % Unknown method was specified.
    error('Unknown METHOD parameter was specified, valid is "frequency" (default) and "time".')
  end
end

if nargin == 5
  if strcmpi(method,'time')
    switch(upper(edge)) % Simple check
      case 'EDGE_ZERO'  % Ok
      case 'EDGE_WRAP'  % Ok
      case 'EDGE_TRUNC' % Usually bad and not recommended, but Ok
      case 'EDGE_NON'   % Usually bad and not recommended, but Ok
      otherwise
        % Wrong input
        error('If you define method as "time" and specify an edge treatment it must be a valid one.')
    end
  else
    % (frequency method, no edge parameter).
    error('EDGE treatment is only a valid parameter for "time" method.')
  end
end

% Check if data is continuous, require 3 us precision
dt=diff(einp(:,1)); maxJitter = max(abs(dt-median(dt)));
if maxJitter>1e-5
  error('data is not evenly sampled, max jitter %.2f us (>10 us)',maxJitter*1e6)
end

if strcmpi(method,'frequency')
  % Default method
  nfft = size(einp,1);
  if nfft/2==fix(nfft/2), nf = nfft/2;
  else, nf = (nfft+1)/2;
  end
  
  fsamp = c_efw_fsample(einp);
  f = fsamp*((1:nf) -1)'/nfft;
  
  tfinp = get_efw_tf(filt);
  
  tf =         interp1(tfinp(:,1),tfinp(:,4),f,'linear','extrap');    % Add real part.
  tf = tf + 1i*interp1(tfinp(:,1),tfinp(:,5),f,'linear','extrap');    % Add imaginary part.
  
  Pxx = fft(einp(:,2));
  if nfft/2==fix(nfft/2)
    Pxy = Pxx./[tf;flipud(conj(tf))];
  else
    Pxy = Pxx./[tf;flipud(conj(tf(1:end-1)))];
  end
  
  eout = einp;
  eout(:,2) = ifft(Pxy,'symmetric')/14.8;
elseif strcmpi(method,'time')
  % Translated to Matlab by T.Nilsson, IRFU, from LPP/Berkeley IDL code.
  
  fsamp = c_efw_fsample(einp);
  nk = 512;        % Number of samples in interpolated transfer function.
  df = fsamp/nk;
  frq = (0:nk-1)*df; frq(nk/2+2:end) = frq(nk/2+2:end)-nk*df;
  
  % Take absolute value: handle negative f
  myf = abs(frq);
  
  % Interpolation of the filter to the fixed size (nk) of the kernel
  % for positive frequencies myf
  tfinp = get_efw_tf(filt);
  c = interp1(tfinp(:,1), complex(tfinp(:,4),tfinp(:,5)), myf,'linear','extrap');
  
  take_conj = find(frq<0);
  n_take_conj = numel(take_conj);
  if(n_take_conj>0)
    c(take_conj) = conj(c(take_conj));
  else
    c = complex(ones(1,size(tfinp2,1)),zeros(1,size(tfinp2,1)));    % BUG: "tfinp2" not defined anywhere.
    warning('c filter had no negative freq. Not valid');
  end
  
  % FFT of real function must be real at f=0
  c(1) = abs( c(1) ) * sign( real( c(1) ) );
  
  s = complex(ones(1,nk),0);
  
  s = s./c;
  
  % ks = fft(s, 1) <- IDL "1" gives direction (+ => reverse, - => forward)
  ks = ifft(s)*length(s);
  
  % If we did everything right, the imaginary part truly is negligible
  pr = sum(real(ks).*real(ks));
  pi = sum(imag(ks).*imag(ks));
  
  if(pi/pr>1e-5)
    sprintf('*** cluster_efw_deconvo: Imag/Real for impulse response is = %d.',pi/pr);
  end
  
  kernel = real(ks);
  
  % zero time of the kernel is at index 0. Now, shift that to index nk/2
  % to get a kernel suitable for linear convolution and
  % to allow application of the window.
  kernel = circshift(kernel, [0 nk/2]);
  
  % As this is a continuous calibration, the window must be applied to the
  % kernel, rather than to the waveform.
  % note: application of window introduces a slight offset, which must be
  % removed from the signal afterwards.
  % Correcting for the offset in the kernel itself
  % would nullify the benefit of the window.
  kernel = kernel .* hann(nk,'periodic')';
  
  % normalize the kernel
  kernel = kernel/nk;
  
  % Preallocate output (including column 1: time)
  eout = einp;
  
  % Before applying the filter remove any DC levels as the Hann window
  % above could interfere with these and create an offset. The signal should
  % be uncoupled before we let the Hann applied kernel do its thing.
  einp_mean = mean(einp(:,2));
  einp(:,2) = einp(:,2) - einp_mean;
  
  % notes on edge behavior:
  % default: zero output when kernel overlaps edge
  % /edge_zero: usually good, but can emphasize low frequency trends, i.e.
  %                           artifacts of despin
  % /edge_wrap; similar to edge zero (based on analysis of cal signal).
  % /edge_truncate: usually bad
  
  % Normalize with the 14.8322 which corresponds to the 23 dB at
  % DC level.
  eout(:,2) = cluster_staff_fastconvol(einp(:,2)', kernel, edge, 0)' / 14.8322;
  
  % Remove any offset caused by Hann window and the kernel.
  eout(:,2) = eout(:,2) - mean(eout(:,2));
  
  % We should now have completely removed any DC levels, inferred from the
  % filter process. Then add the DC levels back again (the DC level were in
  % effect untouched by the filter process).
  eout(:,2) = eout(:,2) + einp_mean;
else
  % Should not happen as we have done checks at input.
  error('Method as defined in user input was incorrect. Method must be either "time" or "frequency", with "frequency" being default.')
end

% For HX data we uncorrect the filter delay put there by ISDAT
if strcmpi(tm,'HX')
  switch upper(filt)
    case 'M'
      eout(:,1) = eout(:,1) + 4.44e-3;
    case 'L'
      eout(:,1) = eout(:,1) + 0.08;
    otherwise
      error('HX data can use only L or M filter')
  end
end

end


function ao = cluster_staff_fastconvol(a, kernel, edge, blk_c)

% This function is a translated version of cluster_eff_fastconvol
% originally written in IDL by ???. Translated into Matlab code by
% T.Nilsson, IRFU, 20130702
% Input:
%        a - data seq.
%   kernel - normalized kernel of the filter
%     edge - edge_zero, zero outside of edges
%          - edge_wrap, wrap the edge outside of edges
%          - edge_truncate, simply truncate the edges
%          - edge_non, do nothing with edges.
%    blk_c - block length factor (if not set it defaults to 8 times the
%            filter length).

% Number of elements in the input signal to be filtered.
nbp = numel(a);

% Number of elements in the kernel (aka filter) to be used.
nk = numel(kernel);

% If blk_c factor is set use this this, if not use a factor 8.
if(blk_c~=0)
  b_length = blk_c * nk;
else
  b_length = 8 * nk;
end

% Pad the edges
switch(upper(edge))
  case 'EDGE_ZERO'
    % Not optimal treatment.
    ao = [zeros(1,nk/2), a, zeros(1,nk/2)];
  case 'EDGE_WRAP'
    % Default treatment
    ao = [wrev(a(1:nk/2)), a, wrev(a(nbp-nk/2+1:nbp))];
    % A possibility would be using fliplr but wrev is slightly quicker.
  case 'EDGE_TRUNC'
    % Not recommended treatment.
    ao = [zeros(1,nk/2) + a(1), a, zeros(1,nk/2)+ a(nbp)];
  case 'EDGE_NON'
    % Not recommended treatment.
    ao = a;
  otherwise
    % This should not happen as we have checked input already.
    error('Unknown EDGE treatment parameter.')
end

% perform the convolution.
if(b_length<nbp)
  % Use fast convolution
  ao = block_conv(kernel, ao, b_length);
  % If no edge padding, then zero out the edge
  if numel(ao)==nbp
    %warning('Block convolution but with no edge treatment.');
    ao(1:nk-1) = 0;
  end
else
  % Use brute force convolution, with correct attribute.
  ao = conv(ao, kernel,'full');
  
end

% shift back to zero delay
ao = circshift(ao, [0 -nk/2]);

% trim any edge padding
if(numel(ao)~=nbp)
  ao = ao(nk/2+1:nbp+nk/2);
end

end



function result = block_conv( filter, signal, blen )
% Linear Convolution computation via the Overlap-and-Save method. Ver 2.
% This function require: length(signal) > blen > length(filter).

% This separate function was originally created by Ilias Konsoulas who
% REQUIRED that the following copyright message remain intact. The code has
% since been modified by T.Nilsson, IRFU, 20130802.

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Copyright (c) 2012, Ilias Konsoulas
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% * Redistributions of source code must retain the above copyright
%   notice, this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright
%   notice, this list of conditions and the following disclaimer in
%   the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Make sure signal and filter are row vectors for subsequent processing:
signal = signal(:).';
filter = filter(:).';

SignLen = length(signal);
FiltLen = length(filter);

if(FiltLen>=blen)
  error('Block length must be longer than filter length.');
end

if(blen>=SignLen)
  error('Block length must be shorter than the signal to be filtered.');
end

% Number of block segments
m = ceil((SignLen+FiltLen-1)/blen);

% Pad signal with zeros to the least integer multiple of
% blen greater than or equal to FiltLen+SignLen-1.
signal = [signal  zeros(1, m*blen-SignLen)];

% This is the pre-allocated size of the output sequence
% result[n] where the convolution values will be stored.
% result may be a large vector therefor use x(1:n)=0 instead of
% x=zeros(1,n) which is not as quick for very large vectors.
result(1, m*blen) = 0;
y_current = zeros(1, blen+FiltLen);

for k=0:m-1
  if(k==0)
    % First x_block processing. We introduce some zeros in the
    % beginning of signal_block as necessary.
    signal_block = [zeros(1,FiltLen) signal(1:blen)];
    % Compute the circular convolution of size (blen + FiltLen).
    y_current = cconv(signal_block, filter, blen+FiltLen);
    result(1:blen) = y_current(FiltLen+1:end);
  else
    % Successive signal_blocks are of size blen+FiltLen overlapping
    % by FiltLen samples.
    signal_block = signal(k*blen - FiltLen+1:(k+1)*blen);
    % Compute the circular convolution of size blen+filtLen.
    y_current = cconv(signal_block, filter, blen+FiltLen);
    % Keep only the last blen samples of y_current.
    result(k*blen+1:(k+1)*blen) = y_current(FiltLen+1:end) ;
  end
end

% Finally, keep only the valid entries of result[n]:
result = result(1:SignLen+FiltLen-1);
end


function tfinp = get_efw_tf(filt)
% filt is one of 'L', 'M', 'H', 'BP', 'U'

switch upper(filt)
  case 'L' % 10 Hz filters
    tfinp = [...
      0.100017    -23.5153705     -2.8560805      0.0666334     -0.0033243
      0.101759    -23.5154116     -2.9063186      0.0666301     -0.0033827
      0.103531    -23.5154254     -2.9428639      0.0666278     -0.0034252
      0.105335    -23.5154142     -2.9945304      0.0666248     -0.0034853
      0.107170    -23.5153723     -3.0467564      0.0666219     -0.0035460
      0.109036    -23.5153477     -3.0997824      0.0666188     -0.0036077
      0.110936    -23.5153589     -3.1537505      0.0666153     -0.0036704
      0.112868    -23.5154155     -3.2084256      0.0666113     -0.0037340
      0.114834    -23.5154918     -3.2643405      0.0666070     -0.0037989
      0.116835    -23.5155662     -3.3209490      0.0666027     -0.0038647
      0.118870    -23.5156350     -3.3780805      0.0665983     -0.0039311
      0.120940    -23.5156345     -3.4364561      0.0665942     -0.0039990
      0.123047    -23.5155675     -3.4960729      0.0665906     -0.0040683
      0.125190    -23.5155186     -3.5574579      0.0665865     -0.0041396
      0.127371    -23.5155524     -3.6200170      0.0665817     -0.0042123
      0.129590    -23.5156436     -3.6834963      0.0665763     -0.0042860
      0.131847    -23.5157139     -3.7479628      0.0665709     -0.0043609
      0.134144    -23.5157367     -3.8129401      0.0665657     -0.0044364
      0.136480    -23.5157117     -3.8787656      0.0665608     -0.0045129
      0.138858    -23.5156393     -3.9458640      0.0665560     -0.0045909
      0.141277    -23.5155735     -4.0142839      0.0665510     -0.0046704
      0.143738    -23.5155532     -4.0845079      0.0665454     -0.0047519
      0.146241    -23.5155875     -4.1560888      0.0665391     -0.0048351
      0.148789    -23.5156552     -4.2284403      0.0665324     -0.0049190
      0.151380    -23.5157060     -4.3020795      0.0665257     -0.0050045
      0.154017    -23.5157061     -4.3767985      0.0665191     -0.0050913
      0.156700    -23.5156856     -4.4529498      0.0665124     -0.0051797
      0.159430    -23.5156990     -4.5307635      0.0665052     -0.0052700
      0.162207    -23.5157694     -4.6095222      0.0664974     -0.0053614
      0.165032    -23.5158496     -4.6899798      0.0664892     -0.0054547
      0.167907    -23.5158931     -4.7719179      0.0664810     -0.0055498
      0.170832    -23.5158862     -4.8550450      0.0664729     -0.0056462
      0.173808    -23.5158479     -4.9397060      0.0664648     -0.0057444
      0.176835    -23.5157906     -5.0255820      0.0664565     -0.0058441
      0.179916    -23.5157520     -5.1130621      0.0664478     -0.0059456
      0.183050    -23.5157904     -5.2022198      0.0664382     -0.0060489
      0.186238    -23.5158961     -5.2925593      0.0664278     -0.0061536
      0.189482    -23.5160222     -5.3848085      0.0664168     -0.0062605
      0.192783    -23.5161185     -5.4788098      0.0664057     -0.0063694
      0.196141    -23.5161475     -5.5742298      0.0663948     -0.0064799
      0.199558    -23.5161399     -5.6716355      0.0663837     -0.0065928
      0.203034    -23.5161433     -5.7700601      0.0663723     -0.0067068
      0.206570    -23.5161832     -5.8699634      0.0663602     -0.0068225
      0.210169    -23.5162862     -5.9719679      0.0663472     -0.0069405
      0.213830    -23.5164388     -6.0759269      0.0663333     -0.0070608
      0.217554    -23.5165290     -6.1822047      0.0663194     -0.0071838
      0.221344    -23.5165028     -6.2899189      0.0663060     -0.0073084
      0.225200    -23.5163976     -6.3987203      0.0662928     -0.0074344
      0.229122    -23.5163325     -6.5097839      0.0662787     -0.0075630
      0.233114    -23.5163631     -6.6231842      0.0662634     -0.0076941
      0.237174    -23.5164223     -6.7387563      0.0662473     -0.0078277
      0.241306    -23.5164649     -6.8566345      0.0662307     -0.0079639
      0.245509    -23.5164669     -6.9756178      0.0662140     -0.0081015
      0.249785    -23.5164919     -7.0968215      0.0661966     -0.0082415
      0.254137    -23.5165709     -7.2202540      0.0661780     -0.0083840
      0.258563    -23.5166648     -7.3455377      0.0661588     -0.0085286
      0.263067    -23.5167463     -7.4738723      0.0661389     -0.0086767
      0.267650    -23.5168146     -7.6042908      0.0661185     -0.0088271
      0.272312    -23.5169009     -7.7370595      0.0660972     -0.0089802
      0.277055    -23.5170194     -7.8720919      0.0660750     -0.0091359
      0.281881    -23.5171299     -8.0087558      0.0660521     -0.0092933
      0.286792    -23.5172288     -8.1480666      0.0660286     -0.0094538
      0.291787    -23.5173477     -8.2898473      0.0660041     -0.0096170
      0.296870    -23.5174847     -8.4338848      0.0659787     -0.0097828
      0.302041    -23.5175846     -8.5806412      0.0659526     -0.0099516
      0.307302    -23.5176344     -8.7293561      0.0659262     -0.0101227
      0.312655    -23.5176607     -8.8806918      0.0658991     -0.0102968
      0.318102    -23.5176646     -9.0357214      0.0658709     -0.0104750
      0.323643    -23.5176728     -9.1933024      0.0658418     -0.0106561
      0.329280    -23.5177178     -9.3546661      0.0658112     -0.0108415
      0.335016    -23.5177983     -9.5189381      0.0657792     -0.0110300
      0.340852    -23.5179037     -9.6852817      0.0657461     -0.0112208
      0.346789    -23.5180533     -9.8545783      0.0657116     -0.0114148
      0.352830    -23.5182470    -10.0257046      0.0656757     -0.0116108
      0.358976    -23.5184273    -10.1997018      0.0656388     -0.0118099
      0.365229    -23.5185441    -10.3769673      0.0656010     -0.0120128
      0.371591    -23.5186140    -10.5566347      0.0655625     -0.0122183
      0.378064    -23.5186560    -10.7398917      0.0655228     -0.0124279
      0.384649    -23.5186762    -10.9268616      0.0654817     -0.0126416
      0.391349    -23.5187024    -11.1172700      0.0654392     -0.0128591
      0.398166    -23.5187787    -11.3118712      0.0653945     -0.0130812
      0.405102    -23.5189073    -11.5089899      0.0653482     -0.0133059
      0.412159    -23.5190568    -11.7089388      0.0653002     -0.0135337
      0.419338    -23.5192307    -11.9127467      0.0652504     -0.0137656
      0.426643    -23.5194108    -12.1194129      0.0651989     -0.0140005
      0.434074    -23.5195744    -12.3306094      0.0651457     -0.0142405
      0.441636    -23.5197125    -12.5454164      0.0650908     -0.0144844
      0.449329    -23.5198356    -12.7635704      0.0650342     -0.0147319
      0.457155    -23.5199685    -12.9864457      0.0649754     -0.0149846
      0.465119    -23.5201243    -13.2126178      0.0649146     -0.0152407
      0.473221    -23.5202951    -13.4426757      0.0648516     -0.0155009
      0.481464    -23.5204568    -13.6768839      0.0647865     -0.0157656
      0.489850    -23.5205819    -13.9144967      0.0647196     -0.0160339
      0.498383    -23.5206726    -14.1571147      0.0646505     -0.0163076
      0.507065    -23.5207444    -14.4043019      0.0645790     -0.0165862
      0.515897    -23.5208415    -14.6547469      0.0645052     -0.0168682
      0.524884    -23.5209979    -14.9103668      0.0644281     -0.0171555
      0.534027    -23.5212032    -15.1696593      0.0643483     -0.0174465
      0.543329    -23.5214272    -15.4335750      0.0642656     -0.0177422
      0.552793    -23.5216746    -15.7031943      0.0641796     -0.0180439
      0.562423    -23.5219456    -15.9763992      0.0640908     -0.0183492
      0.572220    -23.5222141    -16.2552980      0.0639987     -0.0186604
      0.582187    -23.5224590    -16.5388194      0.0639038     -0.0189763
      0.592328    -23.5226600    -16.8266148      0.0638062     -0.0192966
      0.602646    -23.5228502    -17.1201446      0.0637051     -0.0196228
      0.613144    -23.5230797    -17.4179449      0.0636006     -0.0199531
      0.623824    -23.5233283    -17.7209610      0.0634924     -0.0202886
      0.634691    -23.5235443    -18.0302136      0.0633803     -0.0206305
      0.645746    -23.5237288    -18.3441597      0.0632650     -0.0209770
      0.656995    -23.5239533    -18.6645493      0.0631451     -0.0213299
      0.668439    -23.5242241    -18.9905040      0.0630208     -0.0216881
      0.680083    -23.5245327    -19.3204760      0.0628926     -0.0220499
      0.691929    -23.5248758    -19.6574454      0.0627593     -0.0224185
      0.703982    -23.5252553    -19.9994958      0.0626216     -0.0227918
      0.716245    -23.5256803    -20.3478769      0.0624788     -0.0231710
      0.728721    -23.5261182    -20.7032480      0.0623308     -0.0235569
      0.741415    -23.5265108    -21.0633006      0.0621787     -0.0239470
      0.754330    -23.5268919    -21.4307509      0.0620211     -0.0243442
      0.767469    -23.5273260    -21.8043311      0.0618580     -0.0247469
      0.780838    -23.5277639    -22.1841393      0.0616895     -0.0251551
      0.794440    -23.5281365    -22.5721736      0.0615151     -0.0255712
      0.808278    -23.5284908    -22.9657604      0.0613355     -0.0259921
      0.822358    -23.5289052    -23.3651311      0.0611499     -0.0264177
      0.836682    -23.5293526    -23.7721455      0.0609575     -0.0268501
      0.851257    -23.5298039    -24.1848416      0.0607594     -0.0272870
      0.866085    -23.5303012    -24.6060719      0.0605537     -0.0277314
      0.881171    -23.5308327    -25.0352715      0.0603406     -0.0281825
      0.896521    -23.5313914    -25.4712835      0.0601205     -0.0286390
      0.912137    -23.5319803    -25.9170147      0.0598918     -0.0291039
      0.928026    -23.5325803    -26.3686476      0.0596564     -0.0295730
      0.944191    -23.5331777    -26.8282185      0.0594132     -0.0300485
      0.960638    -23.5337579    -27.2968233      0.0591615     -0.0305314
      0.977372    -23.5343633    -27.7721776      0.0589021     -0.0310190
      0.994397    -23.5349850    -28.2569874      0.0586333     -0.0315140
      1.011718    -23.5355706    -28.7498988      0.0583561     -0.0320151
      1.029342    -23.5361543    -29.2502698      0.0580704     -0.0325213
      1.047272    -23.5368239    -29.7613095      0.0577735     -0.0330354
      1.065514    -23.5375678    -30.2807676      0.0574667     -0.0335550
      1.084075    -23.5383236    -30.8092825      0.0571498     -0.0340807
      1.102958    -23.5390800    -31.3489078      0.0568213     -0.0346144
      1.122171    -23.5398226    -31.8954346      0.0564838     -0.0351518
      1.141718    -23.5405430    -32.4522838      0.0561348     -0.0356961
      1.161606    -23.5412737    -33.0181170      0.0557749     -0.0362457
      1.181840    -23.5420796    -33.5917629      0.0554040     -0.0367989
      1.202427    -23.5429483    -34.1783677      0.0550189     -0.0373604
      1.223372    -23.5438081    -34.7742101      0.0546220     -0.0379268
      1.244682    -23.5446855    -35.3814777      0.0542115     -0.0384997
      1.266364    -23.5456277    -36.0007898      0.0537863     -0.0390792
      1.288423    -23.5466346    -36.6283336      0.0533489     -0.0396613
      1.310866    -23.5476651    -37.2681550      0.0528964     -0.0402498
      1.333700    -23.5486689    -37.9196095      0.0524293     -0.0408439
      1.356932    -23.5496654    -38.5814938      0.0519480     -0.0414421
      1.380569    -23.5507915    -39.2577208      0.0514486     -0.0420469
      1.404617    -23.5521227    -39.9436500      0.0509338     -0.0426532
      1.429084    -23.5535599    -40.6401866      0.0504032     -0.0432621
      1.453977    -23.5549830    -41.3507105      0.0498546     -0.0438766
      1.479305    -23.5562804    -42.0710491      0.0492917     -0.0444933
      1.505073    -23.5575111    -42.8070869      0.0487092     -0.0451164
      1.531290    -23.5587763    -43.5558080      0.0481085     -0.0457424
      1.557964    -23.5601174    -44.3165050      0.0474896     -0.0463699
      1.585102    -23.5615265    -45.0936808      0.0468487     -0.0470021
      1.612713    -23.5629802    -45.8815105      0.0461903     -0.0476339
      1.640805    -23.5644607    -46.6829579      0.0455117     -0.0482671
      1.669387    -23.5659248    -47.4999104      0.0448113     -0.0489028
      1.698466    -23.5674063    -48.3284934      0.0440919     -0.0495373
      1.728052    -23.5689975    -49.1732530      0.0433489     -0.0501728
      1.758153    -23.5706881    -50.0336012      0.0425823     -0.0508081
      1.788778    -23.5724854    -50.9070895      0.0417942     -0.0514407
      1.819937    -23.5744392    -51.8003529      0.0409779     -0.0520743
      1.851639    -23.5765117    -52.7072557      0.0401390     -0.0527038
      1.883893    -23.5786747    -53.6300434      0.0392752     -0.0533301
      1.916709    -23.5809279    -54.5720124      0.0383832     -0.0539546
      1.950096    -23.5832121    -55.5260083      0.0374697     -0.0545718
      1.984066    -23.5855335    -56.5001661      0.0365267     -0.0551862
      2.018626    -23.5878822    -57.4909221      0.0355574     -0.0557945
      2.053789    -23.5902633    -58.4975580      0.0345623     -0.0563951
      2.089564    -23.5927470    -59.5256024      0.0335353     -0.0569898
      2.125963    -23.5953809    -60.5695365      0.0324816     -0.0575739
      2.162995    -23.5981496    -61.6319208      0.0313985     -0.0581477
      2.200673    -23.6010060    -62.7154497      0.0302833     -0.0587117
      2.239007    -23.6038473    -63.8143882      0.0291422     -0.0592623
      2.278008    -23.6066498    -64.9356163      0.0279680     -0.0598019
      2.317689    -23.6094547    -66.0778609      0.0267616     -0.0603281
      2.358061    -23.6123103    -67.2369276      0.0255274     -0.0608371
      2.399137    -23.6153508    -68.4221795      0.0242551     -0.0613306
      2.440928    -23.6185877    -69.6243859      0.0229544     -0.0618030
      2.483446    -23.6220170    -70.8478557      0.0216210     -0.0622545
      2.526706    -23.6256882    -72.0964823      0.0202507     -0.0626843
      2.570719    -23.6294837    -73.3630351      0.0188520     -0.0630890
      2.615499    -23.6333642    -74.6575044      0.0174142     -0.0634705
      2.661058    -23.6372509    -75.9740034      0.0159442     -0.0638252
      2.707412    -23.6411983    -77.3114463      0.0144435     -0.0641508
      2.754573    -23.6453851    -78.6778652      0.0129035     -0.0644459
      2.802555    -23.6497810    -80.0665331      0.0113321     -0.0647070
      2.851373    -23.6542897    -81.4793212      0.0097283     -0.0649330
      2.901042    -23.6589958    -82.9217018      0.0080863     -0.0651220
      2.951575    -23.6637428    -84.3840080      0.0064183     -0.0652715
      3.002989    -23.6684624    -85.8775629      0.0047123     -0.0653810
      3.055299    -23.6732507    -87.4006459      0.0029712     -0.0654471
      3.108519    -23.6783191    -88.9482160      0.0012019     -0.0654653
      3.162667    -23.6836557    -90.5309805     -0.0006064     -0.0654333
      3.217758    -23.6890903    -92.1349355     -0.0024362     -0.0653497
      3.273809    -23.6948401    -93.7687183     -0.0042955     -0.0652105
      3.330836    -23.7010230    -95.4369976     -0.0061878     -0.0650116
      3.388856    -23.7073282    -97.1304647     -0.0081004     -0.0647533
      3.447887    -23.7136064    -98.8614964     -0.0100455     -0.0644324
      3.507946    -23.7199926   -100.6249630     -0.0120147     -0.0640457
      3.569052    -23.7263758   -102.4176749     -0.0140021     -0.0635917
      3.631222    -23.7328206   -104.2485806     -0.0160148     -0.0630651
      3.694474    -23.7396044   -106.1089769     -0.0180397     -0.0624631
      3.758829    -23.7468768   -108.0022840     -0.0200767     -0.0617813
      3.824305    -23.7545964   -109.9368168     -0.0221311     -0.0610141
      3.890921    -23.7624026   -111.8977223     -0.0241842     -0.0601670
      3.958697    -23.7706942   -113.9048342     -0.0262515     -0.0592265
      4.027654    -23.7793770   -115.9523133     -0.0283225     -0.0581926
      4.097813    -23.7880833   -118.0351145     -0.0303882     -0.0570675
      4.169193    -23.7967087   -120.1689205     -0.0324597     -0.0558410
      4.241817    -23.8052431   -122.3355380     -0.0345137     -0.0545204
      4.315706    -23.8139989   -124.5448724     -0.0365530     -0.0530958
      4.390882    -23.8234619   -126.8021964     -0.0385739     -0.0515586
      4.467367    -23.8342139   -129.0981024     -0.0405581     -0.0499101
      4.545185    -23.8460954   -131.4467159     -0.0425112     -0.0481403
      4.624358    -23.8577844   -133.8392579     -0.0444239     -0.0462613
      4.704910    -23.8687598   -136.2696452     -0.0462872     -0.0442799
      4.786866    -23.8802456   -138.7605632     -0.0481042     -0.0421706
      4.870249    -23.8925195   -141.3011450     -0.0498557     -0.0399403
      4.955085    -23.9047721   -143.8919526     -0.0515374     -0.0375928
      5.041398    -23.9172045   -146.5452922     -0.0531463     -0.0351164
      5.129215    -23.9298829   -149.2358141     -0.0546563     -0.0325354
      5.218562    -23.9431753   -151.9948980     -0.0560732     -0.0298211
      5.309464    -23.9570498   -154.8119606     -0.0573794     -0.0269861
      5.401951    -23.9724434   -157.6876665     -0.0585571     -0.0240307
      5.496048    -23.9881206   -160.6415350     -0.0596100     -0.0209434
      5.591784    -24.0020104   -163.6471871     -0.0605293     -0.0177606
      5.689188    -24.0162043   -166.7193803     -0.0612939     -0.0144674
      5.788290    -24.0336861   -169.8686177     -0.0618715     -0.0110560
      5.889117    -24.0540907   -173.0760500     -0.0622467     -0.0075591
      5.991700    -24.0746188   -176.3672120     -0.0624303     -0.0039637
      6.096070    -24.0948627   -179.7354777     -0.0624097     -0.0002881
      6.202259    -24.1140131   -183.1707641     -0.0621776      0.0034444
      6.310297    -24.1333815   -186.7063159     -0.0617091      0.0072560
      6.420217    -24.1560048   -190.3163884     -0.0609707      0.0110983
      6.532052    -24.1825699   -194.0135733     -0.0599446      0.0149609
      6.645834    -24.2100105   -197.8336206     -0.0586291      0.0188617
      6.761599    -24.2332521   -201.7336352     -0.0570576      0.0227448
      6.879380    -24.2596112   -205.7481530     -0.0551576      0.0266027
      6.999213    -24.2929254   -209.8440493     -0.0529133      0.0303578
      7.121133    -24.3315299   -214.0392499     -0.0503265      0.0339959
      7.245178    -24.3694903   -218.3910641     -0.0473942      0.0375522
      7.371382    -24.4061696   -222.8662770     -0.0441329      0.0409624
      7.499786    -24.4458870   -227.4676384     -0.0405188      0.0441684
      7.630425    -24.4905220   -232.2238016     -0.0365288      0.0471331
      7.763341    -24.5446908   -237.0550592     -0.0322278      0.0497310
      7.898572    -24.6110168   -242.0505321     -0.0275637      0.0519502
      8.036159    -24.6804407   -247.2651041     -0.0225471      0.0538086
      8.176142    -24.7463425   -252.7220472     -0.0171969      0.0552878
      8.318563    -24.8231786   -258.3566685     -0.0115825      0.0562097
      8.463465    -24.9190125   -264.0203858     -0.0059130      0.0564521
      8.610891    -25.0321425   -269.8893541     -0.0001082      0.0560263
      8.760886    -25.1640682   -276.1029756      0.0058667      0.0548691
      8.913493    -25.3049812   -282.5428601      0.0117910      0.0529981
      9.068758    -25.4729436   -289.2558599      0.0175625      0.0502748
      9.226728    -25.6643252   -296.1605161      0.0229674      0.0467572
      9.387450    -25.8844206   -303.1170212      0.0277492      0.0425396
      9.550972    -26.1501723   -310.2732582      0.0318431      0.0375836
      9.717341    -26.4741640   -317.8596140      0.0351888      0.0318406
      9.886609    -26.8457691   -325.8158185      0.0376132      0.0255468
      10.058825    -27.2513968   -333.8734786      0.0389602      0.0191088
      10.234041    -27.6741819   -341.7861027      0.0392615      0.0129191
      10.412310    -28.2043314   -349.7978800      0.0382703      0.0068874
      10.593683    -28.8235372   -357.9729982      0.0361869      0.0012807
      10.778216    -29.5056570   -366.2138140      0.0332781     -0.0036233
      10.965964    -30.2577785   -374.4230997      0.0297306     -0.0076463
      11.156981    -31.0460625   -382.5170582      0.0258975     -0.0107361
      11.351327    -31.8650050   -390.4142874      0.0220015     -0.0129156
      11.549057    -32.7524431   -397.8571051      0.0181867     -0.0141361
      11.750232    -33.7658038   -404.9080452      0.0145175     -0.0144709
      11.954911    -34.8694322   -411.8270561      0.0111569     -0.0141917
      12.163156    -35.9457687   -419.0039761      0.0082130     -0.0136708
      12.375027    -36.9645397   -426.1765291      0.0057289     -0.0129747
      12.590590    -38.0326747   -432.5332805      0.0037645     -0.0119637
      12.809907    -39.1788655   -438.0029060      0.0022847     -0.0107514
      13.033045    -40.3589605   -443.2640319      0.0011255     -0.0095289
      13.260069    -41.5809098   -448.2676498      0.0002520     -0.0083321
      13.491048    -42.8134873   -453.4982369     -0.0004413     -0.0072196
      13.726050    -44.0047555   -458.7436402     -0.0009586     -0.0062328
      13.965147    -45.1351105   -463.5414757     -0.0012964     -0.0053827
      14.208408    -46.2136913   -466.8423264     -0.0014168     -0.0046803
      14.455906    -47.3880088   -469.6665284     -0.0014376     -0.0040225
      14.707716    -48.7380653   -473.4816080     -0.0014571     -0.0033539
      14.963912    -50.0275465   -477.6967712     -0.0014651     -0.0027911
      15.224570    -51.0270903   -481.6023738     -0.0014723     -0.0023930
      15.489769    -51.8503091   -485.2045771     -0.0014733     -0.0020881
      15.759588    -52.7438861   -488.5041787     -0.0014355     -0.0018044
      16.034107    -53.6880397   -491.2700304     -0.0013642     -0.0015545
      16.313408    -54.8225063   -493.4716058     -0.0012487     -0.0013172
      16.597572    -56.1237724   -496.3646152     -0.0011308     -0.0010782
      16.886688    -57.2421330   -501.2012515     -0.0010706     -0.0008607
      17.180840    -58.0652424   -506.8120843     -0.0010457     -0.0006840
      17.480116    -58.8257217   -510.6389236     -0.0009977     -0.0005613
      17.784605    -59.8803168   -511.5083321     -0.0008911     -0.0004836
      18.094397    -61.0981611   -512.4261566     -0.0007811     -0.0004079
      18.409586    -62.2357130   -515.2789805     -0.0007022     -0.0003233
      18.730265    -63.1661695   -519.5070722     -0.0006506     -0.0002431
      19.056530    -63.9157087   -523.6108220     -0.0006112     -0.0001798
      19.388479    -64.5838405   -527.5479396     -0.0005761     -0.0001272
      19.726210    -65.3074864   -531.4060733     -0.0005367     -0.0000811
      20.069824    -66.2513221   -534.4456565     -0.0004846     -0.0000471
      20.419422    -67.5239257   -538.2573201     -0.0004203     -0.0000128
      20.775112    -69.0006374   -544.1785359     -0.0003538      0.0000259
      21.136997    -70.3380343   -550.9013597     -0.0002987      0.0000575
      21.505186    -71.2746327   -556.6980351     -0.0002616      0.0000785
      21.879787    -72.1256239   -560.4819311     -0.0002319      0.0000866
      22.260914    -73.1772752   -562.2275923     -0.0002030      0.0000830
      22.648682    -74.3102726   -563.1362060     -0.0001770      0.0000756
      23.043201    -75.4435548   -563.8744400     -0.0001545      0.0000684
      23.444595    -76.4928852   -565.5956383     -0.0001351      0.0000647
      23.852980    -77.5710341   -568.2657975     -0.0001165      0.0000626
      24.268478    -78.7319413   -571.2225604     -0.0000990      0.0000600
      24.691214    -80.0205684   -573.4382826     -0.0000833      0.0000550
      25.121315    -81.1639948   -574.6962165     -0.0000719      0.0000498
      25.558907    -81.9364717   -575.6226045     -0.0000650      0.0000466
      26.004120    -82.3631385   -576.0352228     -0.0000616      0.0000448
      26.457090    -82.4579375   -575.1538664     -0.0000616      0.0000434
      26.917950    -82.6154688   -573.3558982     -0.0000618      0.0000407
      27.386837    -83.0388037   -569.7122024     -0.0000612      0.0000349
      27.863894    -83.6517735   -564.9494728     -0.0000595      0.0000277
      28.349258    -84.3825565   -559.5864870     -0.0000569      0.0000202
      28.843079    -85.2259705   -555.3754168     -0.0000528      0.0000145
      29.345501    -86.1306947   -553.6914344     -0.0000480      0.0000117
      29.856674    -87.1315556   -553.3196362     -0.0000428      0.0000101
      30.376751    -88.3904995   -553.2021773     -0.0000371      0.0000087
      30.905890    -89.8893852   -551.5068853     -0.0000314      0.0000064
      31.444242    -91.6437963   -549.5134056     -0.0000258      0.0000043
      31.991976    -93.4763613   -548.5287201     -0.0000210      0.0000031
      32.549248    -95.3057436   -551.3298912     -0.0000168      0.0000034
      33.116230    -96.5843098   -555.0216085     -0.0000143      0.0000038
      33.693085    -96.4691918   -552.2374743     -0.0000147      0.0000032
      34.279991    -95.5124307   -545.2436259     -0.0000167      0.0000015
      34.877117    -95.0215069   -543.0583152     -0.0000177      0.0000009
      35.484646    -95.4306418   -549.0485871     -0.0000167      0.0000027
      36.102760    -95.6513561   -557.2979763     -0.0000158      0.0000049
      36.731640    -95.2134804   -562.7713843     -0.0000160      0.0000067
      37.371471    -94.8636993   -564.2254643     -0.0000165      0.0000074
      38.022453    -94.9112574   -566.1628971     -0.0000161      0.0000079
      38.684772    -95.1047449   -571.3470468     -0.0000150      0.0000091
      39.358627    -95.0462590   -577.7977435     -0.0000140      0.0000108
      40.044220    -94.8075586   -583.6873330     -0.0000131      0.0000126
      40.741756    -95.0095792   -587.9004325     -0.0000119      0.0000132
      41.451443    -95.8195307   -588.3047116     -0.0000108      0.0000121
      42.173492    -96.6512441   -585.0799844     -0.0000104      0.0000104
      42.908119    -97.0097786   -582.7727047     -0.0000104      0.0000096
      43.655540    -96.8255530   -585.4462216     -0.0000101      0.0000103
      44.415985    -96.3171001   -589.1941232     -0.0000100      0.0000116
      45.189674    -95.8076698   -589.9403176     -0.0000104      0.0000124
      45.976837    -95.1362095   -596.5228308     -0.0000097      0.0000146
      46.777718    -94.4643326   -610.3214220     -0.0000064      0.0000178
      47.592545    -94.1032210   -620.7069739     -0.0000032      0.0000195
      48.421566    -94.0416593   -626.0845097     -0.0000014      0.0000198
      49.265026    -93.9634828   -631.3694549      0.0000005      0.0000200
      50.123184    -86.9031567   -602.2547656     -0.0000210      0.0000400
      50.996284    -87.1694030   -590.5927882     -0.0000278      0.0000338
      51.884598    -87.6711356   -578.6925324     -0.0000323      0.0000258
      52.788383    -88.1637903   -568.2262285     -0.0000344      0.0000185
      53.707912    -88.3346497   -559.9190702     -0.0000360      0.0000131
      54.643459    -88.2635419   -560.3177821     -0.0000362      0.0000134
      55.595299    -88.1908936   -567.8440394     -0.0000344      0.0000182
      56.563725    -88.4264235   -581.9054330     -0.0000282      0.0000253
      57.549015    -89.2572730   -598.1668221     -0.0000182      0.0000293
      58.551472    -90.6663135   -611.1113902     -0.0000095      0.0000277
      59.571388    -92.3740698   -618.5955578     -0.0000048      0.0000236
      60.609070    -94.1645059   -618.2102879     -0.0000040      0.0000192
      61.664829    -95.6809434   -614.5919956     -0.0000044      0.0000159
      62.738976    -97.3121278   -612.5551098     -0.0000041      0.0000130
      63.831837    -99.2713339   -614.5544898     -0.0000029      0.0000105
      64.943733   -101.5225266   -623.7872432     -0.0000009      0.0000083
      66.074997   -103.1104292   -639.5363490      0.0000012      0.0000069
      67.225967   -103.4021222   -657.7431348      0.0000031      0.0000060
      68.396988   -102.7716411   -674.6904061      0.0000051      0.0000052
      69.588402   -101.5369798   -687.1856750      0.0000070      0.0000045
      70.800575   -100.7178604   -695.4737568      0.0000084      0.0000038
      72.033859   -100.1226484   -694.9310398      0.0000089      0.0000042
      73.288628    -99.4339028   -688.7213606      0.0000091      0.0000055
      74.565254    -98.4474993   -679.3661536      0.0000091      0.0000078
      75.864120    -97.4572235   -671.6166286      0.0000089      0.0000100
      77.185608    -96.7645906   -669.4673680      0.0000092      0.0000112
      78.530113    -97.0886599   -672.4321642      0.0000094      0.0000103
      79.898041    -99.4960666   -680.3632970      0.0000082      0.0000068
      81.289795   -103.3428003   -688.1476314      0.0000058      0.0000036
      82.705795   -106.7140332   -695.5843249      0.0000042      0.0000019
      84.146461   -108.5429529   -702.8462176      0.0000036      0.0000011
      85.612221   -107.6803761   -708.9275467      0.0000041      0.0000008
      87.103516   -105.0714992   -710.7896330      0.0000055      0.0000009
      88.620781   -102.5498111   -707.5083253      0.0000073      0.0000016
      90.164482   -100.7556556   -698.7482603      0.0000085      0.0000033
      91.735069   -100.0152166   -682.5514591      0.0000079      0.0000061
      93.333015    -99.7065559   -663.6247474      0.0000057      0.0000086
      94.958794   -100.7832412   -652.2377351      0.0000035      0.0000085
      96.612900   -102.9302101   -650.4282059      0.0000025      0.0000067
      98.295815   -105.6865228   -660.3489700      0.0000026      0.0000045];
  case 'M' % 180 Hz filters
    tfinp = [...
      0.100017    -23.4240888     -0.1257509      0.0674209     -0.0001480
      0.101759    -23.4240888     -0.1287127      0.0674209     -0.0001515
      0.103531    -23.4240888     -0.1316745      0.0674209     -0.0001549
      0.105335    -23.4240888     -0.1346363      0.0674209     -0.0001584
      0.107170    -23.4240888     -0.1375981      0.0674209     -0.0001619
      0.109036    -23.4240888     -0.1405599      0.0674209     -0.0001654
      0.110936    -23.4240888     -0.1435217      0.0674208     -0.0001689
      0.112868    -23.4240888     -0.1464835      0.0674208     -0.0001724
      0.114834    -23.4240888     -0.1494453      0.0674208     -0.0001759
      0.116835    -23.4240888     -0.1524071      0.0674208     -0.0001793
      0.118870    -23.4240888     -0.1553689      0.0674208     -0.0001828
      0.120940    -23.4240821     -0.1583307      0.0674209     -0.0001863
      0.123047    -23.4240699     -0.1612734      0.0674209     -0.0001898
      0.125190    -23.4240602     -0.1642458      0.0674210     -0.0001933
      0.127371    -23.4240576     -0.1672811      0.0674210     -0.0001968
      0.129590    -23.4240647     -0.1704001      0.0674209     -0.0002005
      0.131847    -23.4240762     -0.1736063      0.0674208     -0.0002043
      0.134144    -23.4240840     -0.1768407      0.0674208     -0.0002081
      0.136480    -23.4240858     -0.1801294      0.0674207     -0.0002120
      0.138858    -23.4240815     -0.1835189      0.0674208     -0.0002160
      0.141277    -23.4240755     -0.1869882      0.0674208     -0.0002200
      0.143738    -23.4240712     -0.1905262      0.0674208     -0.0002242
      0.146241    -23.4240702     -0.1940903      0.0674208     -0.0002284
      0.148789    -23.4240721     -0.1976411      0.0674208     -0.0002326
      0.151380    -23.4240735     -0.2012274      0.0674208     -0.0002368
      0.154017    -23.4240701     -0.2048626      0.0674208     -0.0002411
      0.156700    -23.4240627     -0.2085818      0.0674208     -0.0002454
      0.159430    -23.4240545     -0.2124156      0.0674209     -0.0002500
      0.162207    -23.4240484     -0.2163159      0.0674209     -0.0002545
      0.165032    -23.4240468     -0.2203269      0.0674209     -0.0002593
      0.167907    -23.4240490     -0.2244210      0.0674208     -0.0002641
      0.170832    -23.4240522     -0.2285759      0.0674208     -0.0002690
      0.173808    -23.4240583     -0.2328380      0.0674207     -0.0002740
      0.176835    -23.4240659     -0.2371702      0.0674207     -0.0002791
      0.179916    -23.4240718     -0.2415856      0.0674206     -0.0002843
      0.183050    -23.4240752     -0.2460926      0.0674205     -0.0002896
      0.186238    -23.4240741     -0.2506449      0.0674205     -0.0002949
      0.189482    -23.4240723     -0.2552397      0.0674205     -0.0003003
      0.192783    -23.4240736     -0.2598669      0.0674205     -0.0003058
      0.196141    -23.4240773     -0.2645158      0.0674204     -0.0003113
      0.199558    -23.4240813     -0.2692838      0.0674204     -0.0003169
      0.203034    -23.4240805     -0.2741749      0.0674203     -0.0003226
      0.206570    -23.4240729     -0.2791899      0.0674204     -0.0003285
      0.210169    -23.4240585     -0.2843270      0.0674205     -0.0003346
      0.213830    -23.4240402     -0.2894919      0.0674206     -0.0003407
      0.217554    -23.4240288     -0.2947441      0.0674206     -0.0003468
      0.221344    -23.4240267     -0.3000922      0.0674206     -0.0003531
      0.225200    -23.4240339     -0.3055344      0.0674205     -0.0003595
      0.229122    -23.4240451     -0.3111055      0.0674204     -0.0003661
      0.233114    -23.4240498     -0.3167436      0.0674203     -0.0003727
      0.237174    -23.4240495     -0.3224590      0.0674203     -0.0003794
      0.241306    -23.4240527     -0.3282689      0.0674202     -0.0003863
      0.245509    -23.4240623     -0.3341358      0.0674201     -0.0003932
      0.249785    -23.4240703     -0.3401669      0.0674200     -0.0004003
      0.254137    -23.4240725     -0.3463526      0.0674200     -0.0004076
      0.258563    -23.4240623     -0.3526644      0.0674200     -0.0004150
      0.263067    -23.4240515     -0.3590835      0.0674200     -0.0004225
      0.267650    -23.4240497     -0.3654914      0.0674200     -0.0004301
      0.272312    -23.4240552     -0.3719561      0.0674199     -0.0004377
      0.277055    -23.4240684     -0.3785573      0.0674197     -0.0004455
      0.281881    -23.4240781     -0.3853094      0.0674196     -0.0004534
      0.286792    -23.4240821     -0.3922724      0.0674195     -0.0004616
      0.291787    -23.4240830     -0.3993660      0.0674195     -0.0004699
      0.296870    -23.4240797     -0.4065313      0.0674194     -0.0004784
      0.302041    -23.4240773     -0.4138439      0.0674194     -0.0004870
      0.307302    -23.4240773     -0.4212593      0.0674193     -0.0004957
      0.312655    -23.4240781     -0.4288124      0.0674193     -0.0005046
      0.318102    -23.4240801     -0.4365339      0.0674192     -0.0005137
      0.323643    -23.4240791     -0.4443482      0.0674191     -0.0005229
      0.329280    -23.4240758     -0.4523267      0.0674191     -0.0005323
      0.335016    -23.4240726     -0.4604429      0.0674190     -0.0005418
      0.340852    -23.4240681     -0.4686667      0.0674190     -0.0005515
      0.346789    -23.4240637     -0.4770720      0.0674189     -0.0005614
      0.352830    -23.4240578     -0.4855885      0.0674189     -0.0005714
      0.358976    -23.4240513     -0.4942521      0.0674188     -0.0005816
      0.365229    -23.4240500     -0.5030987      0.0674188     -0.0005920
      0.371591    -23.4240518     -0.5120539      0.0674187     -0.0006025
      0.378064    -23.4240589     -0.5211912      0.0674185     -0.0006133
      0.384649    -23.4240692     -0.5304758      0.0674183     -0.0006242
      0.391349    -23.4240774     -0.5398805      0.0674182     -0.0006353
      0.398166    -23.4240825     -0.5495140      0.0674180     -0.0006466
      0.405102    -23.4240822     -0.5593068      0.0674179     -0.0006581
      0.412159    -23.4240782     -0.5692749      0.0674178     -0.0006699
      0.419338    -23.4240741     -0.5794254      0.0674177     -0.0006818
      0.426643    -23.4240715     -0.5896599      0.0674176     -0.0006939
      0.434074    -23.4240726     -0.6001112      0.0674175     -0.0007062
      0.441636    -23.4240762     -0.6107764      0.0674173     -0.0007187
      0.449329    -23.4240784     -0.6216426      0.0674172     -0.0007315
      0.457155    -23.4240830     -0.6327805      0.0674170     -0.0007446
      0.465119    -23.4240884     -0.6440546      0.0674168     -0.0007579
      0.473221    -23.4240918     -0.6554941      0.0674166     -0.0007713
      0.481464    -23.4240935     -0.6671351      0.0674165     -0.0007850
      0.489850    -23.4240906     -0.6789029      0.0674163     -0.0007989
      0.498383    -23.4240867     -0.6909221      0.0674162     -0.0008130
      0.507065    -23.4240853     -0.7031597      0.0674160     -0.0008274
      0.515897    -23.4240860     -0.7155846      0.0674158     -0.0008420
      0.524884    -23.4240882     -0.7283165      0.0674156     -0.0008570
      0.534027    -23.4240880     -0.7412331      0.0674154     -0.0008722
      0.543329    -23.4240854     -0.7543729      0.0674152     -0.0008877
      0.552793    -23.4240837     -0.7677747      0.0674150     -0.0009034
      0.562423    -23.4240824     -0.7813228      0.0674148     -0.0009194
      0.572220    -23.4240815     -0.7951454      0.0674146     -0.0009356
      0.582187    -23.4240808     -0.8092069      0.0674144     -0.0009522
      0.592328    -23.4240773     -0.8234727      0.0674142     -0.0009690
      0.602646    -23.4240755     -0.8380772      0.0674139     -0.0009861
      0.613144    -23.4240763     -0.8528886      0.0674137     -0.0010036
      0.623824    -23.4240769     -0.8679582      0.0674134     -0.0010213
      0.634691    -23.4240781     -0.8833425      0.0674131     -0.0010394
      0.645746    -23.4240762     -0.8989035      0.0674129     -0.0010577
      0.656995    -23.4240720     -0.9147938      0.0674126     -0.0010764
      0.668439    -23.4240685     -0.9309601      0.0674123     -0.0010954
      0.680083    -23.4240655     -0.9473539      0.0674120     -0.0011147
      0.691929    -23.4240670     -0.9641223      0.0674117     -0.0011345
      0.703982    -23.4240719     -0.9811181      0.0674113     -0.0011544
      0.716245    -23.4240782     -0.9983994      0.0674109     -0.0011748
      0.728721    -23.4240882     -1.0160200      0.0674105     -0.0011955
      0.741415    -23.4240968     -1.0338361      0.0674100     -0.0012165
      0.754330    -23.4241029     -1.0520472      0.0674096     -0.0012379
      0.767469    -23.4241080     -1.0706020      0.0674091     -0.0012597
      0.780838    -23.4241092     -1.0894508      0.0674087     -0.0012819
      0.794440    -23.4241077     -1.1087339      0.0674083     -0.0013046
      0.808278    -23.4241041     -1.1282466      0.0674079     -0.0013275
      0.822358    -23.4240999     -1.1480797      0.0674074     -0.0013509
      0.836682    -23.4241003     -1.1683228      0.0674070     -0.0013747
      0.851257    -23.4241040     -1.1888322      0.0674064     -0.0013988
      0.866085    -23.4241086     -1.2097755      0.0674059     -0.0014235
      0.881171    -23.4241118     -1.2310652      0.0674053     -0.0014485
      0.896521    -23.4241090     -1.2526224      0.0674048     -0.0014739
      0.912137    -23.4241029     -1.2746707      0.0674043     -0.0014998
      0.928026    -23.4240961     -1.2970480      0.0674037     -0.0015261
      0.944191    -23.4240892     -1.3198367      0.0674032     -0.0015529
      0.960638    -23.4240846     -1.3431194      0.0674026     -0.0015803
      0.977372    -23.4240803     -1.3666723      0.0674019     -0.0016080
      0.994397    -23.4240823     -1.3907336      0.0674012     -0.0016363
      1.011718    -23.4240955     -1.4206269      0.0674003     -0.0016715
      1.029342    -23.4241087     -1.4443881      0.0673995     -0.0016995
      1.047272    -23.4241163     -1.4698202      0.0673987     -0.0017294
      1.065514    -23.4241121     -1.4955608      0.0673979     -0.0017597
      1.084075    -23.4241020     -1.5217653      0.0673972     -0.0017905
      1.102958    -23.4240947     -1.5485367      0.0673964     -0.0018220
      1.122171    -23.4240941     -1.5756853      0.0673955     -0.0018539
      1.141718    -23.4240991     -1.6033630      0.0673946     -0.0018865
      1.161606    -23.4241077     -1.6315102      0.0673936     -0.0019196
      1.181840    -23.4241138     -1.6600340      0.0673926     -0.0019531
      1.202427    -23.4241177     -1.6891550      0.0673915     -0.0019874
      1.223372    -23.4241199     -1.7186794      0.0673905     -0.0020221
      1.244682    -23.4241230     -1.7487394      0.0673894     -0.0020574
      1.266364    -23.4241282     -1.7794495      0.0673882     -0.0020936
      1.288423    -23.4241316     -1.8105970      0.0673871     -0.0021302
      1.310866    -23.4241304     -1.8423967      0.0673859     -0.0021676
      1.333700    -23.4241278     -1.8747497      0.0673847     -0.0022057
      1.356932    -23.4241246     -1.9075485      0.0673834     -0.0022442
      1.380569    -23.4241240     -1.9410544      0.0673821     -0.0022836
      1.404617    -23.4241231     -1.9750225      0.0673807     -0.0023236
      1.429084    -23.4241192     -2.0095785      0.0673794     -0.0023642
      1.453977    -23.4241160     -2.0448342      0.0673779     -0.0024057
      1.479305    -23.4241190     -2.0805758      0.0673764     -0.0024477
      1.505073    -23.4241279     -2.1170505      0.0673747     -0.0024906
      1.531290    -23.4241367     -2.1541557      0.0673730     -0.0025342
      1.557964    -23.4241376     -2.1917681      0.0673714     -0.0025785
      1.585102    -23.4241345     -2.2301985      0.0673696     -0.0026236
      1.612713    -23.4241317     -2.2691826      0.0673679     -0.0026695
      1.640805    -23.4241329     -2.3088662      0.0673660     -0.0027161
      1.669387    -23.4241298     -2.3493371      0.0673641     -0.0027637
      1.698466    -23.4241046     -2.3902629      0.0673623     -0.0028118
      1.728052    -23.4240653     -2.4319837      0.0673605     -0.0028609
      1.758153    -23.4240351     -2.4744955      0.0673586     -0.0029109
      1.788778    -23.4240285     -2.5176905      0.0673564     -0.0029617
      1.819937    -23.4240379     -2.5618441      0.0673541     -0.0030136
      1.851639    -23.4240480     -2.6065921      0.0673516     -0.0030662
      1.883893    -23.4240539     -2.6521302      0.0673491     -0.0031197
      1.916709    -23.4240556     -2.6986239      0.0673465     -0.0031744
      1.950096    -23.4240582     -2.7457508      0.0673439     -0.0032298
      1.984066    -23.4240632     -2.7938419      0.0673411     -0.0032863
      2.018626    -23.4240688     -2.8427734      0.0673382     -0.0033438
      2.053789    -23.4240711     -2.8923862      0.0673353     -0.0034021
      2.089564    -23.4240714     -2.9430488      0.0673323     -0.0034616
      2.125963    -23.4240696     -2.9944012      0.0673292     -0.0035220
      2.162995    -23.4240683     -3.0466676      0.0673259     -0.0035834
      2.200673    -23.4240727     -3.1000374      0.0673225     -0.0036461
      2.239007    -23.4240830     -3.1541532      0.0673190     -0.0037097
      2.278008    -23.4240947     -3.2093915      0.0673153     -0.0037746
      2.317689    -23.4241029     -3.2655964      0.0673115     -0.0038406
      2.358061    -23.4241041     -3.3225704      0.0673076     -0.0039075
      2.399137    -23.4241040     -3.3807466      0.0673036     -0.0039759
      2.440928    -23.4241084     -3.4397113      0.0672994     -0.0040451
      2.483446    -23.4241162     -3.4997096      0.0672951     -0.0041156
      2.526706    -23.4241254     -3.5609498      0.0672906     -0.0041875
      2.570719    -23.4241333     -3.6230386      0.0672860     -0.0042604
      2.615499    -23.4241435     -3.6864262      0.0672811     -0.0043349
      2.661058    -23.4241585     -3.7509285      0.0672761     -0.0044106
      2.707412    -23.4241761     -3.8163457      0.0672709     -0.0044874
      2.754573    -23.4241881     -3.8832022      0.0672655     -0.0045659
      2.802555    -23.4241901     -3.9509931      0.0672600     -0.0046455
      2.851373    -23.4241892     -4.0199301      0.0672544     -0.0047264
      2.901042    -23.4241939     -4.0902331      0.0672485     -0.0048089
      2.951575    -23.4242057     -4.1615047      0.0672424     -0.0048926
      3.002989    -23.4242232     -4.2343090      0.0672360     -0.0049780
      3.055299    -23.4242431     -4.3084178      0.0672293     -0.0050649
      3.108519    -23.4242593     -4.3835609      0.0672225     -0.0051531
      3.162667    -23.4242700     -4.4603097      0.0672155     -0.0052431
      3.217758    -23.4242750     -4.5381022      0.0672082     -0.0053344
      3.273809    -23.4242785     -4.6172466      0.0672008     -0.0054272
      3.330836    -23.4242885     -4.6979985      0.0671930     -0.0055219
      3.388856    -23.4243086     -4.7798689      0.0671849     -0.0056179
      3.447887    -23.4243351     -4.8634451      0.0671764     -0.0057159
      3.507946    -23.4243593     -4.9484682      0.0671677     -0.0058155
      3.569052    -23.4243744     -5.0346724      0.0671587     -0.0059166
      3.631222    -23.4243868     -5.1227443      0.0671495     -0.0060198
      3.694474    -23.4244046     -5.2120438      0.0671399     -0.0061244
      3.758829    -23.4244256     -5.3029134      0.0671299     -0.0062309
      3.824305    -23.4244444     -5.3956465      0.0671196     -0.0063395
      3.890921    -23.4244561     -5.4896674      0.0671090     -0.0064496
      3.958697    -23.4244644     -5.5856751      0.0670980     -0.0065621
      4.027654    -23.4244770     -5.6833938      0.0670866     -0.0066765
      4.097813    -23.4244960     -5.7824712      0.0670748     -0.0067925
      4.169193    -23.4245145     -5.8836358      0.0670626     -0.0069109
      4.241817    -23.4245292     -5.9861492      0.0670500     -0.0070309
      4.315706    -23.4245484     -6.0904654      0.0670370     -0.0071529
      4.390882    -23.4245798     -6.1969686      0.0670233     -0.0072775
      4.467367    -23.4246199     -6.3049742      0.0670092     -0.0074038
      4.545185    -23.4246539     -6.4152156      0.0669945     -0.0075327
      4.624358    -23.4246775     -6.5273696      0.0669795     -0.0076638
      4.704910    -23.4246875     -6.6410217      0.0669641     -0.0077966
      4.786866    -23.4246955     -6.7570665      0.0669481     -0.0079322
      4.870249    -23.4247133     -6.8747210      0.0669315     -0.0080696
      4.955085    -23.4247401     -6.9944967      0.0669143     -0.0082095
      5.041398    -23.4247761     -7.1168167      0.0668963     -0.0083523
      5.129215    -23.4248122     -7.2408423      0.0668778     -0.0084971
      5.218562    -23.4248469     -7.3673954      0.0668586     -0.0086447
      5.309464    -23.4248826     -7.4961052      0.0668388     -0.0087949
      5.401951    -23.4249168     -7.6265317      0.0668183     -0.0089470
      5.496048    -23.4249502     -7.7597262      0.0667971     -0.0091022
      5.591784    -23.4249852     -7.8947983      0.0667751     -0.0092596
      5.689188    -23.4250264     -8.0323086      0.0667524     -0.0094198
      5.788290    -23.4250806     -8.1727082      0.0667287     -0.0095833
      5.889117    -23.4251406     -8.3150345      0.0667042     -0.0097490
      5.991700    -23.4251936     -8.4603186      0.0666789     -0.0099180
      6.096070    -23.4252373     -8.6081482      0.0666528     -0.0100900
      6.202259    -23.4252704     -8.7580541      0.0666259     -0.0102643
      6.310297    -23.4253051     -8.9111911      0.0665979     -0.0104423
      6.420217    -23.4253507     -9.0664175      0.0665690     -0.0106226
      6.532052    -23.4254033     -9.2243135      0.0665391     -0.0108060
      6.645834    -23.4254572     -9.3853822      0.0665081     -0.0109929
      6.761599    -23.4255044     -9.5486648      0.0664761     -0.0111823
      6.879380    -23.4255495     -9.7154158      0.0664429     -0.0113757
      6.999213    -23.4255982     -9.8851611      0.0664086     -0.0115724
      7.121133    -23.4256547    -10.0573291      0.0663731     -0.0117719
      7.245178    -23.4257208    -10.2332456      0.0663361     -0.0119755
      7.371382    -23.4257952    -10.4115964      0.0662979     -0.0121818
      7.499786    -23.4258727    -10.5930292      0.0662584     -0.0123916
      7.630425    -23.4259459    -10.7780729      0.0662175     -0.0126054
      7.763341    -23.4260104    -10.9655707      0.0661754     -0.0128219
      7.898572    -23.4260668    -11.1569485      0.0661318     -0.0130428
      8.036159    -23.4261220    -11.3517281      0.0660867     -0.0132675
      8.176142    -23.4261833    -11.5493043      0.0660400     -0.0134952
      8.318563    -23.4262535    -11.7511587      0.0659916     -0.0137277
      8.463465    -23.4263336    -11.9557739      0.0659415     -0.0139631
      8.610891    -23.4264250    -12.1639520      0.0658896     -0.0142025
      8.760886    -23.4265275    -12.3763812      0.0658358     -0.0144465
      8.913493    -23.4266333    -12.5917469      0.0657802     -0.0146937
      9.068758    -23.4267306    -12.8116470      0.0657226     -0.0149459
      9.226728    -23.4268185    -13.0354484      0.0656630     -0.0152023
      9.387450    -23.4269001    -13.2623842      0.0656017     -0.0154621
      9.550972    -23.4269880    -13.4941569      0.0655379     -0.0157272
      9.717341    -23.4270902    -13.7290824      0.0654721     -0.0159956
      9.886609    -23.4272033    -13.9681522      0.0654040     -0.0162684
      10.058825    -23.4273211    -14.2121636      0.0653332     -0.0165466
      10.234041    -23.4274337    -14.4595392      0.0652603     -0.0168283
      10.412310    -23.4275447    -14.7120030      0.0651847     -0.0171155
      10.593683    -23.4276612    -14.9688746      0.0651064     -0.0174073
      10.778216    -23.4277851    -15.2292988      0.0650257     -0.0177028
      10.965964    -23.4279120    -15.4953242      0.0649419     -0.0180043
      11.156981    -23.4280403    -15.7650924      0.0648554     -0.0183096
      11.351327    -23.4281786    -16.0396995      0.0647659     -0.0186199
      11.549057    -23.4283317    -16.3200152      0.0646729     -0.0189362
      11.750232    -23.4284920    -16.6041114      0.0645770     -0.0192563
      11.954911    -23.4286514    -16.8940058      0.0644776     -0.0195824
      12.163156    -23.4288061    -17.1889728      0.0643747     -0.0199137
      12.375027    -23.4289605    -17.4880793      0.0642688     -0.0202492
      12.590590    -23.4291313    -17.7936477      0.0641586     -0.0205912
      12.809907    -23.4293223    -18.1034706      0.0640449     -0.0209374
      13.033045    -23.4295163    -18.4187787      0.0639273     -0.0212891
      13.260069    -23.4297007    -18.7405793      0.0638054     -0.0216473
      13.491048    -23.4298759    -19.0667534      0.0636798     -0.0220097
      13.726050    -23.4300583    -19.3996143      0.0635495     -0.0223788
      13.965147    -23.4302542    -19.7382465      0.0634147     -0.0227535
      14.208408    -23.4304678    -20.0816093      0.0632757     -0.0231326
      14.455906    -23.4306902    -20.4325296      0.0631312     -0.0235191
      14.707716    -23.4309070    -20.7885202      0.0629823     -0.0239103
      14.963912    -23.4311185    -21.1508289      0.0628283     -0.0243075
      15.224570    -23.4313318    -21.5204977      0.0626686     -0.0247117
      15.489769    -23.4315599    -21.8950765      0.0625041     -0.0251202
      15.759588    -23.4318161    -22.2773627      0.0623332     -0.0255360
      16.034107    -23.4320910    -22.6664076      0.0621564     -0.0259578
      16.313408    -23.4323766    -23.0609856      0.0619742     -0.0263844
      16.597572    -23.4326624    -23.4641156      0.0617850     -0.0268189
      16.886688    -23.4329378    -23.8728560      0.0615901     -0.0272581
      17.180840    -23.4332108    -24.2887720      0.0613887     -0.0277036
      17.480116    -23.4334899    -24.7131584      0.0611799     -0.0281566
      17.784605    -23.4337853    -25.1432787      0.0609647     -0.0286141
      18.094397    -23.4340977    -25.5823322      0.0607414     -0.0290794
      18.409586    -23.4344190    -26.0291639      0.0605106     -0.0295511
      18.730265    -23.4347477    -26.4824227      0.0602726     -0.0300277
      19.056530    -23.4350984    -26.9457114      0.0600254     -0.0305129
      19.388479    -23.4354730    -27.4156026      0.0597706     -0.0310028
      19.726210    -23.4358553    -27.8937721      0.0595072     -0.0314991
      20.069824    -23.4362263    -28.3816257      0.0592343     -0.0320033
      20.419422    -23.4365763    -28.8759549      0.0589536     -0.0325119
      20.775112    -23.4369230    -29.3803481      0.0586628     -0.0330283
      21.136997    -23.4372879    -29.8935617      0.0583621     -0.0335510
      21.505186    -23.4376867    -30.4140443      0.0580523     -0.0340782
      21.879787    -23.4381186    -30.9459963      0.0577305     -0.0346140
      22.260914    -23.4385690    -31.4856289      0.0573990     -0.0351543
      22.648682    -23.4390345    -32.0349963      0.0570562     -0.0357012
      23.043201    -23.4395073    -32.5958434      0.0567009     -0.0362560
      23.444595    -23.4399698    -33.1643268      0.0563354     -0.0368148
      23.852980    -23.4404342    -33.7444404      0.0559568     -0.0373813
      24.268478    -23.4409182    -34.3346280      0.0555657     -0.0379536
      24.691214    -23.4414317    -34.9329748      0.0551630     -0.0385295
      25.121315    -23.4419700    -35.5442828      0.0547454     -0.0391134
      25.558907    -23.4425206    -36.1642913      0.0543155     -0.0397010
      26.004120    -23.4430811    -36.7955382      0.0538714     -0.0402944
      26.457090    -23.4436470    -37.4401491      0.0534112     -0.0408953
      26.917950    -23.4442106    -38.0936721      0.0529378     -0.0414991
      27.386837    -23.4448100    -38.7607169      0.0524475     -0.0421097
      27.863894    -23.4454536    -39.4394161      0.0519411     -0.0427248
      28.349258    -23.4461481    -40.1275526      0.0514202     -0.0433421
      28.843079    -23.4468604    -40.8306918      0.0508802     -0.0439662
      29.345501    -23.4475462    -41.5438765      0.0503251     -0.0445926
      29.856674    -23.4482227    -42.2699677      0.0497521     -0.0452232
      30.376751    -23.4489104    -43.0114282      0.0491588     -0.0458596
      30.905890    -23.4496319    -43.7632160      0.0485488     -0.0464968
      31.444242    -23.4503840    -44.5308628      0.0479174     -0.0471390
      31.991976    -23.4511550    -45.3121660      0.0472659     -0.0477838
      32.549248    -23.4519289    -46.1045192      0.0465965     -0.0484285
      33.116230    -23.4527223    -46.9142338      0.0459032     -0.0490777
      33.693085    -23.4535413    -47.7354887      0.0451908     -0.0497259
      34.279991    -23.4544019    -48.5714161      0.0444562     -0.0503749
      34.877117    -23.4553326    -49.4245879      0.0436965     -0.0510258
      35.484646    -23.4563043    -50.2893201      0.0429166     -0.0516737
      36.102760    -23.4573057    -51.1724677      0.0421102     -0.0523230
      36.731640    -23.4583127    -52.0717205      0.0412790     -0.0529713
      37.371471    -23.4593077    -52.9843984      0.0404254     -0.0536160
      38.022453    -23.4603220    -53.9174936      0.0395423     -0.0542608
      38.684772    -23.4613700    -54.8637743      0.0386361     -0.0548999
      39.358627    -23.4624379    -55.8268433      0.0377033     -0.0555347
      40.044220    -23.4635088    -56.8098369      0.0367405     -0.0561664
      40.741756    -23.4645665    -57.8064415      0.0357537     -0.0567900
      41.451443    -23.4656570    -58.8240680      0.0347351     -0.0574088
      42.173492    -23.4667996    -59.8598715      0.0336872     -0.0580197
      42.908119    -23.4680174    -60.9108047      0.0326128     -0.0586196
      43.655540    -23.4693226    -61.9860744      0.0315023     -0.0592124
      44.415985    -23.4706745    -63.0781406      0.0303633     -0.0597927
      45.189674    -23.4720672    -64.1904812      0.0291921     -0.0603612
      45.976837    -23.4735201    -65.3260003      0.0279855     -0.0609177
      46.777718    -23.4750104    -66.4762995      0.0267524     -0.0614567
      47.592545    -23.4765680    -67.6508227      0.0254824     -0.0619810
      48.421566    -23.4781752    -68.8470917      0.0241784     -0.0624880
      49.265026    -23.4798186    -70.0618873      0.0228439     -0.0629746
      50.123184    -23.4815152    -71.3042460      0.0214689     -0.0634427
      50.996284    -23.4832552    -72.5641489      0.0200648     -0.0638866
      51.884598    -23.4849480    -73.8471032      0.0186257     -0.0643073
      52.788383    -23.4864457    -75.1584126      0.0171462     -0.0647055
      53.707912    -23.4877566    -76.4894588      0.0156362     -0.0650765
      54.643459    -23.4891120    -77.8507496      0.0140835     -0.0654194
      55.595299    -23.4906699    -79.2374764      0.0124940     -0.0657293
      56.563725    -23.4925965    -80.6445762      0.0108738     -0.0660017
      57.549015    -23.4948167    -82.0833086      0.0092108     -0.0662369
      58.551472    -23.4970999    -83.5430737      0.0075185     -0.0664326
      59.571388    -23.4993779    -85.0300342      0.0057905     -0.0665879
      60.609070    -23.5016296    -86.5499286      0.0040213     -0.0667007
      61.664829    -23.5038767    -88.0922796      0.0022239     -0.0667675
      62.738976    -23.5061540    -89.6696052      0.0003851     -0.0667859
      63.831837    -23.5084829    -91.2767542     -0.0014877     -0.0667526
      64.943733    -23.5108863    -92.9085514     -0.0033871     -0.0666647
      66.074997    -23.5133890    -94.5789701     -0.0053274     -0.0665185
      67.225967    -23.5159547    -96.2756860     -0.0072924     -0.0663120
      68.396988    -23.5185861    -98.0053292     -0.0092878     -0.0660416
      69.588402    -23.5212954    -99.7742508     -0.0113185     -0.0657030
      70.800575    -23.5240514   -101.5690891     -0.0133665     -0.0652955
      72.033859    -23.5268723   -103.4054303     -0.0154470     -0.0648126
      73.288628    -23.5297497   -105.2776826     -0.0175505     -0.0642520
      74.565254    -23.5326945   -107.1803628     -0.0196674     -0.0636123
      75.864120    -23.5358143   -109.1296260     -0.0218119     -0.0628839
      77.185608    -23.5391216   -111.1102735     -0.0239632     -0.0620689
      78.530113    -23.5425988   -113.1304088     -0.0261258     -0.0611611
      79.898041    -23.5462379   -115.1984126     -0.0283040     -0.0601533
      81.289795    -23.5499047   -117.2982781     -0.0304762     -0.0590509
      82.705795    -23.5536416   -119.4493046     -0.0326570     -0.0578405
      84.146461    -23.5574648   -121.6440947     -0.0348329     -0.0565225
      85.612221    -23.5613976   -123.8762528     -0.0369912     -0.0550980
      87.103516    -23.5654982   -126.1659783     -0.0391444     -0.0535508
      88.620781    -23.5697225   -128.4953485     -0.0412685     -0.0518903
      90.164482    -23.5740344   -130.8738839     -0.0433650     -0.0501080
      91.735069    -23.5784061   -133.3121609     -0.0454346     -0.0481935
      93.333015    -23.5827995   -135.7895882     -0.0474513     -0.0461612
      94.958794    -23.5874065   -138.3317810     -0.0494259     -0.0439877
      96.612900    -23.5922577   -140.9295136     -0.0513401     -0.0416791
      98.295815    -23.5974497   -143.5755529     -0.0531777     -0.0392410
      100.008041    -23.6031234   -146.2938630     -0.0549430     -0.0366509
      101.750099    -23.6091579   -149.0616210     -0.0566093     -0.0339315
      103.522499    -23.6154700   -151.8925718     -0.0581738     -0.0310716
      105.325775    -23.6220506   -154.8037091     -0.0596316     -0.0280558
      107.160454    -23.6287190   -157.7685058     -0.0609561     -0.0249148
      109.027100    -23.6357034   -160.8209966     -0.0621463     -0.0216161
      110.926262    -23.6430306   -163.9463156     -0.0631791     -0.0181804
      112.858498    -23.6508006   -167.1354294     -0.0640353     -0.0146244
      114.824402    -23.6592160   -170.4224968     -0.0647058     -0.0109180
      116.824547    -23.6681639   -173.7803832     -0.0651670     -0.0071020
      118.859528    -23.6778369   -177.2262336     -0.0654032     -0.0031687
      120.929962    -23.6886191   -180.7841555     -0.0653926      0.0008950
      123.036453    -23.7001622   -184.4147568     -0.0651181      0.0050274
      125.179649    -23.7134580   -188.1710207     -0.0645499      0.0092685
      127.360168    -23.7282208   -192.0322375     -0.0636710      0.0135711
      129.578674    -23.7447557   -195.9908964     -0.0624631      0.0179003
      131.835831    -23.7640365   -200.0967439     -0.0608859      0.0222771
      134.132294    -23.7857611   -204.3121179     -0.0589361      0.0266257
      136.468765    -23.8108573   -208.6620138     -0.0565831      0.0309296
      138.845932    -23.8412784   -213.1879514     -0.0537773      0.0351747
      141.264511    -23.8753590   -217.8265570     -0.0505578      0.0392542
      143.725220    -23.9180879   -222.6670495     -0.0468342      0.0431675
      146.228790    -23.9680576   -227.6739867     -0.0426419      0.0468201
      148.775970    -24.0270015   -232.8395695     -0.0379946      0.0501279
      151.367523    -24.1001541   -238.2383903     -0.0328320      0.0530318
      154.004211    -24.1863051   -243.8119551     -0.0272544      0.0554175
      156.686844    -24.2896876   -249.5894517     -0.0212825      0.0571947
      159.416183    -24.4196258   -255.6292778     -0.0149215      0.0582388
      162.193085    -24.5671109   -261.8292675     -0.0084006      0.0585077
      165.018341    -24.7534738   -268.2899198     -0.0017265      0.0578273
      167.892822    -24.9715945   -274.9560559      0.0048741      0.0562074
      170.817383    -25.2257603   -281.7875904      0.0111930      0.0536359
      173.792877    -25.5323054   -288.8328271      0.0170738      0.0500598
      176.820190    -25.8828033   -296.0024327      0.0222710      0.0456574
      179.900253    -26.2832974   -303.2820729      0.0266206      0.0405537
      183.033966    -26.7530803   -310.6479493      0.0299364      0.0348683
      186.222260    -27.2635287   -318.0465373      0.0322266      0.0289696
      189.466095    -27.8456260   -325.4028580      0.0333584      0.0230100
      192.766434    -28.4813919   -332.7212092      0.0334756      0.0172623
      196.124252    -29.1648562   -339.9307070      0.0327003      0.0119467
      199.540573    -29.9080428   -346.9910225      0.0311391      0.0071942
      203.016403    -30.6915976   -353.8882202      0.0290365      0.0031091
      206.552765    -31.5137461   -360.6163907      0.0265636     -0.0002858
      210.150742    -32.3812003   -367.0766557      0.0238572     -0.0029617
      213.811386    -33.2716601   -373.3963407      0.0211075     -0.0050271
      217.535812    -34.1937074   -379.4313843      0.0184012     -0.0064914
      221.325089    -35.1388519   -385.2742438      0.0158255     -0.0074720
      225.180389    -36.0989846   -390.9174226      0.0134429     -0.0080509
      229.102844    -37.0789584   -396.3290509      0.0112768     -0.0082925
      233.093613    -38.0691697   -401.5535123      0.0093463     -0.0082845
      237.153900    -39.0703407   -406.6134557      0.0076452     -0.0080883
      241.284927    -40.0810653   -411.4407219      0.0061753     -0.0077470
      245.487900    -41.0978684   -416.1610122      0.0049074     -0.0073198
      249.764099    -42.1198705   -420.6772999      0.0038367     -0.0068306
      254.114777    -43.1478214   -425.0602710      0.0029348     -0.0063110
      258.541229    -44.1774650   -429.3108988      0.0021841     -0.0057833
      263.044800    -45.2121701   -433.4143681      0.0015665     -0.0052594
      267.626831    -46.2474011   -437.4004690      0.0010626     -0.0047538
      272.288666    -47.2859666   -441.2858590      0.0006548     -0.0042723
      277.031677    -48.3251229   -445.0257454      0.0003325     -0.0038204
      281.857330    -49.3668391   -448.7005155      0.0000771     -0.0034005
      286.767059    -50.4086717   -452.2458378     -0.0001182     -0.0030146
      291.762299    -51.4533968   -455.7067768     -0.0002660     -0.0026618
      296.844543    -52.4975144   -459.0790945     -0.0003743     -0.0023423
      302.015320    -53.5439089   -462.3522361     -0.0004498     -0.0020542
      307.276184    -54.5888230   -465.5444584     -0.0004997     -0.0017963
      312.628662    -55.6352662   -468.6722679     -0.0005292     -0.0015659
      318.074402    -56.6803050   -471.7057115     -0.0005420     -0.0013616
      323.614990    -57.7273274   -474.6983506     -0.0005428     -0.0011802
      329.252075    -58.7729687   -477.5895121     -0.0005334     -0.0010208
      334.987366    -59.8197552   -480.4095577     -0.0005168     -0.0008805
      340.822571    -60.8647289   -483.1558491     -0.0004951     -0.0007579
      346.759399    -61.9119568   -485.8395060     -0.0004698     -0.0006505
      352.799652    -62.9582253   -488.4830411     -0.0004427     -0.0005568
      358.945129    -64.0060477   -491.0880116     -0.0004144     -0.0004752
      365.197632    -65.0526468   -493.6238698     -0.0003856     -0.0004046
      371.559082    -66.1019877   -496.1206136     -0.0003570     -0.0003433
      378.031311    -67.1492556   -498.5263710     -0.0003290     -0.0002908
      384.616302    -68.1966086   -500.8742180     -0.0003019     -0.0002456
      391.315979    -69.2392849   -503.1583579     -0.0002762     -0.0002070
      398.132385    -70.2812082   -505.3866798     -0.0002520     -0.0001739
      405.067505    -71.3209172   -507.5782899     -0.0002293     -0.0001456
      412.123444    -72.3613376   -509.7144396     -0.0002081     -0.0001215
      419.302277    -73.3990065   -511.7483170     -0.0001883     -0.0001012
      426.606171    -74.4373062   -513.7129397     -0.0001701     -0.0000840
      434.037292    -75.4643927   -515.6197345     -0.0001535     -0.0000696
      441.597839    -76.4854912   -517.5252790     -0.0001385     -0.0000573
      449.290100    -77.4905770   -519.4382764     -0.0001250     -0.0000469
      457.116364    -78.4727535   -521.2806400     -0.0001129     -0.0000383
      465.078918    -79.4331429   -522.9944407     -0.0001021     -0.0000312
      473.180206    -80.3847958   -524.5670095     -0.0000922     -0.0000255
      481.422607    -81.3301535   -525.9406478     -0.0000832     -0.0000208
      489.808563    -82.2814286   -527.2362492     -0.0000750     -0.0000170
      498.340637    -83.2128211   -528.4723114     -0.0000677     -0.0000138
      507.021301    -84.1286635   -529.6954612     -0.0000612     -0.0000111
      515.853149    -85.0084352   -530.9619828     -0.0000555     -0.0000088
      524.838928    -85.8400831   -532.2911850     -0.0000506     -0.0000068
      533.981140    -86.6402818   -533.6540943     -0.0000463     -0.0000051
      543.282654    -87.4320038   -534.9921492     -0.0000423     -0.0000037
      552.746155    -88.2106462   -536.2425849     -0.0000388     -0.0000025
      562.374512    -88.9764095   -537.4352799     -0.0000355     -0.0000016
      572.170593    -89.8218909   -538.2137440     -0.0000323     -0.0000010
      582.137329    -90.7795092   -538.5558619     -0.0000289     -0.0000007
      592.277649    -91.8343567   -538.6870642     -0.0000256     -0.0000006
      602.594666    -92.8298325   -539.3486057     -0.0000228     -0.0000003
      613.091309    -93.6562772   -540.8617832     -0.0000208      0.0000003
      623.770874    -94.3603925   -542.3197682     -0.0000191      0.0000008
      634.636414    -94.9262134   -542.8344480     -0.0000179      0.0000009
      645.691223    -95.4470015   -542.3612247     -0.0000169      0.0000007
      656.938660    -95.8638988   -542.0240639     -0.0000161      0.0000006
      668.381958    -96.2176521   -542.3613715     -0.0000154      0.0000006
      680.024597    -96.5877925   -544.0590599     -0.0000148      0.0000010
      691.870056    -97.0783819   -546.4735555     -0.0000139      0.0000016
      703.921814    -97.6555936   -548.5023476     -0.0000130      0.0000019
      716.183533    -98.1579139   -549.9662231     -0.0000122      0.0000021
      728.658875    -98.5525590   -550.3246345     -0.0000116      0.0000021
      741.351501    -98.9123229   -549.7965109     -0.0000112      0.0000019
      754.265198    -99.3164414   -549.8098744     -0.0000107      0.0000018
      767.403870    -99.6940960   -551.5619805     -0.0000101      0.0000021
      780.771362    -99.9951629   -555.1638159     -0.0000097      0.0000026
      794.371765   -100.1741258   -558.2798781     -0.0000093      0.0000031
      808.209045   -100.2491860   -558.8409462     -0.0000092      0.0000031
      822.287354   -100.2876705   -558.3515578     -0.0000092      0.0000030
      836.610901   -100.3201344   -559.0594032     -0.0000091      0.0000031
      851.183899   -100.4664951   -562.1233009     -0.0000088      0.0000036
      866.010803   -100.7668879   -565.0923246     -0.0000083      0.0000039
      881.096008   -101.1314696   -566.2878203     -0.0000079      0.0000039
      896.443970   -101.4716757   -564.5261302     -0.0000077      0.0000035
      912.059204   -101.6747537   -563.2092575     -0.0000076      0.0000032
      927.946533   -101.6929806   -565.5956549     -0.0000074      0.0000036
      944.110596   -101.5102352   -571.1124228     -0.0000072      0.0000043
      960.556152   -101.1947013   -578.4556334     -0.0000068      0.0000054
      977.288269   -100.9839288   -584.4675497     -0.0000064      0.0000063
      994.311768   -100.8045728   -589.8095030     -0.0000059      0.0000070];
  case 'H' % 4 kHz filters
    tfinp = [...
      10.000000    -23.4050313     -0.1285904      0.0675690     -0.0001516
      10.150914    -23.4020964     -0.1255484      0.0675918     -0.0001481
      10.304106    -23.4005355     -0.1253550      0.0676040     -0.0001479
      10.459609    -23.4006131     -0.1285503      0.0676034     -0.0001517
      10.617459    -23.3991632     -0.1230364      0.0676147     -0.0001452
      10.777692    -23.3983803     -0.1240610      0.0676207     -0.0001464
      10.940342    -23.3988641     -0.1366511      0.0676169     -0.0001613
      11.105447    -23.3997189     -0.1451379      0.0676103     -0.0001713
      11.273044    -23.3999602     -0.1518430      0.0676084     -0.0001792
      11.443170    -23.3986899     -0.1561402      0.0676182     -0.0001843
      11.615864    -23.3984877     -0.1511815      0.0676198     -0.0001784
      11.791163    -23.3983615     -0.1465519      0.0676208     -0.0001730
      11.969108    -23.3974604     -0.1519438      0.0676278     -0.0001793
      12.149739    -23.3974365     -0.1553202      0.0676280     -0.0001833
      12.333096    -23.3975093     -0.1590475      0.0676274     -0.0001877
      12.519219    -23.3960979     -0.1690591      0.0676384     -0.0001996
      12.708152    -23.3967049     -0.1683149      0.0676337     -0.0001987
      12.899936    -23.3981238     -0.1621031      0.0676226     -0.0001913
      13.094614    -23.3955326     -0.1659070      0.0676428     -0.0001959
      13.292230    -23.3949544     -0.1719216      0.0676473     -0.0002030
      13.492829    -23.3960493     -0.1797384      0.0676387     -0.0002122
      13.696454    -23.3968964     -0.1836207      0.0676321     -0.0002167
      13.903153    -23.3967163     -0.1877317      0.0676335     -0.0002216
      14.112971    -23.3953886     -0.1920918      0.0676438     -0.0002268
      14.325956    -23.3947437     -0.1980583      0.0676488     -0.0002338
      14.542154    -23.3936292     -0.2032719      0.0676575     -0.0002400
      14.761616    -23.3917569     -0.2071833      0.0676720     -0.0002447
      14.984389    -23.3929777     -0.2081329      0.0676625     -0.0002458
      15.210525    -23.3940598     -0.2109977      0.0676541     -0.0002491
      15.440073    -23.3938397     -0.2187934      0.0676558     -0.0002584
      15.673086    -23.3925569     -0.2225892      0.0676657     -0.0002629
      15.909614    -23.3915137     -0.2252110      0.0676739     -0.0002660
      16.149713    -23.3920383     -0.2274340      0.0676698     -0.0002686
      16.393435    -23.3918452     -0.2369021      0.0676712     -0.0002798
      16.640835    -23.3915146     -0.2470260      0.0676737     -0.0002918
      16.891968    -23.3918765     -0.2465866      0.0676709     -0.0002912
      17.146892    -23.3915899     -0.2455283      0.0676732     -0.0002900
      17.405662    -23.3911181     -0.2451489      0.0676768     -0.0002896
      17.668338    -23.3921801     -0.2535346      0.0676685     -0.0002994
      17.934978    -23.3923271     -0.2592720      0.0676674     -0.0003062
      18.205642    -23.3917661     -0.2631120      0.0676717     -0.0003108
      18.480391    -23.3912841     -0.2713498      0.0676754     -0.0003205
      18.759286    -23.3911550     -0.2785143      0.0676764     -0.0003290
      19.042390    -23.3913955     -0.2845507      0.0676745     -0.0003361
      19.329766    -23.3899797     -0.2894442      0.0676855     -0.0003419
      19.621480    -23.3894387     -0.2937407      0.0676897     -0.0003470
      19.917595    -23.3902359     -0.2972540      0.0676834     -0.0003511
      20.218180    -23.3903920     -0.3100639      0.0676821     -0.0003663
      20.523301    -23.3904919     -0.3205193      0.0676813     -0.0003786
      20.833026    -23.3906505     -0.3239168      0.0676800     -0.0003826
      21.147426    -23.3905213     -0.3306409      0.0676810     -0.0003906
      21.466570    -23.3904869     -0.3352539      0.0676812     -0.0003960
      21.790531    -23.3909057     -0.3317652      0.0676780     -0.0003919
      22.119380    -23.3904823     -0.3388444      0.0676812     -0.0004003
      22.453193    -23.3898471     -0.3489034      0.0676861     -0.0004122
      22.792043    -23.3896402     -0.3548715      0.0676877     -0.0004192
      23.136007    -23.3897498     -0.3631443      0.0676868     -0.0004290
      23.485162    -23.3901138     -0.3724790      0.0676839     -0.0004400
      23.839586    -23.3910997     -0.3803616      0.0676761     -0.0004493
      24.199359    -23.3912266     -0.3871427      0.0676751     -0.0004573
      24.564561    -23.3907672     -0.3933647      0.0676786     -0.0004647
      24.935275    -23.3909981     -0.4038726      0.0676767     -0.0004771
      25.311583    -23.3910302     -0.4116988      0.0676764     -0.0004863
      25.693570    -23.3908670     -0.4168926      0.0676776     -0.0004924
      26.081322    -23.3900908     -0.4251055      0.0676836     -0.0005022
      26.474926    -23.3893478     -0.4322710      0.0676893     -0.0005107
      26.874470    -23.3886745     -0.4379036      0.0676945     -0.0005174
      27.280043    -23.3898798     -0.4486704      0.0676850     -0.0005300
      27.691738    -23.3908258     -0.4589735      0.0676776     -0.0005421
      28.109645    -23.3909153     -0.4673093      0.0676768     -0.0005520
      28.533859    -23.3901846     -0.4741200      0.0676824     -0.0005601
      28.964475    -23.3893553     -0.4826914      0.0676888     -0.0005703
      29.401589    -23.3887373     -0.4964205      0.0676935     -0.0005865
      29.845301    -23.3897318     -0.5054926      0.0676856     -0.0005972
      30.295708    -23.3906040     -0.5134703      0.0676788     -0.0006065
      30.752913    -23.3889212     -0.5230890      0.0676918     -0.0006180
      31.217017    -23.3884860     -0.5328943      0.0676951     -0.0006296
      31.688126    -23.3886961     -0.5427308      0.0676933     -0.0006412
      32.166344    -23.3891266     -0.5518495      0.0676898     -0.0006520
      32.651780    -23.3889287     -0.5619796      0.0676913     -0.0006640
      33.144541    -23.3883213     -0.5730670      0.0676959     -0.0006771
      33.644739    -23.3879644     -0.5875980      0.0676985     -0.0006943
      34.152485    -23.3872202     -0.6008838      0.0677041     -0.0007101
      34.667894    -23.3861193     -0.6130452      0.0677126     -0.0007245
      35.191081    -23.3858540     -0.6225010      0.0677145     -0.0007357
      35.722164    -23.3860551     -0.6323883      0.0677128     -0.0007474
      36.261261    -23.3868211     -0.6428196      0.0677067     -0.0007597
      36.808495    -23.3870465     -0.6542227      0.0677048     -0.0007731
      37.363987    -23.3867329     -0.6671697      0.0677071     -0.0007884
      37.927862    -23.3855962     -0.6824744      0.0677157     -0.0008066
      38.500246    -23.3859794     -0.6906432      0.0677126     -0.0008162
      39.081269    -23.3861802     -0.6983173      0.0677109     -0.0008253
      39.671060    -23.3852566     -0.7079682      0.0677180     -0.0008368
      40.269752    -23.3849973     -0.7248218      0.0677198     -0.0008567
      40.877479    -23.3847964     -0.7430753      0.0677210     -0.0008783
      41.494378    -23.3840895     -0.7580618      0.0677263     -0.0008961
      42.120586    -23.3836008     -0.7650228      0.0677300     -0.0009044
      42.756245    -23.3833355     -0.7710279      0.0677320     -0.0009115
      43.401497    -23.3837736     -0.7918953      0.0677282     -0.0009361
      44.056487    -23.3838990     -0.8099421      0.0677270     -0.0009575
      44.721361    -23.3838345     -0.8262256      0.0677272     -0.0009767
      45.396269    -23.3838750     -0.8420033      0.0677266     -0.0009954
      46.081362    -23.3840992     -0.8555713      0.0677246     -0.0010114
      46.776795    -23.3844812     -0.8672615      0.0677214     -0.0010251
      47.482722    -23.3844650     -0.8802054      0.0677213     -0.0010404
      48.199303    -23.3845995     -0.8939112      0.0677200     -0.0010566
      48.926698    -23.3849094     -0.9084741      0.0677173     -0.0010738
      49.665071    -23.3849124     -0.9284058      0.0677169     -0.0010974
      50.414587    -23.3848129     -0.9466157      0.0677174     -0.0011189
      51.175414    -23.3845888     -0.9613392      0.0677188     -0.0011363
      51.947722    -23.3837187     -0.9781699      0.0677253     -0.0011563
      52.731686    -23.3833233     -0.9956651      0.0677280     -0.0011771
      53.527482    -23.3842106     -1.0136200      0.0677207     -0.0011982
      54.335286    -23.3841030     -1.0335283      0.0677211     -0.0012217
      55.155282    -23.3836922     -1.0531217      0.0677239     -0.0012449
      55.987653    -23.3832208     -1.0696118      0.0677272     -0.0012645
      56.832585    -23.3833269     -1.0855186      0.0677260     -0.0012833
      57.690269    -23.3836808     -1.1016879      0.0677229     -0.0013023
      58.560896    -23.3840504     -1.1199281      0.0677196     -0.0013238
      59.444662    -23.3849665     -1.1376551      0.0677121     -0.0013447
      60.341765    -23.3860651     -1.1554766      0.0677031     -0.0013655
      61.252407    -23.3857611     -1.1763361      0.0677049     -0.0013902
      62.176792    -23.3854788     -1.1973393      0.0677066     -0.0014151
      63.115127    -23.3852249     -1.2184995      0.0677081     -0.0014402
      64.067623    -23.3855288     -1.2388600      0.0677052     -0.0014642
      65.034493    -23.3858065     -1.2591188      0.0677025     -0.0014881
      66.015955    -23.3860553     -1.2792436      0.0677000     -0.0015118
      67.012229    -23.3862143     -1.3005062      0.0676982     -0.0015369
      68.023537    -23.3861512     -1.3224552      0.0676981     -0.0015628
      69.050108    -23.3857683     -1.3451978      0.0677005     -0.0015898
      70.092171    -23.3858950     -1.3683740      0.0676989     -0.0016171
      71.149960    -23.3859007     -1.3907673      0.0676982     -0.0016436
      72.223713    -23.3855060     -1.4111948      0.0677007     -0.0016678
      73.313670    -23.3854701     -1.4343323      0.0677003     -0.0016952
      74.420077    -23.3856380     -1.4592564      0.0676982     -0.0017246
      75.543180    -23.3861258     -1.4869097      0.0676936     -0.0017571
      76.683233    -23.3855853     -1.5102742      0.0676971     -0.0017849
      77.840490    -23.3847216     -1.5329210      0.0677031     -0.0018118
      79.015213    -23.3842262     -1.5592467      0.0677061     -0.0018430
      80.207663    -23.3842568     -1.5860963      0.0677050     -0.0018747
      81.418109    -23.3845986     -1.6130400      0.0677014     -0.0019065
      82.646823    -23.3850976     -1.6373157      0.0676967     -0.0019351
      83.894080    -23.3853773     -1.6643151      0.0676936     -0.0019669
      85.160159    -23.3854596     -1.6933902      0.0676920     -0.0020012
      86.445345    -23.3846814     -1.7203507      0.0676971     -0.0020333
      87.749927    -23.3847891     -1.7486566      0.0676953     -0.0020667
      89.074197    -23.3858048     -1.7783376      0.0676863     -0.0021015
      90.418451    -23.3855128     -1.8077358      0.0676875     -0.0021363
      91.782993    -23.3851567     -1.8373340      0.0676891     -0.0021714
      93.168127    -23.3847858     -1.8670861      0.0676909     -0.0022066
      94.574165    -23.3850562     -1.8977448      0.0676876     -0.0022428
      96.001422    -23.3852198     -1.9282329      0.0676851     -0.0022787
      97.450218    -23.3850244     -1.9578884      0.0676854     -0.0023138
      98.920879    -23.3853688     -1.9899477      0.0676814     -0.0023516
      100.413733    -23.3856947     -2.0231683      0.0676775     -0.0023908
      101.929118    -23.3856297     -2.0574816      0.0676766     -0.0024313
      103.467371    -23.3858393     -2.0885797      0.0676736     -0.0024680
      105.028839    -23.3860300     -2.1199227      0.0676708     -0.0025049
      106.613872    -23.3857417     -2.1561696      0.0676714     -0.0025478
      108.222825    -23.3853971     -2.1914974      0.0676725     -0.0025897
      109.856059    -23.3850672     -2.2267776      0.0676735     -0.0026314
      111.513941    -23.3850551     -2.2637513      0.0676719     -0.0026751
      113.196843    -23.3851404     -2.2999261      0.0676695     -0.0027178
      114.905143    -23.3853201     -2.3357375      0.0676664     -0.0027600
      116.639222    -23.3859481     -2.3728772      0.0676597     -0.0028037
      118.399472    -23.3858200     -2.4115274      0.0676588     -0.0028494
      120.186286    -23.3849628     -2.4516629      0.0676635     -0.0028971
      122.000066    -23.3855177     -2.4890074      0.0676572     -0.0029410
      123.841218    -23.3856813     -2.5283364      0.0676539     -0.0029874
      125.710156    -23.3852941     -2.5701772      0.0676547     -0.0030369
      127.607299    -23.3853878     -2.6113357      0.0676518     -0.0030855
      129.533072    -23.3853087     -2.6531219      0.0676502     -0.0031348
      131.487908    -23.3848324     -2.6957984      0.0676515     -0.0031854
      133.472245    -23.3852427     -2.7370839      0.0676460     -0.0032340
      135.486529    -23.3854855     -2.7795487      0.0676417     -0.0032840
      137.531211    -23.3848510     -2.8250537      0.0676440     -0.0033380
      139.606750    -23.3852486     -2.8705797      0.0676382     -0.0033916
      141.713612    -23.3857754     -2.9164233      0.0676314     -0.0034455
      143.852269    -23.3855154     -2.9624093      0.0676306     -0.0034999
      146.023202    -23.3853242     -3.0095657      0.0676292     -0.0035556
      148.226897    -23.3852339     -3.0576152      0.0676269     -0.0036124
      150.463849    -23.3855626     -3.1061192      0.0676213     -0.0036695
      152.734560    -23.3854494     -3.1553304      0.0676190     -0.0037276
      155.039539    -23.3851187     -3.2051960      0.0676183     -0.0037866
      157.379303    -23.3858606     -3.2548171      0.0676092     -0.0038448
      159.754378    -23.3860004     -3.3053804      0.0676047     -0.0039044
      162.165296    -23.3855964     -3.3568789      0.0676043     -0.0039654
      164.612598    -23.3854811     -3.4113127      0.0676014     -0.0040297
      167.096833    -23.3855531     -3.4656310      0.0675970     -0.0040937
      169.618559    -23.3858507     -3.5196267      0.0675908     -0.0041573
      172.178341    -23.3855068     -3.5745166      0.0675895     -0.0042222
      174.776754    -23.3855363     -3.6317153      0.0675850     -0.0042896
      177.414380    -23.3862775     -3.6921806      0.0675746     -0.0043606
      180.091812    -23.3864996     -3.7480172      0.0675686     -0.0044263
      182.809651    -23.3864934     -3.8052161      0.0675642     -0.0044938
      185.568505    -23.3861843     -3.8667920      0.0675618     -0.0045666
      188.368994    -23.3860565     -3.9274879      0.0675579     -0.0046382
      191.211747    -23.3861134     -3.9885100      0.0675525     -0.0047101
      194.097401    -23.3866052     -4.0503032      0.0675435     -0.0047827
      197.026603    -23.3871385     -4.1131148      0.0675341     -0.0048564
      200.000011    -23.3874934     -4.1770353      0.0675259     -0.0049316
      203.018292    -23.3867751     -4.2425788      0.0675258     -0.0050092
      206.082123    -23.3864824     -4.3092117      0.0675222     -0.0050879
      209.192192    -23.3864424     -4.3770075      0.0675164     -0.0051679
      212.349195    -23.3862762     -4.4468824      0.0675114     -0.0052503
      215.553843    -23.3864272     -4.5166344      0.0675037     -0.0053324
      218.806853    -23.3868430     -4.5864745      0.0674940     -0.0054144
      222.108956    -23.3867661     -4.6593964      0.0674876     -0.0055003
      225.460892    -23.3866964     -4.7325078      0.0674811     -0.0055865
      228.863413    -23.3866353     -4.8056963      0.0674744     -0.0056727
      232.317283    -23.3867069     -4.8799286      0.0674664     -0.0057601
      235.823277    -23.3865334     -4.9557740      0.0674601     -0.0058495
      239.382182    -23.3859672     -5.0335227      0.0674565     -0.0059414
      242.994795    -23.3866254     -5.1124761      0.0674431     -0.0060339
      246.661927    -23.3872659     -5.1918354      0.0674297     -0.0061269
      250.384402    -23.3874327     -5.2707231      0.0674199     -0.0062196
      254.163054    -23.3875100     -5.3531797      0.0674103     -0.0063166
      257.998732    -23.3875440     -5.4374118      0.0674007     -0.0064156
      261.892295    -23.3875196     -5.5225437      0.0673912     -0.0065158
      265.844617    -23.3879322     -5.6072005      0.0673783     -0.0066151
      269.856586    -23.3884769     -5.6927623      0.0673642     -0.0067152
      273.929101    -23.3887773     -5.7812705      0.0673514     -0.0068191
      278.063075    -23.3885017     -5.8712862      0.0673427     -0.0069251
      282.259437    -23.3879869     -5.9626885      0.0673356     -0.0070329
      286.519129    -23.3883600     -6.0548291      0.0673213     -0.0071409
      290.843105    -23.3887043     -6.1496575      0.0673067     -0.0072520
      295.232335    -23.3890233     -6.2468625      0.0672918     -0.0073659
      299.687806    -23.3892056     -6.3427570      0.0672780     -0.0074784
      304.210516    -23.3894850     -6.4403446      0.0672630     -0.0075927
      308.801479    -23.3898684     -6.5396650      0.0672468     -0.0077090
      313.461727    -23.3890948     -6.6399716      0.0672392     -0.0078274
      318.192305    -23.3890000     -6.7419579      0.0672259     -0.0079471
      322.994273    -23.3899779     -6.8457586      0.0672038     -0.0080680
      327.868711    -23.3905608     -6.9514530      0.0671843     -0.0081914
      332.816710    -23.3910068     -7.0587411      0.0671654     -0.0083168
      337.839381    -23.3912876     -7.1675543      0.0671473     -0.0084440
      342.937852    -23.3911979     -7.2801127      0.0671313     -0.0085760
      348.113265    -23.3912212     -7.3938294      0.0671139     -0.0087092
      353.366783    -23.3918270     -7.5062802      0.0670920     -0.0088403
      358.699584    -23.3922289     -7.6225098      0.0670709     -0.0089760
      364.112864    -23.3924777     -7.7405643      0.0670503     -0.0091139
      369.607839    -23.3923924     -7.8573272      0.0670322     -0.0092506
      375.185740    -23.3923625     -7.9797724      0.0670126     -0.0093939
      380.847820    -23.3924188     -8.1056470      0.0669913     -0.0095410
      386.595348    -23.3929096     -8.2294195      0.0669668     -0.0096852
      392.429615    -23.3926425     -8.3561988      0.0669472     -0.0098336
      398.351928    -23.3918918     -8.4857022      0.0669306     -0.0099858
      404.363618    -23.3930709     -8.6165499      0.0668986     -0.0101372
      410.466033    -23.3940487     -8.7486466      0.0668675     -0.0102903
      416.660542    -23.3948239     -8.8820168      0.0668374     -0.0104450
      422.948534    -23.3952100     -9.0200661      0.0668091     -0.0106055
      429.331421    -23.3954856     -9.1593791      0.0667809     -0.0107676
      435.810635    -23.3956282     -9.2995722      0.0667533     -0.0109308
      442.387630    -23.3954849     -9.4434575      0.0667267     -0.0110985
      449.063880    -23.3957776     -9.5891572      0.0666961     -0.0112678
      455.840885    -23.3969420     -9.7360493      0.0666580     -0.0114372
      462.720164    -23.3969050     -9.8863050      0.0666281     -0.0116121
      469.703261    -23.3971193    -10.0396278      0.0665951     -0.0117900
      476.791742    -23.3986971    -10.1966654      0.0665505     -0.0119703
      483.987199    -23.3990252    -10.3538793      0.0665149     -0.0121524
      491.291245    -23.3989698    -10.5125511      0.0664814     -0.0123367
      498.705520    -23.3991282    -10.6730324      0.0664453     -0.0125226
      506.231687    -23.3997025    -10.8379524      0.0664046     -0.0127130
      513.871433    -23.4004176    -11.0062625      0.0663615     -0.0129069
      521.626475    -23.4006489    -11.1764072      0.0663212     -0.0131036
      529.498551    -23.4013674    -11.3486389      0.0662760     -0.0133018
      537.489427    -23.4024193    -11.5232439      0.0662271     -0.0135020
      545.600897    -23.4033149    -11.7022365      0.0661778     -0.0137075
      553.834781    -23.4038129    -11.8817156      0.0661307     -0.0139139
      562.192925    -23.4039353    -12.0618401      0.0660857     -0.0141215
      570.677206    -23.4039057    -12.2482230      0.0660397     -0.0143365
      579.289526    -23.4044314    -12.4382608      0.0659878     -0.0145546
      588.031818    -23.4056665    -12.6321440      0.0659288     -0.0147757
      596.906044    -23.4062809    -12.8244966      0.0658741     -0.0149959
      605.914194    -23.4068777    -13.0213177      0.0658177     -0.0152210
      615.058290    -23.4075496    -13.2245678      0.0657582     -0.0154532
      624.340383    -23.4081445    -13.4281216      0.0656984     -0.0156856
      633.762556    -23.4087042    -13.6352664      0.0656370     -0.0159220
      643.326923    -23.4092085    -13.8481130      0.0655736     -0.0161648
      653.035630    -23.4099172    -14.0594604      0.0655082     -0.0164053
      662.890854    -23.4108222    -14.2739572      0.0654395     -0.0166486
      672.894808    -23.4121751    -14.4968808      0.0653640     -0.0169005
      683.049735    -23.4131482    -14.7223269      0.0652897     -0.0171556
      693.357915    -23.4139402    -14.9510097      0.0652148     -0.0174145
      703.821659    -23.4146869    -15.1844595      0.0651377     -0.0176785
      714.443316    -23.4155033    -15.4169957      0.0650593     -0.0179411
      725.225269    -23.4164463    -15.6507707      0.0649785     -0.0182044
      736.169936    -23.4183881    -15.8956336      0.0648856     -0.0184778
      747.279774    -23.4196025    -16.1408916      0.0647968     -0.0187528
      758.557275    -23.4201739    -16.3869665      0.0647115     -0.0190296
      770.004970    -23.4218201    -16.6451949      0.0646128     -0.0193174
      781.625426    -23.4235777    -16.9043797      0.0645117     -0.0196055
      793.421251    -23.4254614    -17.1639934      0.0644082     -0.0198933
      805.395091    -23.4264436    -17.4303697      0.0643078     -0.0201903
      817.549634    -23.4277804    -17.7015011      0.0642016     -0.0204912
      829.887606    -23.4298046    -17.9775220      0.0640872     -0.0207954
      842.411775    -23.4312941    -18.2542543      0.0639751     -0.0211011
      855.124952    -23.4328381    -18.5350611      0.0638595     -0.0214106
      868.029988    -23.4347014    -18.8212521      0.0637381     -0.0217246
      881.129779    -23.4360771    -19.1142902      0.0636161     -0.0220468
      894.427264    -23.4374467    -19.4130038      0.0634903     -0.0223747
      907.925428    -23.4392388    -19.7178829      0.0633572     -0.0227075
      921.627297    -23.4408224    -20.0270560      0.0632223     -0.0230448
      935.535947    -23.4425694    -20.3397864      0.0630828     -0.0233849
      949.654498    -23.4455893    -20.6521875      0.0629325     -0.0237202
      963.986117    -23.4479719    -20.9735410      0.0627813     -0.0240662
      978.534021    -23.4499991    -21.3024045      0.0626275     -0.0244205
      993.301474    -23.4522656    -21.6373369      0.0624674     -0.0247797
      1008.291788    -23.4547310    -21.9765690      0.0623019     -0.0251420
      1023.508326    -23.4573661    -22.3202652      0.0621311     -0.0255075
      1038.954504    -23.4598650    -22.6665746      0.0619579     -0.0258751
      1054.633786    -23.4627337    -23.0193039      0.0617771     -0.0262474
      1070.549691    -23.4660128    -23.3786760      0.0615880     -0.0266243
      1086.705788    -23.4693043    -23.7442336      0.0613936     -0.0270064
      1103.105705    -23.4722328    -24.1156385      0.0611966     -0.0273946
      1119.753118    -23.4745916    -24.4930693      0.0609982     -0.0277896
      1136.651765    -23.4790453    -24.8777147      0.0607791     -0.0281840
      1153.805436    -23.4834771    -25.2685282      0.0605546     -0.0285833
      1171.217980    -23.4871232    -25.6654971      0.0603297     -0.0289900
      1188.893304    -23.4908248    -26.0704327      0.0600977     -0.0294031
      1206.835373    -23.4945177    -26.4815646      0.0598598     -0.0298209
      1225.048213    -23.4980742    -26.8975887      0.0596172     -0.0302424
      1243.535911    -23.5028872    -27.3179083      0.0593609     -0.0306619
      1262.302613    -23.5081391    -27.7453807      0.0590947     -0.0310852
      1281.352532    -23.5129422    -28.1867120      0.0588210     -0.0315220
      1300.689941    -23.5182802    -28.6366235      0.0585357     -0.0319632
      1320.319178    -23.5239287    -29.0937938      0.0582409     -0.0324082
      1340.244648    -23.5293999    -29.5525741      0.0579430     -0.0328528
      1360.470822    -23.5351930    -30.0228674      0.0576330     -0.0333051
      1381.002237    -23.5412627    -30.5035041      0.0573115     -0.0337638
      1401.843499    -23.5475905    -30.9821251      0.0569859     -0.0342164
      1422.999286    -23.5542356    -31.4718035      0.0566480     -0.0346756
      1444.474343    -23.5612113    -31.9728509      0.0562974     -0.0351414
      1466.273489    -23.5691052    -32.4799727      0.0559333     -0.0356060
      1488.401614    -23.5765770    -32.9960847      0.0555625     -0.0360773
      1510.863685    -23.5833493    -33.5219578      0.0551860     -0.0365572
      1533.664739    -23.5927087    -34.0598159      0.0547813     -0.0370337
      1556.809893    -23.6020718    -34.6052222      0.0543677     -0.0375131
      1580.304340    -23.6106601    -35.1567036      0.0539507     -0.0379970
      1604.153351    -23.6202220    -35.7176995      0.0535172     -0.0384811
      1628.362278    -23.6303238    -36.2886136      0.0530693     -0.0389671
      1652.936551    -23.6411019    -36.8713272      0.0526050     -0.0394558
      1677.881684    -23.6524220    -37.4657315      0.0521248     -0.0399473
      1703.203274    -23.6638264    -38.0689837      0.0516335     -0.0404407
      1728.907003    -23.6744586    -38.6763745      0.0511393     -0.0409357
      1754.998637    -23.6872163    -39.2985315      0.0506174     -0.0414277
      1781.484030    -23.7011867    -39.9323161      0.0500754     -0.0419176
      1808.369126    -23.7152834    -40.5704939      0.0495250     -0.0424038
      1835.659955    -23.7300872    -41.2235687      0.0489549     -0.0428924
      1863.362641    -23.7454770    -41.8902269      0.0483668     -0.0433822
      1891.483399    -23.7612274    -42.5655175      0.0477655     -0.0438696
      1920.028539    -23.7777806    -43.2510676      0.0471472     -0.0443533
      1949.004466    -23.7951341    -43.9470302      0.0465120     -0.0448331
      1978.417679    -23.8148256    -44.6594891      0.0458468     -0.0453051
      2008.274780    -23.8343237    -45.3841593      0.0451687     -0.0457785
      2038.582466    -23.8533899    -46.1214484      0.0444781     -0.0462543
      2069.347537    -23.8746382    -46.8652400      0.0437667     -0.0467133
      2100.576897    -23.8967852    -47.6204242      0.0430374     -0.0471657
      2132.277552    -23.9199097    -48.3882690      0.0422887     -0.0476113
      2164.456614    -23.9447020    -49.1695047      0.0415169     -0.0480461
      2197.121303    -23.9706831    -49.9640843      0.0407247     -0.0484720
      2230.278948    -23.9983896    -50.7735906      0.0399083     -0.0488864
      2263.936989    -24.0256927    -51.5950609      0.0390803     -0.0492983
      2298.102976    -24.0537565    -52.4293161      0.0382346     -0.0497012
      2332.784577    -24.0846166    -53.2779505      0.0373613     -0.0500838
      2367.989571    -24.1164348    -54.1356567      0.0364736     -0.0504523
      2403.725859    -24.1493869    -55.0040932      0.0355695     -0.0508062
      2440.001457    -24.1854792    -55.8831561      0.0346416     -0.0511330
      2476.824504    -24.2226914    -56.7825676      0.0336900     -0.0514496
      2514.203264    -24.2609618    -57.6995908      0.0327178     -0.0517537
      2552.146121    -24.3014759    -58.6190453      0.0317348     -0.0520288
      2590.661589    -24.3436535    -59.5548639      0.0307312     -0.0522856
      2629.758310    -24.3874293    -60.5070800      0.0297079     -0.0525237
      2669.445055    -24.4340163    -61.4671316      0.0286695     -0.0527306
      2709.730729    -24.4817697    -62.4429426      0.0276151     -0.0529194
      2750.624370    -24.5307799    -63.4351527      0.0265444     -0.0530893
      2792.135153    -24.5840123    -64.4325804      0.0254597     -0.0532162
      2834.272393    -24.6389792    -65.4434622      0.0243622     -0.0533186
      2877.045543    -24.6957734    -66.4684982      0.0232520     -0.0533956
      2920.464199    -24.7552258    -67.5110547      0.0221246     -0.0534427
      2964.538104    -24.8167222    -68.5655140      0.0209883     -0.0534610
      3009.277147    -24.8810106    -69.6259820      0.0198478     -0.0534433
      3054.691364    -24.9478469    -70.6998579      0.0186983     -0.0533934
      3100.790945    -25.0171491    -71.7873504      0.0175410     -0.0533116
      3147.586234    -25.0908203    -72.8851257      0.0163769     -0.0531849
      3195.087730    -25.1669169    -73.9824476      0.0152214     -0.0530222
      3243.306089    -25.2451628    -75.0896712      0.0140667     -0.0528284
      3292.252132    -25.3268624    -76.2196701      0.0129003     -0.0525984
      3341.936839    -25.4110500    -77.3556141      0.0117406     -0.0523341
      3392.371358    -25.4975499    -78.5005656      0.0105866     -0.0520375
      3443.567005    -25.5883440    -79.6493574      0.0094420     -0.0516960
      3495.535265    -25.6821964    -80.8064014      0.0083059     -0.0513186
      3548.287800    -25.7789228    -81.9730764      0.0071789     -0.0509070
      3601.836445    -25.8804316    -83.1442687      0.0060656     -0.0504501
      3656.193214    -25.9842770    -84.3227037      0.0049670     -0.0499632
      3711.370303    -26.0906198    -85.5068612      0.0038855     -0.0494461
      3767.380091    -26.2023332    -86.6818960      0.0028341     -0.0488826
      3824.235146    -26.3171577    -87.8691914      0.0017967     -0.0482883
      3881.948224    -26.4354911    -89.0691114      0.0007744     -0.0476615
      3940.532273    -26.5564873    -90.2692997     -0.0002209     -0.0470079
      4000.000437    -26.6808484    -91.4741362     -0.0011921     -0.0463248
      4060.366060    -26.8101037    -92.6745678     -0.0021304     -0.0456059
      4121.642685    -26.9428064    -93.8700140     -0.0030347     -0.0448609
      4183.844060    -27.0781362    -95.0723413     -0.0039139     -0.0440950
      4246.984141    -27.2161317    -96.2783954     -0.0047649     -0.0433093
      4311.077094    -27.3584116    -97.4820544     -0.0055814     -0.0424977
      4376.137300    -27.5042263    -98.6883245     -0.0063670     -0.0416655
      4442.179356    -27.6545378    -99.8784551     -0.0071070     -0.0408118];
  case 'BP' % 50 Hz - 8 kHz bandpass filters
    tfinp = [...
      10.000000    -25.6655250    -95.1080000     -0.0046374     -0.0518795
      10.150914    -25.5596949    -95.2685891     -0.0048414     -0.0525021
      10.304106    -25.4369356    -95.3963975     -0.0050291     -0.0532383
      10.459609    -25.2943567    -95.4848809     -0.0051959     -0.0541114
      10.617459    -25.1608328    -95.5306109     -0.0053203     -0.0549455
      10.777692    -25.0287688    -95.6145626     -0.0054835     -0.0557793
      10.940342    -24.8988029    -95.7645548     -0.0057144     -0.0566054
      11.105447    -24.7758541    -95.7692146     -0.0058005     -0.0574119
      11.273044    -24.6563060    -95.7731986     -0.0058850     -0.0582071
      11.443170    -24.5429684    -95.8332703     -0.0060241     -0.0589654
      11.615864    -24.4204281    -95.8996806     -0.0061790     -0.0597960
      11.791163    -24.2945054    -95.9814397     -0.0063558     -0.0606602
      11.969108    -24.1688397    -96.1052910     -0.0065814     -0.0615301
      12.149739    -24.0421771    -96.2176383     -0.0068005     -0.0624207
      12.333096    -23.9142504    -96.3150116     -0.0070091     -0.0633351
      12.519219    -23.7857331    -96.3575240     -0.0071612     -0.0642738
      12.708152    -23.6616541    -96.4242368     -0.0073402     -0.0651901
      12.899936    -23.5394440    -96.5070212     -0.0075397     -0.0661030
      13.094614    -23.4148961    -96.6014781     -0.0077591     -0.0670450
      13.292230    -23.2897432    -96.6885060     -0.0079750     -0.0680060
      13.492829    -23.1638279    -96.7694613     -0.0081890     -0.0689875
      13.696454    -23.0416565    -96.8436120     -0.0083955     -0.0699539
      13.903153    -22.9169908    -96.9275908     -0.0086209     -0.0709526
      14.112971    -22.7897273    -97.0224256     -0.0088673     -0.0719852
      14.325956    -22.6680003    -97.1420949     -0.0091449     -0.0729822
      14.542154    -22.5468298    -97.2362953     -0.0093950     -0.0739921
      14.761616    -22.4266917    -97.2896059     -0.0095957     -0.0750138
      14.984389    -22.3022451    -97.3873355     -0.0098639     -0.0760795
      15.210525    -22.1764678    -97.4887405     -0.0101444     -0.0771714
      15.440073    -22.0507176    -97.5823287     -0.0104202     -0.0782799
      15.673086    -21.9287538    -97.6898292     -0.0107165     -0.0793668
      15.909614    -21.8072544    -97.7889431     -0.0110066     -0.0804659
      16.149713    -21.6863106    -97.8503870     -0.0112485     -0.0815821
      16.393435    -21.5642677    -97.9475768     -0.0115480     -0.0827171
      16.640835    -21.4417578    -98.0608956     -0.0118779     -0.0838687
      16.891968    -21.3223880    -98.1800742     -0.0122191     -0.0850040
      17.146892    -21.2003620    -98.2793592     -0.0125413     -0.0861850
      17.405662    -21.0769814    -98.3698503     -0.0128588     -0.0873978
      17.668338    -20.9598880    -98.4763840     -0.0131980     -0.0885596
      17.934978    -20.8397997    -98.5749786     -0.0135362     -0.0897693
      18.205642    -20.7171528    -98.6680195     -0.0138766     -0.0910235
      18.480391    -20.5986188    -98.7707709     -0.0142327     -0.0922488
      18.759286    -20.4803318    -98.8829602     -0.0146109     -0.0934852
      19.042390    -20.3623609    -99.0049836     -0.0150125     -0.0947318
      19.329766    -20.2420285    -99.0966360     -0.0153756     -0.0960289
      19.621480    -20.1228045    -99.1988627     -0.0157617     -0.0973281
      19.917595    -20.0058632    -99.3174188     -0.0161795     -0.0986141
      20.218180    -19.8840591    -99.4272321     -0.0165996     -0.0999751
      20.523301    -19.7639883    -99.5419504     -0.0170336     -0.1013328
      20.833026    -19.6496835    -99.6672762     -0.0174838     -0.1026371
      21.147426    -19.5315527    -99.7798515     -0.0179276     -0.1040075
      21.466570    -19.4121414    -99.8906246     -0.0183796     -0.1054119
      21.790531    -19.2936688   -100.0030757     -0.0188417     -0.1068228
      22.119380    -19.1779599   -100.1098973     -0.0192961     -0.1082196
      22.453193    -19.0625288   -100.2170770     -0.0197594     -0.1096306
      22.792043    -18.9467915   -100.3317878     -0.0202468     -0.1110608
      23.136007    -18.8312830   -100.4505462     -0.0207511     -0.1125049
      23.485162    -18.7153266   -100.5705916     -0.0212688     -0.1139726
      23.839586    -18.5995059   -100.6800955     -0.0217751     -0.1154611
      24.199359    -18.4842123   -100.7965462     -0.0223038     -0.1169588
      24.564561    -18.3689802   -100.9180914     -0.0228530     -0.1184726
      24.935275    -18.2556233   -101.0326959     -0.0233933     -0.1199823
      25.311583    -18.1412286   -101.1500422     -0.0239523     -0.1215241
      25.693570    -18.0257590   -101.2701396     -0.0245310     -0.1230993
      26.081322    -17.9123229   -101.3857690     -0.0251051     -0.1246671
      26.474926    -17.7989347   -101.5011599     -0.0256893     -0.1262537
      26.874470    -17.6860000   -101.6159267     -0.0262816     -0.1278536
      27.280043    -17.5744552   -101.7361357     -0.0268930     -0.1294500
      27.691738    -17.4632625   -101.8571224     -0.0275162     -0.1310600
      28.109645    -17.3533869   -101.9773521     -0.0281450     -0.1326696
      28.533859    -17.2402530   -102.0961315     -0.0287926     -0.1343496
      28.964475    -17.1267449   -102.2152483     -0.0294542     -0.1360558
      29.401589    -17.0157140   -102.3342880     -0.0301194     -0.1377439
      29.845301    -16.9070111   -102.4481931     -0.0307759     -0.1394177
      30.295708    -16.7979661   -102.5621776     -0.0314455     -0.1411167
      30.752913    -16.6870956   -102.6804826     -0.0321445     -0.1428635
      31.217017    -16.5774483   -102.7988017     -0.0328516     -0.1446108
      31.688126    -16.4681025   -102.9169199     -0.0335696     -0.1463739
      32.166344    -16.3605928   -103.0296886     -0.0342793     -0.1481297
      32.651780    -16.2521727   -103.1434696     -0.0350077     -0.1499211
      33.144541    -16.1429018   -103.2582243     -0.0357548     -0.1517478
      33.644739    -16.0365773   -103.3705134     -0.0364962     -0.1535455
      34.152485    -15.9293688   -103.4833150     -0.0372554     -0.1553794
      34.667894    -15.8211941   -103.5967499     -0.0380337     -0.1572516
      35.191081    -15.7151467   -103.7018784     -0.0387929     -0.1591124
      35.722164    -15.6089608   -103.8067997     -0.0395650     -0.1609973
      36.261261    -15.5029014   -103.9112621     -0.0403481     -0.1629019
      36.808495    -15.3980591   -104.0245996     -0.0411642     -0.1647990
      37.363987    -15.2923690   -104.1355458     -0.0419911     -0.1667356
      37.927862    -15.1858642   -104.2401601     -0.0428172     -0.1687147
      38.500246    -15.0828160   -104.3427229     -0.0436338     -0.1706504
      39.081269    -14.9792830   -104.4446745     -0.0444642     -0.1726178
      39.671060    -14.8744136   -104.5447592     -0.0453094     -0.1746357
      40.269752    -14.7712910   -104.6460846     -0.0461630     -0.1766400
      40.877479    -14.6679226   -104.7460611     -0.0470275     -0.1786730
      41.494378    -14.5639469   -104.8380231     -0.0478840     -0.1807480
      42.120586    -14.4614002   -104.9316750     -0.0487516     -0.1828152
      42.756245    -14.3588525   -105.0265472     -0.0496368     -0.1849044
      43.401497    -14.2558370   -105.1210224     -0.0505375     -0.1870273
      44.056487    -14.1544293   -105.2127446     -0.0514338     -0.1891416
      44.721361    -14.0535044   -105.3025567     -0.0523348     -0.1912703
      45.396269    -13.9513484   -105.3859599     -0.0532356     -0.1934559
      46.081362    -13.8493316   -105.4699038     -0.0541513     -0.1956623
      46.776795    -13.7472268   -105.5544677     -0.0550837     -0.1978948
      47.482722    -13.6482292   -105.6315671     -0.0559845     -0.2000881
      48.199303    -13.5488072   -105.7095429     -0.0569043     -0.2023143
      48.926698    -13.4491125   -105.7883658     -0.0578427     -0.2045704
      49.665071    -13.3479628   -105.8678634     -0.0588074     -0.2068852
      50.414587    -13.2472230   -105.9405279     -0.0597588     -0.2092230
      51.175414    -13.1479756   -106.0018276     -0.0606719     -0.2115626
      51.947722    -13.0488897   -106.0694411     -0.0616205     -0.2139173
      52.731686    -12.9492716   -106.1372901     -0.0625875     -0.2163109
      53.527482    -12.8496222   -106.2024938     -0.0635586     -0.2187346
      54.335286    -12.7497041   -106.2640153     -0.0645315     -0.2211962
      55.155282    -12.6495687   -106.3239436     -0.0655138     -0.2236926
      55.987653    -12.5508628   -106.3810105     -0.0664878     -0.2261830
      56.832585    -12.4521844   -106.4274751     -0.0674330     -0.2287127
      57.690269    -12.3530852   -106.4729818     -0.0683905     -0.2312828
      58.560896    -12.2546361   -106.5346282     -0.0694217     -0.2338446
      59.444662    -12.1567495   -106.5827478     -0.0704071     -0.2364358
      60.341765    -12.0586626   -106.6232265     -0.0713756     -0.2390705
      61.252407    -11.9598516   -106.6651832     -0.0723692     -0.2417528
      62.176792    -11.8605246   -106.7075281     -0.0733822     -0.2444791
      63.115127    -11.7605605   -106.7501098     -0.0744154     -0.2472537
      64.067623    -11.6630415   -106.7832753     -0.0754004     -0.2500018
      65.034493    -11.5648767   -106.8128286     -0.0763878     -0.2528039
      66.015955    -11.4661188   -106.8384052     -0.0773754     -0.2556601
      67.012229    -11.3673969   -106.8669076     -0.0783884     -0.2585435
      68.023537    -11.2684236   -106.8950443     -0.0794152     -0.2614674
      69.050108    -11.1696325   -106.9222362     -0.0804491     -0.2644201
      70.092171    -11.0708264   -106.9451720     -0.0814765     -0.2674126
      71.149960    -10.9712999   -106.9659878     -0.0825137     -0.2704644
      72.223713    -10.8713803   -106.9835527     -0.0835523     -0.2735681
      73.313670    -10.7724331   -106.9950855     -0.0845652     -0.2766853
      74.420077    -10.6731149   -107.0058287     -0.0855902     -0.2798511
      75.543180    -10.5736891   -107.0186157     -0.0866387     -0.2830536
      76.683233    -10.4735905   -107.0290419     -0.0876951     -0.2863185
      77.840490    -10.3728576   -107.0382393     -0.0887645     -0.2896441
      79.015213    -10.2730704   -107.0458128     -0.0898289     -0.2929790
      80.207663    -10.1727307   -107.0482943     -0.0908854     -0.2963792
      81.418109    -10.0717312   -107.0480477     -0.0919471     -0.2998460
      82.646823     -9.9718934   -107.0484818     -0.0930124     -0.3033117
      83.894080     -9.8709592   -107.0507807     -0.0941119     -0.3068531
      85.160159     -9.7689122   -107.0540781     -0.0952419     -0.3104740
      86.445345     -9.6680956   -107.0459849     -0.0963095     -0.3141122
      87.749927     -9.5666607   -107.0381167     -0.0973971     -0.3178153
      89.074197     -9.4646070   -107.0304805     -0.0985054     -0.3215846
      90.418451     -9.3631090   -107.0230497     -0.0996210     -0.3253774
      91.782993     -9.2605462   -107.0145492     -0.1007555     -0.3292572
      93.168127     -9.1569587   -107.0045960     -0.1019064     -0.3332251
      94.574165     -9.0547080   -106.9933882     -0.1030471     -0.3371912
      96.001422     -8.9514545   -106.9814597     -0.1042084     -0.3412452
      97.450218     -8.8469613   -106.9685833     -0.1053920     -0.3453989
      98.920879     -8.7435978   -106.9541306     -0.1065655     -0.3495607
      100.413733     -8.6395516   -106.9381676     -0.1077511     -0.3538032
      101.929118     -8.5345934   -106.9193920     -0.1089437     -0.3581401
      103.467371     -8.4303144   -106.8990022     -0.1101305     -0.3625049
      105.028839     -8.3255119   -106.8795525     -0.1113428     -0.3669432
      106.613872     -8.2201022   -106.8667054     -0.1126190     -0.3714487
      108.222825     -8.1142861   -106.8477277     -0.1138748     -0.3760393
      109.856059     -8.0076784   -106.8264716     -0.1151399     -0.3807258
      111.513941     -7.9008601   -106.8121014     -0.1164679     -0.3854661
      113.196843     -7.7940596   -106.7951916     -0.1177937     -0.3902698
      114.905143     -7.6868762   -106.7762723     -0.1191258     -0.3951549
      116.639222     -7.5796173   -106.7548191     -0.1204561     -0.4001099
      118.399472     -7.4712310   -106.7339267     -0.1218209     -0.4051784
      120.186286     -7.3616756   -106.7135589     -0.1232213     -0.4103651
      122.000066     -7.2540761   -106.6958034     -0.1246285     -0.4155189
      123.841218     -7.1456435   -106.6794362     -0.1260739     -0.4207747
      125.710156     -7.0364821   -106.6648575     -0.1275599     -0.4261287
      127.607299     -6.9271108   -106.6460113     -0.1290343     -0.4315708
      129.533072     -6.8169761   -106.6266270     -0.1305330     -0.4371220
      131.487908     -6.7064564   -106.6072724     -0.1320549     -0.4427642
      133.472245     -6.5958978   -106.5946112     -0.1336474     -0.4484655
      135.486529     -6.4847973   -106.5830687     -0.1352764     -0.4542659
      137.531211     -6.3739459   -106.5709930     -0.1369169     -0.4601293
      139.606750     -6.2630415   -106.5569941     -0.1385624     -0.4660759
      141.713612     -6.1512965   -106.5430661     -0.1402418     -0.4721449
      143.852269     -6.0389642   -106.5320315     -0.1419752     -0.4783180
      146.023202     -5.9265036   -106.5215953     -0.1437371     -0.4845775
      148.226897     -5.8135295   -106.5120008     -0.1455366     -0.4909458
      150.463849     -5.7015184   -106.5061203     -0.1473745     -0.4973330
      152.734560     -5.5879230   -106.5010113     -0.1492697     -0.5038933
      155.039539     -5.4729405   -106.4966754     -0.1512202     -0.5106195
      157.379303     -5.3597847   -106.4962173     -0.1531989     -0.5173164
      159.754378     -5.2454056   -106.4941064     -0.1552104     -0.5241794
      162.165296     -5.1297316   -106.4904985     -0.1572577     -0.5312167
      164.612598     -5.0151869   -106.4951815     -0.1593893     -0.5382555
      167.096833     -4.8997629   -106.5005051     -0.1615722     -0.5454409
      169.618559     -4.7835947   -106.5065281     -0.1638057     -0.5527676
      172.178341     -4.6680607   -106.5147406     -0.1660794     -0.5601455
      174.776754     -4.5519375   -106.5244798     -0.1684112     -0.5676558
      177.414380     -4.4355998   -106.5363509     -0.1708012     -0.5752747
      180.091812     -4.3195121   -106.5510018     -0.1732484     -0.5829706
      182.809651     -4.2025730   -106.5662190     -0.1757535     -0.5908256
      185.568505     -4.0850691   -106.5813548     -0.1783055     -0.5988256
      188.368994     -3.9683008   -106.6004387     -0.1809209     -0.6068700
      191.211747     -3.8508941   -106.6210685     -0.1836044     -0.6150627
      194.097401     -3.7329386   -106.6424687     -0.1863477     -0.6234027
      197.026603     -3.6154231   -106.6693680     -0.1891826     -0.6318056
      200.000011     -3.4972372   -106.6996183     -0.1921125     -0.6403598
      203.018292     -3.3779358   -106.7338808     -0.1951576     -0.6490993
      206.082123     -3.2588616   -106.7670175     -0.1982319     -0.6579444
      209.192192     -3.1394260   -106.8005473     -0.2013669     -0.6669363
      212.349195     -3.0202716   -106.8445392     -0.2046674     -0.6759915
      215.553843     -2.9008379   -106.8868451     -0.2080071     -0.6851974
      218.806853     -2.7808989   -106.9278162     -0.2113960     -0.6945736
      222.108956     -2.6621187   -106.9683473     -0.2148047     -0.7039854
      225.460892     -2.5420784   -107.0177036     -0.2184088     -0.7135943
      228.863413     -2.4208265   -107.0770561     -0.2222286     -0.7233958
      232.317283     -2.3006978   -107.1328748     -0.2260379     -0.7332503
      235.823277     -2.1799177   -107.1903821     -0.2299491     -0.7432871
      239.382182     -2.0587712   -107.2505242     -0.2339698     -0.7534816
      242.994795     -1.9383094   -107.3172118     -0.2381263     -0.7637275
      246.661927     -1.8174484   -107.3848973     -0.2423776     -0.7741430
      250.384402     -1.6968911   -107.4516421     -0.2466795     -0.7846760
      254.163054     -1.5753175   -107.5240243     -0.2511615     -0.7954194
      257.998732     -1.4529310   -107.6003817     -0.2558003     -0.8063662
      261.892295     -1.3311385   -107.6827996     -0.2605884     -0.8173786
      265.844617     -1.2097439   -107.7690176     -0.2655030     -0.8284840
      269.856586     -1.0877839   -107.8569949     -0.2705471     -0.8397845
      273.929101     -0.9656665   -107.9434993     -0.2756632     -0.8512594
      278.063075     -0.8432172   -108.0369736     -0.2809850     -0.8628878
      282.259437     -0.7200822   -108.1352100     -0.2864969     -0.8747176
      286.519129     -0.5977660   -108.2358181     -0.2921172     -0.8866111
      290.843105     -0.4747282   -108.3390405     -0.2979042     -0.8987243
      295.232335     -0.3508174   -108.4447977     -0.3038666     -0.9110779
      299.687806     -0.2290329   -108.5567556     -0.3099620     -0.9233381
      304.210516     -0.1064057   -108.6711328     -0.3162378     -0.9358368
      308.801479      0.0170200   -108.7880074     -0.3226993     -0.9485736
      313.461727      0.1403526   -108.9102402     -0.3293659     -0.9614382
      318.192305      0.2637761   -109.0342056     -0.3361886     -0.9744723
      322.994273      0.3866780   -109.1596127     -0.3431418     -0.9876100
      327.868711      0.5097110   -109.2950138     -0.3504031     -1.0008735
      332.816710      0.6332162   -109.4335209     -0.3578743     -1.0143445
      337.839381      0.7563477   -109.5739390     -0.3655039     -1.0279336
      342.937852      0.8799918   -109.7167383     -0.3733416     -1.0417437
      348.113265      1.0043643   -109.8636820     -0.3814349     -1.0557929
      353.366783      1.1283656   -110.0182351     -0.3898069     -1.0699261
      358.699584      1.2516079   -110.1767792     -0.3983786     -1.0841170
      364.112864      1.3754532   -110.3382363     -0.4071967     -1.0985422
      369.607839      1.5000368   -110.5016740     -0.4162565     -1.1132295
      375.185740      1.6238974   -110.6727762     -0.4256052     -1.1279519
      380.847820      1.7479522   -110.8491528     -0.4352476     -1.1428429
      386.595348      1.8717215   -111.0274731     -0.4450997     -1.1578647
      392.429615      1.9959426   -111.2085206     -0.4552202     -1.1731103
      398.351928      2.1208270   -111.3926001     -0.4656337     -1.1886091
      404.363618      2.2435698   -111.5858908     -0.4763247     -1.2039248
      410.466033      2.3672296   -111.7804254     -0.4872980     -1.2195401
      416.660542      2.4918276   -111.9762356     -0.4985637     -1.2354637
      422.948534      2.6152957   -112.1836058     -0.5102335     -1.2513124
      429.331421      2.7394062   -112.3944831     -0.5222447     -1.2674069
      435.810635      2.8638944   -112.6086473     -0.5345856     -1.2837133
      442.387630      2.9886941   -112.8248130     -0.5472314     -1.3002356
      449.063880      3.1135479   -113.0472630     -0.5602714     -1.3168953
      455.840885      3.2372837   -113.2788831     -0.5737052     -1.3334812
      462.720164      3.3611574   -113.5140026     -0.5874915     -1.3502354
      469.703261      3.4852732   -113.7524814     -0.6016423     -1.3671756
      476.791742      3.6080587   -113.9940680     -0.6160490     -1.3840542
      483.987199      3.7319263   -114.2427029     -0.6309837     -1.4012084
      491.291245      3.8566765   -114.4963964     -0.6463992     -1.4186301
      498.705520      3.9805169   -114.7543898     -0.6621543     -1.4360343
      506.231687      4.1044694   -115.0184996     -0.6783789     -1.4535627
      513.871433      4.2290074   -115.2878730     -0.6951004     -1.4713022
      521.626475      4.3527081   -115.5624293     -0.7122139     -1.4890100
      529.498551      4.4760515   -115.8414302     -0.7297456     -1.5067699
      537.489427      4.5996095   -116.1252796     -0.7477631     -1.5246713
      545.600897      4.7230875   -116.4206124     -0.7664305     -1.5425706
      553.834781      4.8470501   -116.7183513     -0.7855679     -1.5606825
      562.192925      4.9715991   -117.0186732     -0.8052012     -1.5790239
      570.677206      5.0943760   -117.3292192     -0.8253319     -1.5970523
      579.289526      5.2174735   -117.6441366     -0.8460024     -1.6152216
      588.031818      5.3405868   -117.9632818     -0.8671910     -1.6334737
      596.906044      5.4639855   -118.2833222     -0.8888399     -1.6519067
      605.914194      5.5873528   -118.6125452     -0.9111671     -1.6703284
      615.058290      5.7096172   -118.9548706     -0.9341882     -1.6884553
      624.340383      5.8327885   -119.2976443     -0.9577585     -1.7068699
      633.762556      5.9563157   -119.6455750     -0.9819721     -1.7253866
      643.326923      6.0785612   -120.0009931     -1.0067256     -1.7436301
      653.035630      6.2005907   -120.3654303     -1.0321958     -1.7617698
      662.890854      6.3231227   -120.7363900     -1.0584068     -1.7799839
      672.894808      6.4451003   -121.1124280     -1.0851993     -1.7980734
      683.049735      6.5671467   -121.4942844     -1.1126837     -1.8161414
      693.357915      6.6896857   -121.8822004     -1.1409373     -1.8342623
      703.821659      6.8110033   -122.2773139     -1.1697842     -1.8520388
      714.443316      6.9321256   -122.6775884     -1.1993020     -1.8697131
      725.225269      7.0535561   -123.0836951     -1.2295945     -1.8873679
      736.169936      7.1740661   -123.4998928     -1.2606414     -1.9046292
      747.279774      7.2949623   -123.9228420     -1.2925324     -1.9218356
      758.557275      7.4164312   -124.3525853     -1.3253155     -1.9390148
      770.004970      7.5361718   -124.7891499     -1.3586526     -1.9556351
      781.625426      7.6566033   -125.2326885     -1.3929307     -1.9722154
      793.421251      7.7775450   -125.6833744     -1.4281485     -1.9886964
      805.395091      7.8969719   -126.1440925     -1.4640860     -2.0045210
      817.549634      8.0167501   -126.6111572     -1.5009335     -2.0201867
      829.887606      8.1364593   -127.0838709     -1.5386093     -2.0355970
      842.411775      8.2557482   -127.5664718     -1.5772131     -2.0505340
      855.124952      8.3753092   -128.0562520     -1.6167863     -2.0652096
      868.029988      8.4941646   -128.5520985     -1.6571194     -2.0794007
      881.129779      8.6125086   -129.0570314     -1.6983630     -2.0930404
      894.427264      8.7311741   -129.5698889     -1.7406485     -2.1063354
      907.925428      8.8490782   -130.0899518     -1.7837446     -2.1190190
      921.627297      8.9669133   -130.6194123     -1.8278797     -2.1311621
      935.535947      9.0850131   -131.1576562     -1.8731151     -2.1428354
      949.654498      9.2011904   -131.7046697     -1.9189838     -2.1534668
      963.986117      9.3177217   -132.2590804     -1.9659303     -2.1636314
      978.534021      9.4349140   -132.8213255     -2.0140588     -2.1733624
      993.301474      9.5512362   -133.3918954     -2.0630459     -2.1822281
      1008.291788      9.6671335   -133.9722609     -2.1130517     -2.1902501
      1023.508326      9.7829555   -134.5623539     -2.1641632     -2.1974797
      1038.954504      9.8979444   -135.1611597     -2.2161562     -2.2037240
      1054.633786     10.0132462   -135.7673902     -2.2692733     -2.2092862
      1070.549691     10.1287152   -136.3809870     -2.3234863     -2.2140968
      1086.705788     10.2419463   -137.0078757     -2.3783756     -2.2172600
      1103.105705     10.3555017   -137.6437214     -2.4344547     -2.2195584
      1119.753118     10.4690845   -138.2880184     -2.4916296     -2.2208959
      1136.651765     10.5813842   -138.9462811     -2.5497334     -2.2206499
      1153.805436     10.6939424   -139.6131592     -2.6089980     -2.2193983
      1171.217980     10.8062053   -140.2859998     -2.6691567     -2.2170786
      1188.893304     10.9179517   -140.9663941     -2.7301957     -2.2135216
      1206.835373     11.0297874   -141.6581634     -2.7924452     -2.2086530
      1225.048213     11.1403213   -142.3657198     -2.8556161     -2.2018434
      1243.535911     11.2503179   -143.0814988     -2.9196402     -2.1936024
      1262.302613     11.3605538   -143.8058612     -2.9847803     -2.1840601
      1281.352532     11.4699568   -144.5354301     -3.0505301     -2.1730771
      1300.689941     11.5784699   -145.2816238     -3.1172737     -2.1599821
      1320.319178     11.6868986   -146.0415383     -3.1851607     -2.1450599
      1340.244648     11.7944395   -146.8068643     -3.2535626     -2.1285150
      1360.470822     11.9016840   -147.5829526     -3.3228698     -2.1101444
      1381.002237     12.0089647   -148.3701982     -3.3932008     -2.0899448
      1401.843499     12.1143855   -149.1709800     -3.4638650     -2.0672566
      1422.999286     12.2198783   -149.9830867     -3.5354973     -2.0426121
      1444.474343     12.3253869   -150.8066606     -3.6080540     -2.0159228
      1466.273489     12.4291341   -151.6380061     -3.6806248     -1.9869523
      1488.401614     12.5330857   -152.4818595     -3.7541498     -1.9557977
      1510.863685     12.6369189   -153.3387015     -3.8284713     -1.9222814
      1533.664739     12.7387662   -154.2042083     -3.9025640     -1.8862195
      1556.809893     12.8405377   -155.0834021     -3.9773776     -1.8476379
      1580.304340     12.9416725   -155.9781852     -4.0526596     -1.8062095
      1604.153351     13.0419786   -156.8805578     -4.1279991     -1.7623967
      1628.362278     13.1422409   -157.7961522     -4.2038811     -1.7159015
      1652.936551     13.2413134   -158.7283924     -4.2797807     -1.6661741
      1677.881684     13.3393708   -159.6702133     -4.3554835     -1.6137169
      1703.203274     13.4373325   -160.6240881     -4.4314423     -1.5584634
      1728.907003     13.5340658   -161.5902604     -4.5070066     -1.5001325
      1754.998637     13.6301332   -162.5687880     -4.5823707     -1.4387698
      1781.484030     13.7261505   -163.5606261     -4.6577948     -1.3743424
      1808.369126     13.8207866   -164.5653042     -4.7324585     -1.3066200
      1835.659955     13.9143025   -165.5843257     -4.8064178     -1.2354800
      1863.362641     14.0072656   -166.6178768     -4.8798709     -1.1609404
      1891.483399     14.0989560   -167.6615822     -4.9522097     -1.0832352
      1920.028539     14.1905164   -168.7171754     -5.0240065     -1.0023295
      1949.004466     14.2819842   -169.7849364     -5.0951864     -0.9181517
      1978.417679     14.3711374   -170.8723389     -5.1644304     -0.8297652
      2008.274780     14.4600505   -171.9733087     -5.2327119     -0.7378957
      2038.582466     14.5483673   -173.0870192     -5.2996788     -0.6425499
      2069.347537     14.6351886   -174.2141461     -5.3646490     -0.5435830
      2100.576897     14.7214325   -175.3566954     -5.4280510     -0.4408600
      2132.277552     14.8061786   -176.5143440     -5.4891458     -0.3343512
      2164.456614     14.8898360   -177.6824970     -5.5479999     -0.2245283
      2197.121303     14.9732231   -178.8663930     -5.6050068     -0.1109105
      2230.278948     15.0552974   -180.0672226     -5.6593241      0.0066398
      2263.936989     15.1366624   -181.2810183     -5.7111631      0.1277114
      2298.102976     15.2177589   -182.5111340     -5.7606396      0.2526365
      2332.784577     15.2969260   -183.7588857     -5.8064545      0.3814795
      2367.989571     15.3750850   -185.0208747     -5.8490402      0.5138721
      2403.725859     15.4528541   -186.2987203     -5.8886142      0.6499759
      2440.001457     15.5286460   -187.5898175     -5.9239406      0.7893509
      2476.824504     15.6037495   -188.8979265     -5.9556510      0.9324092
      2514.203264     15.6785335   -190.2237563     -5.9837806      1.0792136
      2552.146121     15.7511273   -191.5642044     -6.0068914      1.2291273
      2590.661589     15.8232535   -192.9184633     -6.0259959      1.3821812
      2629.758310     15.8950392   -194.2872926     -6.0409790      1.5383995
      2669.445055     15.9642800   -195.6785357     -6.0498826      1.6981001
      2709.730729     16.0329288   -197.0870829     -6.0539719      1.8609510
      2750.624370     16.1006552   -198.5123089     -6.0528249      2.0266929
      2792.135153     16.1667489   -199.9511317     -6.0458574      2.1946743
      2834.272393     16.2321628   -201.4099976     -6.0332893      2.3656347
      2877.045543     16.2961809   -202.8893587     -6.0143702      2.5392532
      2920.464199     16.3589759   -204.3831455     -5.9892752      2.7147316
      2964.538104     16.4212574   -205.8976561     -5.9580009      2.8927486
      3009.277147     16.4820558   -207.4343627     -5.9195725      3.0729205
      3054.691364     16.5411763   -208.9832624     -5.8741952      3.2538767
      3100.790945     16.5998608   -210.5504099     -5.8222135      3.4364493
      3147.586234     16.6576455   -212.1353652     -5.7631500      3.6201761
      3195.087730     16.7131477   -213.7429743     -5.6955982      3.8046669
      3243.306089     16.7676739   -215.3732296     -5.6202228      3.9901336
      3292.252132     16.8208862   -217.0222379     -5.5368895      4.1757158
      3341.936839     16.8730189   -218.6889750     -5.4456793      4.3610916
      3392.371358     16.9245438   -220.3764678     -5.3465137      4.5464548
      3443.567005     16.9743800   -222.0909110     -5.2380665      4.7314471
      3495.535265     17.0232371   -223.8224540     -5.1214333      4.9151323
      3548.287800     17.0713371   -225.5726056     -4.9965232      5.0974038
      3601.836445     17.1168818   -227.3502268     -4.8614196      5.2775452
      3656.193214     17.1615873   -229.1506649     -4.7174254      5.4556857
      3711.370303     17.2052026   -230.9736379     -4.5643462      5.6311985
      3767.380091     17.2468749   -232.8198760     -4.4016215      5.8031013
      3824.235146     17.2880221   -234.6912888     -4.2297548      5.9719715
      3881.948224     17.3282975   -236.5872466     -4.0485879      6.1370330
      3940.532273     17.3666612   -238.5079247     -3.8576260      6.2970284
      4000.000437     17.4039249   -240.4548200     -3.6571250      6.4520703
      4060.366060     17.4390450   -242.4265357     -3.4468777      6.6007129
      4121.642685     17.4729304   -244.4270612     -3.2269192      6.7432713
      4183.844060     17.5062705   -246.4572418     -2.9974906      6.8797118
      4246.984141     17.5383685   -248.5172226     -2.7584321      7.0088655
      4311.077094     17.5682622   -250.6009860     -2.5103871      7.1290222
      4376.137300     17.5970270   -252.7148592     -2.2531694      7.2407073
      4442.179356     17.6239943   -254.8693507     -1.9855236      7.3430575
      4509.218079     17.6494192   -257.0497495     -1.7097097      7.4350138
      4577.268511     17.6738320   -259.2597058     -1.4257347      7.5165092
      4646.345919     17.6962157   -261.5084449     -1.1326203      7.5861880
      4716.465803     17.7171070   -263.7929323     -0.8313222      7.6436689
      4787.643894     17.7367826   -266.1134107     -0.5223385      7.6884568
      4859.896162     17.7536691   -268.4719203     -0.2058995      7.7184302
      4933.238818     17.7691244   -270.8683776      0.1172264      7.7340387
      5007.688318     17.7829837   -273.3035554      0.4464445      7.7344047
      5083.261366     17.7940359   -275.7834502      0.7816790      7.7176580
      5159.974917     17.8034043   -278.3056243      1.1217544      7.6840664
      5237.846183     17.8104606   -280.8722917      1.4659259      7.6323214
      5316.892636     17.8144908   -283.4872346      1.8134541      7.5610006
      5397.132010     17.8163413   -286.1472032      2.1628559      7.4702856
      5478.582310     17.8146443   -288.8556531      2.5129479      7.3582998
      5561.261809     17.8094932   -291.6218307      2.8634338      7.2241652
      5645.189058     17.8020054   -294.4377332      3.2121076      7.0686741
      5730.382886     17.7904363   -297.3065843      3.5571278      6.8898661
      5816.862410     17.7744645   -300.2353208      3.8973383      6.6868125
      5904.647030     17.7553780   -303.2187505      4.2307785      6.4606916
      5993.756444     17.7305196   -306.2652203      4.5551034      6.2089177
      6084.210644     17.6999500   -309.3745351      4.8680136      5.9317899
      6176.029925     17.6651185   -312.5428777      5.1676571      5.6310420
      6269.234888     17.6226970   -315.7820994      5.4508957      5.3040789
      6363.846444     17.5736176   -319.0870871      5.7152297      4.9529422
      6459.885822     17.5188716   -322.4554905      5.9586949      4.5796262
      6557.374569     17.4524441   -325.8973132      6.1754754      4.1815335
      6656.334557     17.3789513   -329.4030284      6.3655238      3.7641062
      6756.787991     17.2981752   -332.9738198      6.5266223      3.3292375
      6858.757409     17.2013419   -336.6230858      6.6507307      2.8748472
      6962.265688     17.0960085   -340.3345763      6.7406366      2.4089090
      7067.336052     16.9802601   -344.1104977      6.7935056      1.9338353
      7173.992075     16.8466978   -347.9481137      6.8022956      1.4523122
      7282.257688     16.7021789   -351.8448751      6.7716553      0.9703972
      7392.157181     16.5425687   -355.8016072      6.6982516      0.4916999
      7503.715211     16.3638253   -359.8014553      6.5794358      0.0227995
      7616.956809     16.1723071   -363.8492963      6.4214708     -0.4320632
      7731.907382     15.9614608   -367.9355812      6.2214865     -0.8672409
      7848.592720     15.7299192   -372.0506072      5.9816176     -1.2769551
      7967.039003     15.4839969   -376.2012362      5.7095410     -1.6589086
      8087.272807     15.2158438   -380.3552851      5.4049087     -2.0052674
      8209.321109     14.9278009   -384.5151254      5.0741194     -2.3140274
      8333.211290     14.6244350   -388.6950887      4.7240458     -2.5858132
      8458.971148     14.2965046   -392.8435098      4.3569708     -2.8125611
      8586.628899     13.9513024   -396.9770968      3.9814797     -2.9977656
      8716.213185     13.5911488   -401.1110674      3.6025004     -3.1438870
      8847.753080     13.2060120   -405.1697980      3.2247402     -3.2439103
      8981.278096     12.8059469   -409.2085564      2.8537441     -3.3070927
      9116.818192     12.3911309   -413.2306269      2.4928146     -3.3359302
      9254.403778     11.9546141   -417.1479608      2.1483626     -3.3269661
      9394.065725     11.5058102   -421.0373271      1.8211718     -3.2905347
      9535.835366     11.0433802   -424.8784175      1.5138685     -3.2285966
      9679.744510     10.5646097   -428.6083485      1.2308786     -3.1421795
      9825.825445     10.0750962   -432.3110728      0.9691981     -3.0389265
      9974.110946      9.5738221   -435.9555828      0.7306584     -2.9208627];
  case 'U' % Unfiltered
    tfinp = [...
      100.000000    -23.3861313     -0.1979361      0.0677159     -0.0002339
      101.509140    -23.3861313     -0.2031420      0.0677159     -0.0002401
      103.041056    -23.3858179     -0.2088265      0.0677183     -0.0002468
      104.596090    -23.3852907     -0.2145165      0.0677224     -0.0002536
      106.174592    -23.3848689     -0.2182648      0.0677257     -0.0002580
      107.776916    -23.3841739     -0.2242174      0.0677311     -0.0002651
      109.403421    -23.3834060     -0.2306945      0.0677370     -0.0002727
      111.054472    -23.3827016     -0.2365140      0.0677425     -0.0002796
      112.730440    -23.3829922     -0.2419003      0.0677402     -0.0002860
      114.431700    -23.3831735     -0.2472293      0.0677388     -0.0002923
      116.158635    -23.3830134     -0.2524319      0.0677400     -0.0002984
      117.911632    -23.3813687     -0.2553526      0.0677528     -0.0003020
      119.691084    -23.3818006     -0.2618119      0.0677494     -0.0003096
      121.497391    -23.3831967     -0.2698601      0.0677385     -0.0003190
      123.330957    -23.3834561     -0.2757137      0.0677364     -0.0003260
      125.192194    -23.3835262     -0.2823791      0.0677358     -0.0003338
      127.081520    -23.3833334     -0.2879249      0.0677373     -0.0003404
      128.999359    -23.3827649     -0.2910394      0.0677417     -0.0003441
      130.946140    -23.3833086     -0.2996901      0.0677374     -0.0003543
      132.922301    -23.3835041     -0.3085033      0.0677358     -0.0003647
      134.928286    -23.3834297     -0.3167948      0.0677364     -0.0003745
      136.964543    -23.3845053     -0.3218497      0.0677279     -0.0003805
      139.031530    -23.3842466     -0.3281191      0.0677299     -0.0003879
      141.129711    -23.3836533     -0.3350566      0.0677345     -0.0003961
      143.259556    -23.3839876     -0.3421758      0.0677318     -0.0004045
      145.421544    -23.3837263     -0.3504758      0.0677338     -0.0004143
      147.616160    -23.3834317     -0.3590924      0.0677361     -0.0004245
      149.843895    -23.3832641     -0.3678040      0.0677373     -0.0004348
      152.105249    -23.3836965     -0.3716345      0.0677339     -0.0004393
      154.400731    -23.3839521     -0.3785948      0.0677318     -0.0004476
      156.730855    -23.3841084     -0.3881948      0.0677306     -0.0004589
      159.096144    -23.3845081     -0.3988200      0.0677273     -0.0004714
      161.497128    -23.3846981     -0.4052913      0.0677258     -0.0004791
      163.934346    -23.3844792     -0.4112829      0.0677275     -0.0004862
      166.408346    -23.3834974     -0.4198188      0.0677351     -0.0004963
      168.919681    -23.3837793     -0.4318836      0.0677327     -0.0005106
      171.468916    -23.3836345     -0.4426118      0.0677338     -0.0005233
      174.056623    -23.3828232     -0.4510868      0.0677400     -0.0005333
      176.683382    -23.3830951     -0.4575085      0.0677378     -0.0005409
      179.349782    -23.3836187     -0.4641231      0.0677337     -0.0005487
      182.056422    -23.3840557     -0.4722711      0.0677302     -0.0005583
      184.803909    -23.3834755     -0.4861974      0.0677346     -0.0005748
      187.592860    -23.3831403     -0.4942696      0.0677371     -0.0005844
      190.423899    -23.3830073     -0.5016704      0.0677381     -0.0005931
      193.297663    -23.3831174     -0.5109932      0.0677371     -0.0006041
      196.214797    -23.3828229     -0.5223042      0.0677393     -0.0006175
      199.175953    -23.3829901     -0.5359734      0.0677379     -0.0006337
      202.181798    -23.3835414     -0.5507446      0.0677334     -0.0006511
      205.233005    -23.3830895     -0.5545043      0.0677369     -0.0006556
      208.330259    -23.3820247     -0.5631850      0.0677451     -0.0006659
      211.474256    -23.3809609     -0.5729494      0.0677533     -0.0006775
      214.665699    -23.3808773     -0.5782376      0.0677539     -0.0006838
      217.905306    -23.3831206     -0.5933286      0.0677362     -0.0007015
      221.193803    -23.3845396     -0.6078878      0.0677249     -0.0007186
      224.531928    -23.3839488     -0.6184051      0.0677294     -0.0007310
      227.920430    -23.3849863     -0.6314586      0.0677211     -0.0007464
      231.360069    -23.3851204     -0.6436975      0.0677199     -0.0007608
      234.851617    -23.3846355     -0.6556900      0.0677236     -0.0007751
      238.395858    -23.3849869     -0.6704401      0.0677206     -0.0007925
      241.993586    -23.3849465     -0.6799055      0.0677208     -0.0008037
      245.645609    -23.3846585     -0.6871661      0.0677229     -0.0008123
      249.352746    -23.3841523     -0.6941738      0.0677268     -0.0008206
      253.115829    -23.3840263     -0.7132382      0.0677275     -0.0008431
      256.935703    -23.3839698     -0.7296356      0.0677277     -0.0008625
      260.813223    -23.3839071     -0.7406864      0.0677280     -0.0008756
      264.749261    -23.3830959     -0.7572514      0.0677341     -0.0008953
      268.744699    -23.3827059     -0.7708865      0.0677369     -0.0009114
      272.800434    -23.3824855     -0.7832222      0.0677384     -0.0009260
      276.917375    -23.3817998     -0.7981796      0.0677435     -0.0009438
      281.096447    -23.3807702     -0.8111781      0.0677514     -0.0009593
      285.338587    -23.3802227     -0.8239678      0.0677554     -0.0009745
      289.644747    -23.3808218     -0.8374687      0.0677505     -0.0009904
      294.015893    -23.3804567     -0.8534903      0.0677531     -0.0010093
      298.453006    -23.3804040     -0.8704525      0.0677532     -0.0010294
      302.957081    -23.3808125     -0.8878630      0.0677497     -0.0010499
      307.529128    -23.3819271     -0.9018778      0.0677407     -0.0010664
      312.170175    -23.3813581     -0.9152609      0.0677449     -0.0010823
      316.881261    -23.3802540     -0.9295497      0.0677533     -0.0010993
      321.663444    -23.3803243     -0.9477915      0.0677523     -0.0011209
      326.517797    -23.3812323     -0.9670424      0.0677449     -0.0011435
      331.445409    -23.3820177     -0.9845660      0.0677384     -0.0011641
      336.447386    -23.3823042     -0.9989520      0.0677359     -0.0011811
      341.524849    -23.3813548     -1.0150910      0.0677429     -0.0012003
      346.678938    -23.3808966     -1.0343709      0.0677461     -0.0012232
      351.910810    -23.3809422     -1.0553203      0.0677453     -0.0012479
      357.221639    -23.3814208     -1.0699278      0.0677413     -0.0012651
      362.612615    -23.3827679     -1.0888162      0.0677303     -0.0012873
      368.084948    -23.3838446     -1.1085890      0.0677215     -0.0013105
      373.639867    -23.3832944     -1.1258878      0.0677254     -0.0013310
      379.278617    -23.3832188     -1.1453511      0.0677255     -0.0013540
      385.002464    -23.3832052     -1.1629704      0.0677252     -0.0013749
      390.812692    -23.3831842     -1.1779107      0.0677250     -0.0013925
      396.710604    -23.3836357     -1.2007734      0.0677209     -0.0014195
      402.697524    -23.3833372     -1.2227466      0.0677227     -0.0014455
      408.774795    -23.3827754     -1.2442028      0.0677265     -0.0014709
      414.943780    -23.3833818     -1.2659101      0.0677212     -0.0014965
      421.205865    -23.3828236     -1.2866018      0.0677250     -0.0015211
      427.562452    -23.3823672     -1.3073308      0.0677281     -0.0015456
      434.014970    -23.3828914     -1.3286978      0.0677234     -0.0015708
      440.564865    -23.3821939     -1.3480971      0.0677283     -0.0015939
      447.213608    -23.3822597     -1.3706711      0.0677271     -0.0016205
      453.962689    -23.3832060     -1.3965916      0.0677190     -0.0016510
      460.813623    -23.3834736     -1.4197911      0.0677163     -0.0016784
      467.767948    -23.3835234     -1.4431117      0.0677152     -0.0017059
      474.827223    -23.3835615     -1.4662125      0.0677142     -0.0017332
      481.993032    -23.3839749     -1.4881263      0.0677103     -0.0017590
      489.266984    -23.3834539     -1.5099729      0.0677137     -0.0017849
      496.650710    -23.3832802     -1.5343205      0.0677143     -0.0018138
      504.145866    -23.3840516     -1.5627117      0.0677074     -0.0018471
      511.754135    -23.3837399     -1.5888295      0.0677089     -0.0018781
      519.477224    -23.3836411     -1.6147067      0.0677089     -0.0019087
      527.316864    -23.3838012     -1.6405412      0.0677067     -0.0019392
      535.274816    -23.3842358     -1.6669043      0.0677025     -0.0019702
      543.352865    -23.3839541     -1.6926893      0.0677038     -0.0020008
      551.552822    -23.3837363     -1.7203388      0.0677045     -0.0020335
      559.876529    -23.3845328     -1.7534416      0.0676971     -0.0020724
      568.325851    -23.3845977     -1.7773834      0.0676957     -0.0021007
      576.902686    -23.3844964     -1.8050643      0.0676955     -0.0021334
      585.608958    -23.3843455     -1.8397350      0.0676953     -0.0021744
      594.446619    -23.3845077     -1.8667890      0.0676930     -0.0022063
      603.417653    -23.3849973     -1.8966842      0.0676881     -0.0022415
      612.524073    -23.3854516     -1.9281099      0.0676833     -0.0022785
      621.767921    -23.3845596     -1.9568713      0.0676891     -0.0023127
      631.151272    -23.3847785     -1.9921218      0.0676859     -0.0023543
      640.676231    -23.3848490     -2.0272418      0.0676839     -0.0023958
      650.344935    -23.3837964     -2.0571605      0.0676909     -0.0024314
      660.159553    -23.3844968     -2.0904744      0.0676840     -0.0024706
      670.122288    -23.3848027     -2.1232892      0.0676802     -0.0025093
      680.235374    -23.3844828     -2.1555348      0.0676813     -0.0025475
      690.501081    -23.3844038     -2.1940365      0.0676802     -0.0025930
      700.921711    -23.3848439     -2.2285176      0.0676751     -0.0026336
      711.499604    -23.3852331     -2.2613157      0.0676706     -0.0026722
      722.237132    -23.3843405     -2.2961412      0.0676759     -0.0027136
      733.136704    -23.3842259     -2.3342485      0.0676750     -0.0027586
      744.200767    -23.3846123     -2.3730556      0.0676701     -0.0028043
      755.431801    -23.3855070     -2.4114534      0.0676612     -0.0028494
      766.832327    -23.3847698     -2.4491962      0.0676651     -0.0028942
      778.404904    -23.3853177     -2.4869491      0.0676589     -0.0029386
      790.152127    -23.3868159     -2.5248898      0.0676453     -0.0029829
      802.076632    -23.3859818     -2.5639028      0.0676497     -0.0030292
      814.181094    -23.3856646     -2.6069680      0.0676499     -0.0030802
      826.468230    -23.3857141     -2.6504045      0.0676471     -0.0031315
      838.940796    -23.3861233     -2.6893175      0.0676418     -0.0031773
      851.601590    -23.3863687     -2.7285005      0.0676377     -0.0032234
      864.453454    -23.3860674     -2.7707948      0.0676377     -0.0032735
      877.499270    -23.3850521     -2.8172736      0.0676429     -0.0033287
      890.741966    -23.3864984     -2.8657048      0.0676288     -0.0033853
      904.184513    -23.3866008     -2.9099403      0.0676254     -0.0034375
      917.829927    -23.3861010     -2.9523343      0.0676267     -0.0034878
      931.681269    -23.3877001     -3.0011920      0.0676112     -0.0035448
      945.741648    -23.3833103     -3.0503800      0.0676424     -0.0036046
      960.014217    -23.3775504     -3.1001394      0.0676841     -0.0036658
      974.502179    -23.3735873     -3.1505964      0.0677117     -0.0037271
      989.208786    -23.3741214     -3.2012323      0.0677042     -0.0037867
      1004.137335    -23.3749308     -3.2503083      0.0676946     -0.0038443
      1019.291177    -23.3753900     -3.2978978      0.0676878     -0.0039004
      1034.673712    -23.3749427     -3.3521021      0.0676876     -0.0039646
      1050.288391    -23.3753358     -3.4078712      0.0676807     -0.0040303
      1066.138718    -23.3760129     -3.4636377      0.0676714     -0.0040959
      1082.228248    -23.3758327     -3.5152879      0.0676691     -0.0041569
      1098.560591    -23.3760035     -3.5701969      0.0676638     -0.0042217
      1115.139413    -23.3763133     -3.6244922      0.0676573     -0.0042857
      1131.968432    -23.3766966     -3.6758760      0.0676505     -0.0043462
      1149.051425    -23.3767807     -3.7359762      0.0676452     -0.0044171
      1166.392225    -23.3767601     -3.7972135      0.0676406     -0.0044894
      1183.994721    -23.3766896     -3.8583666      0.0676363     -0.0045616
      1201.862864    -23.3769003     -3.9134206      0.0676303     -0.0046265
      1220.000662    -23.3768294     -3.9793866      0.0676255     -0.0047044
      1238.412185    -23.3766806     -4.0479067      0.0676209     -0.0047853
      1257.101563    -23.3767026     -4.1079384      0.0676157     -0.0048562
      1276.072991    -23.3767501     -4.1771990      0.0676094     -0.0049379
      1295.330724    -23.3766640     -4.2400872      0.0676046     -0.0050121
      1314.879083    -23.3763773     -4.2909558      0.0676024     -0.0050723
      1334.722454    -23.3772792     -4.3582009      0.0675894     -0.0051511
      1354.865290    -23.3776675     -4.4324003      0.0675796     -0.0052384
      1375.312110    -23.3777619     -4.5093112      0.0675718     -0.0053291
      1396.067500    -23.3784778     -4.5711357      0.0675604     -0.0054015
      1417.136119    -23.3783633     -4.6436171      0.0675544     -0.0054871
      1438.522693    -23.3782997     -4.7180895      0.0675478     -0.0055749
      1460.232020    -23.3790835     -4.7873304      0.0675349     -0.0056560
      1482.268971    -23.3794073     -4.8601794      0.0675251     -0.0057417
      1504.638491    -23.3791340     -4.9357043      0.0675196     -0.0058309
      1527.345598    -23.3782841     -5.0138464      0.0675182     -0.0059235
      1550.395388    -23.3788147     -5.0938322      0.0675057     -0.0060174
      1573.793031    -23.3790512     -5.1770382      0.0674951     -0.0061153
      1597.543777    -23.3791601     -5.2610344      0.0674852     -0.0062141
      1621.652956    -23.3794986     -5.3393111      0.0674740     -0.0063061
      1646.125976    -23.3802238     -5.4163183      0.0674599     -0.0063962
      1670.968328    -23.3805428     -5.4978856      0.0674482     -0.0064920
      1696.185586    -23.3798995     -5.5882117      0.0674429     -0.0065988
      1721.783408    -23.3800852     -5.6722877      0.0674317     -0.0066976
      1747.767537    -23.3807559     -5.7613774      0.0674160     -0.0068020
      1774.143803    -23.3817155     -5.8551411      0.0673973     -0.0069115
      1800.918124    -23.3813995     -5.9421998      0.0673892     -0.0070142
      1828.096507    -23.3819923     -6.0340763      0.0673733     -0.0071217
      1855.685050    -23.3827875     -6.1269282      0.0673555     -0.0072302
      1883.689943    -23.3827875     -6.2145451      0.0673443     -0.0073332
      1912.117469    -23.3827405     -6.3130504      0.0673320     -0.0074490
      1940.974006    -23.3826128     -6.4130424      0.0673199     -0.0075667
      1970.266029    -23.3823791     -6.5115145      0.0673086     -0.0076825
      2000.000109    -23.3829045     -6.6122727      0.0672909     -0.0078004
      2030.182919    -23.3835014     -6.7153399      0.0672721     -0.0079209
      2060.821230    -23.3840823     -6.8200895      0.0672530     -0.0080434
      2091.921915    -23.3842712     -6.9234079      0.0672370     -0.0081644
      2123.491954    -23.3834179     -7.0289087      0.0672284     -0.0082890
      2155.538429    -23.3828403     -7.1382273      0.0672169     -0.0084179
      2188.068530    -23.3839364     -7.2539338      0.0671913     -0.0085525
      2221.089557    -23.3849290     -7.3655752      0.0671669     -0.0086824
      2254.608916    -23.3855617     -7.4790001      0.0671446     -0.0088147
      2288.634130    -23.3857489     -7.5951118      0.0671252     -0.0089506
      2323.172833    -23.3853148     -7.7078294      0.0671108     -0.0090831
      2358.232772    -23.3857248     -7.8252516      0.0670889     -0.0092202
      2393.821816    -23.3866072     -7.9464499      0.0670624     -0.0093611
      2429.947948    -23.3871478     -8.0697762      0.0670380     -0.0095048
      2466.619274    -23.3875132     -8.1928706      0.0670146     -0.0096484
      2503.844022    -23.3877203     -8.3197278      0.0669914     -0.0097966
      2541.630544    -23.3877247     -8.4533645      0.0669684     -0.0099528
      2579.987317    -23.3889550     -8.5792151      0.0669369     -0.0100984
      2618.922948    -23.3895263     -8.7092427      0.0669094     -0.0102496
      2658.446172    -23.3893529     -8.8439357      0.0668864     -0.0104071
      2698.565858    -23.3894670     -8.9762203      0.0668613     -0.0105614
      2739.291005    -23.3901545     -9.1131166      0.0668306     -0.0107202
      2780.630752    -23.3911071     -9.2556634      0.0667964     -0.0108853
      2822.594375    -23.3919340     -9.4075476      0.0667610     -0.0110612
      2865.191287    -23.3923660     -9.5462923      0.0667307     -0.0112223
      2908.431046    -23.3927763     -9.6893207      0.0666993     -0.0113883
      2952.323354    -23.3932921     -9.8433206      0.0666645     -0.0115669
      2996.878058    -23.3933915    -10.0005485      0.0666317     -0.0117496
      3042.105156    -23.3936119    -10.1518413      0.0665988     -0.0119252
      3088.014794    -23.3939311    -10.2993584      0.0665654     -0.0120962
      3134.617272    -23.3941726    -10.4543365      0.0665306     -0.0122759
      3181.923048    -23.3952379    -10.6171931      0.0664873     -0.0124634
      3229.942734    -23.3964374    -10.7842631      0.0664415     -0.0126555
      3278.687105    -23.3970977    -10.9524898      0.0663990     -0.0128495
      3328.167097    -23.3966307    -11.1124883      0.0663664     -0.0130356
      3378.393811    -23.3973953    -11.2796655      0.0663223     -0.0132280
      3429.378517    -23.3998971    -11.4575881      0.0662618     -0.0134300
      3481.132653    -23.3994762    -11.6409337      0.0662217     -0.0136427
      3533.667833    -23.4006128    -11.8254797      0.0661687     -0.0138541
      3586.995842    -23.4022014    -12.0101239      0.0661117     -0.0140647
      3641.128645    -23.4009753    -12.1899850      0.0660765     -0.0142742
      3696.078388    -23.4020168    -12.3714460      0.0660230     -0.0144816
      3751.857401    -23.4031945    -12.5578863      0.0659666     -0.0146944
      3808.478196    -23.4031593    -12.7523245      0.0659166     -0.0149182
      3865.953479    -23.4036568    -12.9551804      0.0658596     -0.0151506
      3924.296145    -23.4053205    -13.1595389      0.0657926     -0.0153825
      3983.519284    -23.4079507    -13.3645785      0.0657172     -0.0156131
      4043.636183    -23.4078717    -13.5694535      0.0656616     -0.0158482
      4104.660330    -23.4083605    -13.7774334      0.0655999     -0.0160855
      4166.605417    -23.4093210    -13.9880576      0.0655331     -0.0163247
      4229.485343    -23.4108250    -14.2001494      0.0654609     -0.0165643
      4293.314215    -23.4119516    -14.4199318      0.0653884     -0.0168131
      4358.106354    -23.4133627    -14.6444244      0.0653114     -0.0170664
      4423.876298    -23.4153816    -14.8728001      0.0652277     -0.0173226
      4490.638802    -23.4159754    -15.1064412      0.0651520     -0.0175872
      4558.408847    -23.4173111    -15.3412121      0.0650694     -0.0178513
      4627.201636    -23.4192440    -15.5772950      0.0649809     -0.0181152
      4697.032605    -23.4205041    -15.8159857      0.0648954     -0.0183831
      4767.917422    -23.4219690    -16.0574191      0.0648064     -0.0186533
      4839.871990    -23.4236302    -16.3031360      0.0647135     -0.0189274
      4912.912454    -23.4255622    -16.5552250      0.0646152     -0.0192077
      4987.055200    -23.4276865    -16.8142663      0.0645119     -0.0194948
      5062.316865    -23.4299093    -17.0766886      0.0644055     -0.0197850
      5138.714334    -23.4322135    -17.3415996      0.0642962     -0.0200773
      5216.264748    -23.4343220    -17.6046064      0.0641878     -0.0203673
      5294.985507    -23.4364149    -17.8789156      0.0640741     -0.0206693
      5374.894272    -23.4385681    -18.1596942      0.0639562     -0.0209779
      5456.008973    -23.4410594    -18.4309374      0.0638379     -0.0212743
      5538.347809    -23.4433849    -18.7189495      0.0637131     -0.0215892
      5621.929253    -23.4457244    -19.0136737      0.0635840     -0.0219107
      5706.772059    -23.4481996    -19.3088041      0.0634523     -0.0222316
      5792.895261    -23.4511356    -19.6062923      0.0633146     -0.0225531
      5880.318184    -23.4543072    -19.9070871      0.0631722     -0.0228769
      5969.060442    -23.4576644    -20.2118056      0.0630253     -0.0232035
      6059.141944    -23.4607658    -20.5275662      0.0628740     -0.0235421
      6150.582903    -23.4637791    -20.8403840      0.0627228     -0.0238767
      6243.403835    -23.4667201    -21.1563944      0.0625689     -0.0242141
      6337.625564    -23.4695562    -21.4874118      0.0624076     -0.0245672
      6433.269232    -23.4736595    -21.8221225      0.0622336     -0.0249196
      6530.356297    -23.4780660    -22.1585827      0.0620547     -0.0252718
      6628.908542    -23.4824161    -22.4951449      0.0618742     -0.0256230
      6728.948079    -23.4860996    -22.8466309      0.0616897     -0.0259911
      6830.497353    -23.4901145    -23.2011272      0.0614993     -0.0263601
      6933.579148    -23.4945209    -23.5583588      0.0613026     -0.0267294
      7038.216592    -23.4994376    -23.9197088      0.0610982     -0.0271002
      7144.433162    -23.5045187    -24.2866088      0.0608878     -0.0274748
      7252.252689    -23.5094360    -24.6607846      0.0606727     -0.0278560
      7361.699365    -23.5136514    -25.0452814      0.0604551     -0.0282489
      7472.797744    -23.5185838    -25.4312025      0.0602292     -0.0286392
      7585.572754    -23.5237477    -25.8213968      0.0599971     -0.0290314
      7700.049697    -23.5290427    -26.2164670      0.0597591     -0.0294265
      7816.254257    -23.5352182    -26.6206650      0.0595077     -0.0298261
      7934.212508    -23.5409629    -27.0313095      0.0592532     -0.0302318
      8053.950915    -23.5466697    -27.4485880      0.0589927     -0.0306424
      8175.496342    -23.5538174    -27.8731469      0.0587157     -0.0310531
      8298.876060    -23.5601075    -28.3009995      0.0584398     -0.0314679
      8424.117751    -23.5665765    -28.7343659      0.0581568     -0.0318853
      8551.249516    -23.5740757    -29.1746400      0.0578601     -0.0323033
      8680.299877    -23.5817953    -29.6298547      0.0575504     -0.0327329
      8811.297789    -23.5891623    -30.0877758      0.0572384     -0.0331636
      8944.272644    -23.5961144    -30.5471588      0.0569251     -0.0335946
      9079.254276    -23.6034173    -31.0191887      0.0565988     -0.0340338
      9216.272970    -23.6118669    -31.4930396      0.0562607     -0.0344672
      9355.359469    -23.6209699    -31.9730200      0.0559113     -0.0349007
      9496.544978    -23.6298814    -32.4682941      0.0555505     -0.0353464
      9639.861175    -23.6394180    -32.9679517      0.0551796     -0.0357902
      9785.340214    -23.6490985    -33.4746852      0.0547998     -0.0362364
      9933.014737    -23.6586921    -33.9897237      0.0544117     -0.0366869
      10082.917875    -23.6692294    -34.5147125      0.0540077     -0.0371389
      10235.083262    -23.6802106    -35.0466283      0.0535928     -0.0375911
      10389.545039    -23.6914652    -35.5861385      0.0531675     -0.0380447
      10546.337860    -23.7018987    -36.1391865      0.0527344     -0.0385099
      10705.496906    -23.7137964    -36.6962464      0.0522859     -0.0389673
      10867.057884    -23.7262936    -37.2591578      0.0518259     -0.0394224
      11031.057045    -23.7383971    -37.8291258      0.0513595     -0.0398804
      11197.531184    -23.7515430    -38.4107653      0.0508750     -0.0403386
      11366.517651    -23.7646338    -39.0020285      0.0503800     -0.0407999
      11538.054361    -23.7773214    -39.6025573      0.0498767     -0.0412654
      11712.179800    -23.7925690    -40.2157637      0.0493455     -0.0417235
      11888.933037    -23.8079455    -40.8366692      0.0488040     -0.0421810
      12068.353729    -23.8234664    -41.4655761      0.0482518     -0.0426379
      12250.482131    -23.8392270    -42.1029636      0.0476879     -0.0430938
      12435.359106    -23.8551227    -42.7509798      0.0471112     -0.0435506
      12623.026134    -23.8719132    -43.4103253      0.0465169     -0.0440047
      12813.525321    -23.8905331    -44.0822832      0.0458991     -0.0444518
      13006.899408    -23.9088755    -44.7591596      0.0452751     -0.0448960
      13203.191782    -23.9281089    -45.4445033      0.0446359     -0.0453339
      13402.446483    -23.9484425    -46.1393201      0.0439798     -0.0457646
      13604.708218    -23.9693218    -46.8530941      0.0433020     -0.0461977
      13810.022366    -23.9898773    -47.5741069      0.0426163     -0.0466285
      14018.434993    -24.0109496    -48.3041112      0.0419169     -0.0470533
      14229.992858    -24.0347660    -49.0462791      0.0411908     -0.0474620
      14444.743430    -24.0592223    -49.7986308      0.0404500     -0.0478638
      14662.734888    -24.0841783    -50.5652201      0.0396918     -0.0482619
      14884.016144    -24.1096017    -51.3491604      0.0389137     -0.0486578
      15108.636845    -24.1358585    -52.1319396      0.0381299     -0.0490364
      15336.647388    -24.1632374    -52.9260325      0.0373288     -0.0494042
      15568.098929    -24.1917282    -53.7334170      0.0365090     -0.0497618
      15803.043400    -24.2210009    -54.5620900      0.0356651     -0.0501154
      16041.533512    -24.2517098    -55.3974547      0.0348074     -0.0504514
      16283.622775    -24.2836668    -56.2409069      0.0339359     -0.0507711
      16529.365505    -24.3169939    -57.0921068      0.0330508     -0.0510733
      16778.816838    -24.3505270    -57.9590435      0.0321499     -0.0513689
      17032.032741    -24.3861453    -58.8418632      0.0312263     -0.0516459
      17289.070028    -24.4247275    -59.7412427      0.0302770     -0.0518986
      17549.986369    -24.4624203    -60.6524693      0.0293203     -0.0521468
      17814.840303    -24.5018071    -61.5690971      0.0283535     -0.0523712
      18083.691256    -24.5427977    -62.4950339      0.0273740     -0.0525738
      18356.599546    -24.5850702    -63.4472685      0.0263678     -0.0527641
      18633.626406    -24.6293263    -64.4033308      0.0253543     -0.0529263
      18914.833990    -24.6750045    -65.3709378      0.0243286     -0.0530671
      19200.285391    -24.7219063    -66.3561542      0.0232865     -0.0531896
      19490.044655    -24.7704316    -67.3549250      0.0222312     -0.0532889
      19784.176793    -24.8207773    -68.3650036      0.0211653     -0.0533623
      20082.747798    -24.8730870    -69.3859749      0.0200898     -0.0534083
      20385.824658    -24.9275600    -70.4216246      0.0190016     -0.0534266
      20693.475374    -24.9840913    -71.4709094      0.0179032     -0.0534168
      21005.768971    -25.0425310    -72.5336235      0.0167960     -0.0533793
      21322.775517    -25.1029866    -73.6075223      0.0156831     -0.0533124
      21644.566137    -25.1646033    -74.6926055      0.0145670     -0.0532209
      21971.213029    -25.2282748    -75.7918125      0.0134444     -0.0530999
      22302.789481    -25.2948994    -76.9062308      0.0123143     -0.0529437
      22639.369887    -25.3639905    -78.0314870      0.0111829     -0.0527540
      22981.029764    -25.4354831    -79.1706611      0.0100488     -0.0525318
      23327.845769    -25.5093940    -80.3235866      0.0089137     -0.0522766
      23679.895714    -25.5856413    -81.4839757      0.0077845     -0.0519880
      24037.258587    -25.6645294    -82.6569712      0.0066579     -0.0516651
      24400.014567    -25.7457743    -83.8435040      0.0055346     -0.0513097
      24768.245045    -25.8295111    -85.0427299      0.0044167     -0.0509210
      25142.032637    -25.9167846    -86.2480608      0.0033112     -0.0504927
      25521.461209    -26.0068309    -87.4666318      0.0022136     -0.0500304
      25906.615891    -26.0996881    -88.7001315      0.0011240     -0.0495340
      26297.583097    -26.1946200    -89.9417382      0.0000498     -0.0490082
      26694.450548    -26.2927157    -91.1931488     -0.0010090     -0.0484474
      27097.307286    -26.3940947    -92.4547301     -0.0020514     -0.0478516
      27506.243697    -26.4997145    -93.7242799     -0.0030735     -0.0472168
      27921.351534    -26.6077240    -95.0046553     -0.0040767     -0.0465538
      28342.723930    -26.7187812    -96.2947424     -0.0050587     -0.0458601
      28770.455427    -26.8342275    -97.5889231     -0.0060128     -0.0451303
      29204.641993    -26.9523591    -98.8878356     -0.0069392     -0.0443747
      29645.381044    -27.0729986   -100.1949020     -0.0078400     -0.0435952
      30092.771466    -27.1963049   -101.5099529     -0.0087139     -0.0427920
      30546.913638    -27.3264550   -102.8347789     -0.0095566     -0.0419458
      31007.909453    -27.4584094   -104.1660918     -0.0103699     -0.0410835
      31475.862341    -27.5929244   -105.5050193     -0.0111530     -0.0402026
      31950.877296    -27.7325654   -106.8482378     -0.0118994     -0.0392933
      32433.060894    -27.8740655   -108.1848621     -0.0126057     -0.0383747
      32922.521319    -28.0186655   -109.5276889     -0.0132787     -0.0374403
      33419.368389    -28.1675154   -110.8801604     -0.0139181     -0.0364858
      33923.713578    -28.3215983   -112.2211318     -0.0145085     -0.0355144
      34435.670045    -28.4778894   -113.5747223     -0.0150697     -0.0345348
      34955.352654    -28.6364101   -114.9411673     -0.0156016     -0.0335477
      35482.878003    -28.8009047   -116.2925479     -0.0160810     -0.0325482
      36018.364450    -28.9639837   -117.6773595     -0.0165493     -0.0315520
      36561.932139    -29.1296976   -119.0743570     -0.0169864     -0.0305507
      37113.703028    -29.3056264   -120.4391203     -0.0173541     -0.0295332
      37673.800914    -29.4841504   -121.7985691     -0.0176827     -0.0285209
      38242.351464    -29.6659540   -123.1586330     -0.0179745     -0.0275113
      38819.482240    -29.8514614   -124.5173069     -0.0182283     -0.0265052
      39405.322729    -30.0396977   -125.8607511     -0.0184407     -0.0255116];
  otherwise
    error(['not implemented: ' filt])
end

end