function [res,detailed_spec] = caa_comp_noise_spec(data,whip)

if c_efw_fsample(data)~=450
  error('DATA must be 450 Hz')
end

irf_log('proc','blanking Whisper pulses')
data = caa_rm_blankt(data,whip);

STEP=4*450;
FREQS=50:100;
pos = 0;
count = 0;

if nargout > 1, ncol = 2;
else, ncol = 1;
end
detailed_spec = zeros(fix(size(data,1)/STEP),ncol);

while size(data,1) >=pos+STEP
  % check for nans and data gaps
  if any(isnan(data(pos+1:pos+STEP,2))) || ...
      data(pos+STEP,1) - data(pos+1) > STEP/450.0
    pos = pos + STEP;
    continue
  end
  Pxx = irf_psd(data(pos+1:pos+STEP,:),256,450);
  count = count + 1;
  if nargout > 1
    detailed_spec(count,1) = data(pos+fix(STEP/2),1);
  end
  detailed_spec(count,ncol) = median(Pxx(FREQS));

  pos = pos + STEP;
end
if count < size(detailed_spec,1), detailed_spec(count+1:end,:) = []; end

res = sort(detailed_spec(:,ncol));
res = res(1:fix(length(res)/10));
res = median(res);