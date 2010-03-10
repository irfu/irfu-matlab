function res = caa_comp_noise_spec(data,whip)

if c_efw_fsample(data)~=450,
	error('DATA must be 450 Hz')
end

irf_log('proc','blanking Whisper pulses')
data = caa_rm_blankt(data,whip);

STEP=450;
FREQS=50:100;
pos = 0;
count = 0;

while length(data) >=pos+STEP
	% check for nans and data gaps
	if any(isnan(data(pos+1:pos+STEP,2))) || ...
		data(pos+STEP,1) - data(pos+1) > STEP/450.0
		pos = pos + STEP;
		continue
	end
	Pxx = irf_psd(data(pos+1:pos+STEP,:));
	count = count + 1;
	res(count,1) = data(pos+fix(STEP/2),1);
	res(count,2) = median(Pxx(FREQS));
	pos = pos + STEP;
end