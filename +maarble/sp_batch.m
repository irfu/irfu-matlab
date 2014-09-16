dp = pwd;
d = dir;

cd /Users/yuri/plots/
!find . -name C\*_MAARBLE_\*.png | xargs -J % rm -f %
for i=3:length(d), maarble.summary_plot(d(i).name,dp); close all; end

if strcmp(dp(2:8),'Volumes'), dp=['/data' dp(9:end)]; end
[s,w] = unix(['find . -name C\*_MAARBLE_\*.png | xargs -J % scp % yuri@db.irfu.se:' dp(1:end-3) 'PNG']);
disp(w)
