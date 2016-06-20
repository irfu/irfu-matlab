% A script with defined orbit parameters, that gives the orbit in variables
% t - time in seconds, x - 'x_GSE' in m, y - 'y_GSE' in m

RE = 6371*1e3; % m, Earth radius

% define orbits
r_aps = ([15 26 61]+0); % m, apogee
r_pers = ([4 4 14]+0); % m, perigee
TtotOrb = [1 1 1]*60*60*24*365; % the time to spend in that orbit

dt = 120; % s, timestep, 120s=2min 
t=[]; x=[]; y=[]; tend = 0;
npl = 1;
for nOrb = [1 2 3]; 
    [tOrb,xOrb,yOrb] = orbit(r_pers(nOrb),r_aps(nOrb),round(TtotOrb(nOrb)),'E','dt',dt);
    % add all orbits in same vectors
    t = [t tOrb+tend NaN];
    x = [x xOrb NaN];
    y = [y yOrb NaN];
    tend = tend + tOrb(end);
    if 0
        hca=subplot(3,2,npl); npl = npl+1;
        plot(hca,tOrb/60/60/24,[xOrb' yOrb']/RE)
        hca.XLabel.String = 'Time [days]';
        hca=subplot(3,2,npl); npl = npl+1;
        plot(hca,t/60/60/24,[x' y']/RE)
        hca.XLabel.String = 'Time [days]';        
    end
end
%t=t(2:end);
x = x/RE;
y = y/RE;