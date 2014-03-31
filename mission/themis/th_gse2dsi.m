load mRth.mat
gseR = Rtha;

x = irf_gse2gei(gseR);
inp = x(:,[2 3 4]); % assuming first column is time
t = x(:,1);

d=dataobj('tha_l1_state_20070804.cdf');
spinras=getmat(d,'tha_spinras');
spindec=getmat(d,'tha_spindec');

latlong = [spindec spinras(:,2)];
[Rx,Ry,Rz] = sph2cart(latlong(3)*pi/180,latlong(2)*pi/180,1);
a = 1/sqrt(Ry^2+Rz^2);
M = [[a*(Ry^2+Rz^2) -a*Rx*Ry -a*Rx*Rz];[0 a*Rz	-a*Ry];[Rx	Ry	Rz]];
        
out = M*inp';
out = [t out'];
