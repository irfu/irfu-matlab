function dist = dist_to_MP_shue(swp, bz, xgsm, ygsm, zgsm) 
%
%Input:
%
%Output:.
%
%Descrition of the function:
%
%Using:
% 
%Work method:
%
%Error:
% 
%Discription of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
 %[bz]=nT, [swp]=nPa, [teta]=radians
      r = norm([xgsm ygsm zgsm]);
	  cosTheta = xgsm/r;
	  
      cc=-1.0/6.6;
      alfa=(0.58-0.007*bz)*(1.0+0.024*log(swp));
      r0=(10.22+1.29*tanh(0.184*(bz+8.14)))*swp^cc;
      rmpause = r0 *(2.0/(1.0+cosTheta))^alfa;
	  dist = rmpause - r;

