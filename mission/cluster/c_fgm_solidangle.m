function angle=c_fgm_solidangle(var1,var2,var3)
%purpose: calculate the solid angle of three vectors in a unit sphere
%
%form: angle=solidangle(var1,var2,var3)
%
%input: var1 to var3 are three vectors 
%output: angle is the solid angle of three input vectors
%
%  See also c_fgm_poincare_index
%
%--------written by Y.H.Hu--------------------------------------------
trilength1=dot(var1-var2,var1-var2,2);
trilength2=dot(var2-var3,var2-var3,2);
trilength3=dot(var3-var1,var3-var1,2);
%
alpha1=acos((sqrt(trilength1))/2.0);
beta1=acos((sqrt(trilength3))/2.0);
theta1=acos((trilength1+trilength3-trilength2)...
      ./(2*sqrt(trilength1.*trilength3)));
phi1=acos((cos(theta1)-cos(alpha1).*cos(beta1))./(sin(alpha1).*sin(beta1)));
%    
beta2=acos((sqrt(trilength2))/2.0);
theta2=acos((trilength1+trilength2-trilength3)...
          ./(2*sqrt(trilength1.*trilength2)));
phi2=acos((cos(theta2)-cos(alpha1).*cos(beta2))./(sin(alpha1).*sin(beta2)));
%
alpha2=acos((sqrt(trilength3))/2.0);
theta3=pi-theta1-theta2;
phi3=acos((cos(theta3)-cos(alpha2).*cos(beta2))./(sin(alpha2).*sin(beta2)));
%
angle=(phi1+phi2+phi3-pi);
%
 var=cross(var2-var1,var3-var1,2);
 div=dot(var,-var1,2);
 temp=div./abs(div);   
angle=(angle.*temp)/(4*pi);