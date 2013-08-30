%
% This function implements transformation from ECI_2000/GEI_2000 to GSE
% coordinates. This implementation is heavily based on Mike Hapgood's
% excellent work on coordinate transformations. see Space physics
% coordinate transformations: A user guide, 
% Planet. Space Sci. 40, 711-717 (but see correction note in Hapgood (1997)).
% M.A. Hapgood, (1995) "Space physics coordinate transformations: the role of precession"
% , Annales Geophysicae, 13, 713-716.
% M.A. Hapgood, (1997) "Corrigendum to Space Physics Coordinate Transformations: 
% A User Guide", Planet. Space Sci. 45, 1047. 
%
%
% inpt, is the input matrix or vector.
% time is given in modified julian seconds, since 1950-01-01T00:00:00Z, UT
%
% Ex: 
%      cefTransform([1 0 0], '1950-01-01T00:00:00Z','GSE', 'GEI')
%

% 20061127, Verified and confirmed with fortran routines
%
function out=cefTransform(inpt,time,from, to)


if(not(isnumeric(inpt))) 
 error('Arg 1: not a numeric matrix')
end

if(ischar(time) || iscell(time)) 
    time = cefMjsToMjd2k(cefTimeToMjs(time));
elseif(not(isnumeric(time)))
        error('Arg 2: time must be given either as a CEF formated date string or as a MJS value')
end
    
if(not(ischar(from)))
    disp('You must specify which coordinate system to transform from')
    disp('Available systems are: GEI2000, GEI, GSE')
    error('Arg 3: wrong format')
 end

if(not(ischar(to)))
    disp('You must specify which coordinate system to transform to')
    disp('Available systems are: GEI2000, GEI, GSE')
    error('Arg 4: wrong format')
end


[n,m]=size(inpt);

if(not( (n==3) || (m==3)))
    error('wrong dimension, either row or col must have length 3')
end

if(m==3)
    inpt = inpt';
end

[tr,tc]=size(time);
tlen=length(time);
[ir, ic]=size(inpt);

if(not(tlen==ic))
    error('Arg1 and 2:Vectors must have the same length')
end



ut=[];
for l=1:tlen

%
% Step 0, calculate T0 time and H hour.
%
tid=floor(time(l));
T0 = (floor(time(l))-0.5)/36525.0D0;



%[yy,mm,dd, hr,mn,sec]=dj2000(time(l))
%H = hr+mn/60.0D0+sec/3600.0D0
%H=((cefMjd2kToMjs(time(l))./86400)-floor(cefMjd2kToMjs(time(l))./86400))*24
H=(time(l)-floor(time(l)))*24;

%
% Step 1, calculate variables for the precission transformation
%
%z_a= 0.64062.*T0 + 0.00030.*T0.^2 % degres
%theta_a = 0.55675.*T0 - 0.00012.*T0.^2 % degrees
%zeta_a = 0.64062.*T0 + 0.00008.*T0.^2 % degrees
T=time(l)-0.5D0;
zeta_a= T*(0.3061153D-6 + T*(0.10976D-14 + T*0.179D-20));
z_a = zeta_a + T*T*(0.2881D-14 + T*0.358D-22);
theta_a = T*(0.2660417D-6 - T*(0.1550D-14 + T*0.41549D-20));

%
% Step 2, define the precission matrix
%{
P_1 = [ cos(-z_a)  sin(-z_a)      0
        -sin(-z_a) cos(-z_a)      0  
         0         0              1 ];
     
P_2 = [ cos(-theta_a)        0         sin(-theta_a)   
        0                    1         0
        -sin(-theta_a)       0         cos(-theta_a)  ];

P_3 = [ cos(-zeta_a)  sin(-zeta_a)    0
        -sin(-zeta_a) cos(-zeta_a)    0  
        0             0               1 ];     
%}
P=pr2000(time(l));

%
% Step 3, calculate the variables needed for T2 transformation matrix
%
M = (357.528 + 35999.050*T0 + 0.04107*H)*pi/180; %radians
lambda = 280.460 + 36000.772*T0 + 0.04107*H; %radians
lambda_0 = (lambda + (1.915 - 0.0048*T0)*sin(M) + 0.020*sin(2*M))*pi/180; %degrees

epsilon = (23.439-0.013*T0)*pi/180; %radians
%
% Step 4, define the T2 matrix.
%
T2_1 = [ cos(lambda_0) sin(lambda_0) 0
        -sin(lambda_0) cos(lambda_0) 0
         0             0             1 ];
T2_2 = [ 1   0             0 
         0   cos(epsilon)  sin(epsilon)
         0   -sin(epsilon) cos(epsilon) ];

     
%T2_1
%T2_2
T2 = T2_1*T2_2;

%
% Perform the matrix product and return results
%

switch [lower(from), lower(to)]
    % GEI2000
    case ['gei2000', 'gse']
    %disp('from gei2000 to gse')
%    intmed1 = T2_2*inpt(:,l)
%    intmed2 = T2_1*T2_2*inpt(:,l)
%    intmed3 = T2_1*inpt(:,l)
    out = T2*P*inpt(:,l);       
    case ['gei2000', 'gei']
    %disp('from gei2000 to gei')
    out = P*inpt(:,l);
   
    % GEI
    case ['gei', 'gse']
    %disp('from gei to gse')
    out = T2*inpt(:,l);
    case ['gei', 'gei2000']
    %disp('from gei to gei2000')
    out = inv(P)*inpt(:,l); 
    % GSE
    case ['gse', 'gei']
    %disp('from gse to gei')
    out = inv(T2)*inpt(:,l);
    case ['gse', 'gei2000']
    %disp('from gse to gei2000')
    out = inv(P)*inv(T2)*inpt(:,l);
    otherwise
        error('Unknow coordinate system')
end

ut = [ut; out'];

end

out = ut;