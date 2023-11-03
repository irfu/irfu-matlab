function [angles, ampl1, ampl2] = c_ri_angles(B1,B2)
%
%Input:
%
%Output:
%
%Descrition of the function:
%
%Using:
%
%Work method:
%
%Error:
%
%Description of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------

[i_end ,c ] = size(B1);

for i=1:i_end
  v1 = B1(i,2:4);
  v2 = B2(i,2:4);

  if v1== [0 0 0] || v2 == [0 0 0]
    angles(i) = 0;
    ampl1(i) = norm(v1);
    ampl2(i) = norm(v2);

  else
    n1 = norm(v1);
    n2 = norm(v2);
    ampl1(i) = n1;
    ampl2(i) =n2;
    t1 = v1/n1;
    t2 = v2/n2;

    t3 = dot(t1,t2);
    if t3 > 1
      t3 = 1;
    end

    a = acos(t3)/pi*180;

    if a>180  %to avoid angles greater than 180 degrees
      ang = 360 - a;
    end

    angles(i) = a;
  end
end

[r,c] = size(angles)

if r == 1 && c ~= 1
  angles= angles';
end

ampl1 = ampl1';
ampl2 = ampl2';
