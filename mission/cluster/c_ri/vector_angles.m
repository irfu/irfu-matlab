function ang = vector_angles(v1 , v2)
%
%Input:
% Two vectors of any size
%
%Output:
% The angle between the vectors, in degres and in absolute value. All angles are
% in the range 0-180 degrees.
%
%Descrition of the function:
% This function get the angles between to vectors
% the return value is in absolute angle, so it doesnt matter
% which point is which in the inparameter to the function. The return value
% is in degrees (0-360).
%
%Work method:
% The first vector, v1 ,is projected onto the second vector, v2, and the scalar
% a is how much the have in common.(the dot product). Then the a*v2 is subtracted
% from the vector v1 and the this is the vector v3. Vector v3 is perpendicular
% to vector v2. The scalar b is the length of vector v3. The angle between the
% vectors can be calculated using arctan of the to scalars, who are representing
% the length in to perpendicular directions.
%
%Error:
% If a  zero-vector is input the result will be zero angle.
%
%Description of variables:
% a and b are the constants infront of the unit length vectors
% v1, v2 are the input
% v3 is vector perpendicular to v2
% angles is the absolute vaule of the angles in radians
% ang is is the return value in degrees
% ang_pi is the value in fractions of pi
%
%Written by Robert Isaksson in the summer of -03
%--------------------- the beginning --------------------------
% error check
if norm(v1) == 0 || norm(v2) == 0
  %'Error, zero-vector as input, returning 0-angle'
  ang = 0;

  %the angle calculation
else
  % normalising, since the length of the vector is not of importance
  v1 = v1 / norm(v1);
  v2 = v2 / norm(v2);


  a = dot(v1, v2);
  angles = abs(acos(a));
  ang_pi = angles/pi;
  ang = (angles)/pi*180;

  if ang>180  %to avoid angles greater than 180 degrees
    ang = 360 - ang;
  end

end
