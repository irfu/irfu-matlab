function max_angles = find_max_angles(A,AP, min_ampl)
%
%Input:
% A - [time| angle 1 -> angle 6 ]
% AP - [ampl 1 -> ampl 4]
% min_ampl -mininum amplitude
%
%Output:
% max_angles -[time| angle]
%
%Descrition of the function:
% Finds the maximum angle in a row.
%
%Using:
% ind2nr
%
%Work method:
%
%Error:
%
%Description of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------

[i_end,cA] = size( A );
[rAP, cAP] = size (AP);

if cA == 7 && cAP == 4 && i_end == rAP
  for i = 1:i_end

    A_i = A(i,2:7);
    max_angles(i,1) = A(i,1);
    max_angles(i,2) = 0;
    ampl_i = AP(i,:);

    A_max = 0;
    for k = 1:6
      [nr1, nr2] = ind2nr(k);

      if ampl_i(nr1) > min_ampl && ampl_i(nr2) > min_ampl && A_i(k) > A_max
        max_angles(i,2) = A_i(k);
        A_max = A_i(k);
      end

    end

  end

else
  max_angles = 0 ;
end
