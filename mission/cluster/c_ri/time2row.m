function row = time2row(time,M,start_r)
%
% row = time2row(time,M,start_r)
%
%Input:
% time - in epoch
% M -time in epoch in the left column
% start_r -the row to start look (Optional)
%
%Output:
% row -the row where the time i found.
%
%Descrition of the function:
% Matches the input time to a row
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
if nargin ==2
  start_r = 1;
end

[r,c] = size(M);

if start_r > r
  start_r = r;
end

for i = start_r:r

  if M(i,1) >= time
    row = i;
    break
  end

  if i == r
    row = i,
  end

end
