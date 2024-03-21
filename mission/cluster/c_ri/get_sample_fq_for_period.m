function s_fq = get_sample_fq_for_period(download_logg,B)
%
%function s_fq = get_sample_fq_for_period(download_logg,B)
%
%Input:
% download_logg - [start | duration | download successfull]
% B -[start | Bx | By | Bz]
%
%Output:
% s_fq - an array with the sampling fq for each download
%
%Descrition of the function:
% Calculates the sampling frequncy for each download
%
%Using:
% time2row
%
%Work method:
% Finds the first point and the last point for a download segment in B.
% Then the sampling frequncy calulates as:
% (nr of samples in the period) / (the lenght of the period)
%
%Error:
% s_fq is -1 if the intervall is missing.
%
%Description of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
[i_end,c] = size(download_logg);

for i =1:i_end
  if download_logg(i,1) ~=-1

    first_time = download_logg(i,1);
    last_time = first_time + download_logg(i,2)
    f = time2row(first_time,B,1);
    e = time2row(last_time,B,f);

    s_fq(i,1) = first_time;
    s_fq(i,2) = last_time;
    dt = B(e,1) - B(f,1);

    if dt ~= 0
      s_fq(i,3) = (e-f)/dt;
    else
      s_fq(i) = -1;
    end

  else

    s_fq(i) = -1;

  end
end

