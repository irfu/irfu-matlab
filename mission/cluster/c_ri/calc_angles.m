function calc_angles(path_input,path_output,filename,sample_fraction)
%
% calc_angles(path_input,path_output,filename,sample_fraction)
%
%Input:
% path_input -
% path_output -
% filename -
% This file should contain B1,B2,B3,B4 which are
% measurements of the magnetic field from the four cluster satellites,
% in GSE-coordinates.
% B = [timestamp| Bx | By | Bz ]
% sample_fraction -which resolution to use
%
%Output:
% the angles and amplitud of each B-vector are saved to file:
% /share/robert/angle/"filename" with "angle" and "ampl"
%
%Descrition of the function:
% Calculates the angles between the B-vectors and the amplitude of the angels.
% The results are saves to file. B1,..,B4 only has to exist, no other specification.
%
%Using:
% vector_angles
% fromepoch
% toepoch
% Ut2number
% time_synch
% four_vector_angles
%
%Work method:
%
%Error:
%
%Description of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
fn =sprintf('%s%s',path_input,filename);
load(fn);

[B1_max,c] = size(B1);
[B2_max,c] = size(B2);
[B3_max,c] = size(B3);
[B4_max,c] = size(B4);

%calculating the sample fq
per1 = (B1(B1_max,1) - B1(1,1))/(B1_max-1);
per2 = (B2(B2_max,1) - B2(1,1))/(B2_max-1);
per3 = (B3(B3_max,1) - B3(1,1))/(B3_max-1);
per4 = (B4(B4_max,1) - B4(1,1))/(B4_max-1);

% The one with the highest sample fq should have lovest packet loss.
sample_period = min([per1 per2 per3 per4]);
samplingrate = 1/sample_period;


beginning = fromepoch(min([B1(1,1) B2(1,1) B3(1,1) B4(1,1)]));
ending = fromepoch(max([B1(B1_max,1) B2(B2_max,1) B3(B3_max,1) B4(B4_max,1)]));


%Converting the starttime to the startpoint in the matrix
%Converting the endtime to the endpoint in the matrix
n_pos(1) = Ut2number(B1(1,1), max(size(B1)), beginning, samplingrate);
n_pos(2) = Ut2number(B2(1,1), max(size(B2)), beginning, samplingrate);
n_pos(3) = Ut2number(B3(1,1), max(size(B3)), beginning, samplingrate);
n_pos(4) = Ut2number(B4(1,1), max(size(B4)), beginning, samplingrate);

disp('The calculations beginns. This may take some time');

%creating a time between the samples and the constructing a time line
%the samplingrate is from cluster one
% dt contains the time between samples with the wanted resolutin, in
% terms of fractions of the samplingrate
% the timeline contains the timestamps to be investigated
start_time = toepoch(beginning);
end_time = toepoch(ending);
dt = sample_period*sample_fraction;
time_line = start_time:dt:end_time;

%goes thrugh the timelinje. For every point in the timeline one vector from
%B1, B2, B3 and B4 is searched for. If the vector a vector within the intervall
% time(i)+-dt/2 is found the that vector is included.
% time the position in the time linje
end_of_i = max(size(time_line));
length_of_B1 = max(size(B1));
tic %starts the clock
for i = 1:end_of_i

  if rem(i,10000) == 0
    est_time = end_of_i/i*toc;
    disp(['the running will take ' num2str(est_time) ' seconds'])
    disp(['calculation so far ' num2str(i/end_of_i*100) ' %'])
  end

  time = time_line(i);
  %writing the time
  angles(i,1) = time;

  [B_vectors, n_pos] = time_synch( time, dt, B1, n_pos(1), B2, n_pos(2), B3, n_pos(3), B4, n_pos(4));
  %the calculation of the angles

  % the amplitude of each vector
  for d = 1:4
    ampl(i,d) =  norm(B_vectors(d,:));
  end

  angles(i,2:7) = four_vector_angles(B_vectors);
end

%gives the total time
disp(['This running took ' num2str(toc) ' seconds'])

a = length(filename);

n_file = sprintf('%sA%s',path_output,filename(2:a))
save(n_file,'angles','ampl');
