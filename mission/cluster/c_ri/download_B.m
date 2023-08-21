function download_B(filename)
%
% download_B(filename)
%
%Input:
% filename	-the path and the filename to with file with the information
%         	about the crossing of the MP. [entering | leaving], where
%        	the times are in epoch.
% block_time - the length of a download block.
%
%Output:
% Saving the B-matrixes to file: /share/robert/B_[from]_to_[to]
%
%Descrition of the function:
% Download the B-matrixes and saves all the
%
%Using:
% isGetContent_lite
% get_sample_fq_for_period
% t_and_dt
% download_B_4_cl
% R_datestring
% saves to: /share/robert/B_data/B_[from]_to_[to]
%
%
%Work method:
%
%Error:
%
%Description of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
% loads the variable: passing_MP
if filename ~= 0
  load(filename)
end

%loads the B-content for satellite 1:4
db = Mat_DbOpen('disco:20');
[B_time1, dur_time1] = isGetContentLite(db, 'Cluster', '1', 'fgm', 'Bprim', ' ', ' ', ' ');
[B_time2, dur_time2] = isGetContentLite(db, 'Cluster', '2', 'fgm', 'Bprim', ' ', ' ', ' ');
[B_time3, dur_time3] = isGetContentLite(db, 'Cluster', '3', 'fgm', 'Bprim', ' ', ' ', ' ');
[B_time4, dur_time4] = isGetContentLite(db, 'Cluster', '4', 'fgm', 'Bprim', ' ', ' ', ' ');
Mat_DbClose(db);

%convering the time array into epochs
%and making four matrixes that have time and duration
%in each column
B_t_dt1 = t_and_dt(B_time1, dur_time1);
B_t_dt2 = t_and_dt(B_time2, dur_time2);
B_t_dt3 = t_and_dt(B_time3, dur_time3);
B_t_dt4 = t_and_dt(B_time4, dur_time4);
save /share/robert/B_data/Btdt B_t_dt1 B_t_dt2 B_t_dt3 B_t_dt4


[nr_crossings , c] = size(passing_MP);

for i =1:nr_crossings
  pack % to reduce gargage in the memory
  [B1, dl1] = download_B_4_cl(passing_MP(i,:),1, B_t_dt1);
  [B2, dl2] = download_B_4_cl(passing_MP(i,:),2, B_t_dt2);
  [B3, dl3] = download_B_4_cl(passing_MP(i,:),3, B_t_dt3);
  [B4, dl4] = download_B_4_cl(passing_MP(i,:),4, B_t_dt4);

  s_fq1 = get_sample_fq_for_period(dl1,B1);
  s_fq2 = get_sample_fq_for_period(dl2,B2);
  s_fq3 = get_sample_fq_for_period(dl3,B3);
  s_fq4 = get_sample_fq_for_period(dl4,B4);

  disp( ' ')
  disp(['size of B-fields: B1 = ' num2str(size(B1)) ' B2 = ' num2str(size(B2)) ' B3 = ' num2str(size(B3)) ' B4 = ' num2str(size(B4)) ])

  date1=R_datestring(fromepoch(passing_MP(i,1)));
  date2=R_datestring(fromepoch(passing_MP(i,2)));
  n_filename = sprintf('/share/robert/B_data/B_%s_to_%s',date1,date2);

  disp( ' ')
  disp(['saving as: ' n_filename])

  save(n_filename, 'B1','B2','B3','B4', 'dl1', 'dl2', 'dl3', 'dl4', 's_fq1' , 's_fq2', 's_fq3', 's_fq4');

  %to avoid memory crasch
  clear B1 B2 B3 B4

end

figure
subplot(2,2,1)
plot(dl1(:,1), dl1(:,2));
