function [B1,B2,B3,B4]=c_get_bfgm(tint,sc_list)
%function [B1,B2,B3,B4]=c_get_bfgm(tint)
%
% tint given in isdat_epoch
% sc_list .. vector of satellite unmbers for which to dowload data, default sc_list=1:4
if nargin ==1, sc_list=1:4;end
for ic=sc_list,eval(av_ssub('B?=[];',ic)),end

for ic=sc_list,
  [b_table, n_table] = create_timetable(tint(1),tint(2),ic);

  if b_table ~= 0
    [b_s, col] = size(b_table);
    for k = 1:b_s
      B=c_ri_get_B( b_table(k,1), b_table(k,2), ic, 'b');
      if B(1,1)>0,
        eval(av_ssub('B?=[B?;B];',ic));
      end
    end
  end

  if n_table ~= 0
    [n_s, col] = size(n_table);
    for i = 1:2:n_s
      B=c_ri_get_B( n_table(i,1), n_table(i,2), cl_nr, 'n');
      if B(1,1)>0,
        eval(av_ssub('B?=[B?;B];',ic));
      end
    end
  end
end
