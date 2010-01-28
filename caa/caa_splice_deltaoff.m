function deltaOff=caa_splice_deltaoff(deltaOff_in,iso_t)
%CAA_SPLICE_DELTAOFF
% Use to splice together discontinuous delta offsets in such a way
% that the spline interpolation doesn't go crazy.

deltaOff=deltaOff_in;

for i=1:8
    idx=find(deltaOff(:,1)>iso2epoch(iso_t), 1 );
    n_pts=4;
    t1=deltaOff(idx-2,1);
    t2=deltaOff(idx+1,1);
    dt=(t2-t1)/(n_pts+1);
    t1=t1+dt;
    t2=t2-dt;
    splice=[(t1:dt:t2)' ...
        [deltaOff(idx-2,2)*ones(1,n_pts/2) deltaOff(idx+1,2)*ones(1,n_pts/2)]' ...
        [deltaOff(idx-2,3)*ones(1,n_pts/2) deltaOff(idx+1,3)*ones(1,n_pts/2)]'];
    deltaOff=[deltaOff(1:idx-2,:)' splice' deltaOff(idx+1:end,:)']';
end

%plot
t=iso2epoch(iso_t)-60*86400:3600:iso2epoch(iso_t)+60*86400;
deltaint=[t' interp1(deltaOff(:,1),deltaOff(:,2:3),t,'spline')];
idx=find(deltaOff_in(:,1)> min(t) & deltaOff_in(:,1) < max(t));
idx2=find(deltaOff(:,1)> min(t) & deltaOff(:,1) < max(t));
clf
h=irf_plot({deltaOff_in(idx,:),deltaOff(idx2,:),deltaint},'linestyle',{'-','*','-'});
ylabel(h(1),'Old delta offsets')
ylabel(h(2),'New offsets')
ylabel(h(3),'Interpolation on new')
