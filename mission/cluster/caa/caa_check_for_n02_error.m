function caa_check_for_n02_error(year,month)
%CAA_CHECK_FOR_N02_ERROR: check data to find n*0.2 second errors resulting
%                         from EFW-DWP synchronization failures
%
%  caa_check_for_n02_error(year,month)

old_pwd = pwd;
fighandle=figure('OuterPosition',[0 0 1000 1000]);
irf_log('log_lev',0); % Suppress error messages
pageno=1;
plotno=1;
output_string='# start time             end time                sec  spacecraft\n';
fprintf(1,'Checking %04.0f-%02.0f\nDay ',year,month);
for day=1:eomday(year,month)
  fprintf(1,'%02.0f...',day);
  if mod(day,10)==0,fprintf(1,'\n    ');end
  for cli=1:4
    for t3=0:3:21
      sdir = ['/data/caa/l1/' num2str(year) '/' ...
        num2str(year) num2str(month,'%02.0f') num2str(day,'%02.0f') '_' ...
        num2str(t3,'%02.0f') '00/' ];
      cdir = [sdir 'C' num2str(cli)];

      if ~exist(cdir, 'dir'), continue, end
      d = dir([cdir '/2*_*']);
      if isempty(d), continue, end
      d={d.name};
      d=sort(d);
      for jj=1:length(d)
        curdir = [cdir '/' d{jj}];
        if ~exist([curdir '/.interval'],'file')
          disp(['Error: no .interval file in ' curdir])
          continue
        end
        cd(curdir)

        dphi=caa_comp_sun_angle(cli);
        if length(dphi)>1
          % Blank dphi just after long gaps
          dphi2=dphi;
          idx=find(diff(dphi2(:,1))>60);
          for i=1:length(idx)
            i1=find(dphi2(:,1) > dphi2(idx(i)+1,1)-60, 1 );
            i2=find(dphi2(:,1) > dphi2(idx(i)+1,1)+300, 1 );
            if isempty(i2), i2=length(dphi2(:,1)); end
            dphi2(i1:i2,2)=0;
          end

          % Check for out-of-bounds dphi
          idx= find(abs(dphi2(:,2))>10);
          if length(idx)>10
            if plotno > 4
              plotno=1;
              cd (old_pwd)
              filename=sprintf('comp_sun_angle_%2.2i%2.2i_%3.3i.pdf',year,month,pageno);
              print( fighandle, '-dpdf', filename, '-fillpage' );
              pageno=pageno+1;
              clf
            end
            irf_subplot(4,1,plotno)
            irf_plot(dphi)
            plotno=plotno+1;
            set(gca,'Ylim',[-180 180]);
            txt=[epoch2iso(dphi(idx(1),1)) ' ' ...
              epoch2iso(dphi(idx(end),1)) ' '...
              num2str(mean(dphi(idx,2))*4.1/360) ' '...
              num2str(cli) ];
            title(txt);
            xlabel(' ');
            ylabel('error (deg)')
            output_string=[output_string txt '\n'];
            drawnow;
            pause(0.1);
          end
        end
      end %loop over subintervals
    end % loop over 3-hour intervals
  end %loop over sat
end %loop over days
fprintf(1,'Done.\nPDF files are in the directory %s\n\n',old_pwd);
fprintf(1,output_string);
cd (old_pwd)
filename=sprintf('comp_sun_angle_%2.2i%2.2i_%3.3i.pdf',year,month,pageno);
print( fighandle, '-dpdf', filename, '-fillpage' );
end
