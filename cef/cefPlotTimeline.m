% 
% Function for parsing directories and estimate the timeline from
% files in that directory.
% Ex:
%
%  dir1 =../Interpolate_CEFfiles/efw/
%  dir2 = ../Interpolate_CEFfiles/fgm/
%  cefPlotTimeline(dir1,dir2)
%
%  Current version only handles efw and fgm files.
%
%
function cefPlotTimeline(varargin)


% xaxis length
xlen =1000;
date=[];
startdate=[];
enddate=[];
totfiles={};
colormodes={};
maxy=0;


%loop over all the files
for p=1:length(varargin)

    if(not(ischar(cell2mat(varargin(p)))))
        error('argument %d is not a string, %s',p, cell2mat(varargin(p)))
    end

    files=dir(cell2mat(varargin(p)));
    % start at 3, skipping '.' and '..' 
    for i=3:length(files)
    
        
        n=files(i).name;
        totfiles=[totfiles, {n}];
    
        [foo,n]=cefSplit(n,'__','g');
        splt=cefSplit(n,'_','g');
        [len, foo]=size(splt);

        % First file format type
        if(len==4)
           
            test= splt(len-2,:);

            if(str2num(test(1:2))==0)
                test(1:2)='24';
                splt(len-2,:)= test;
            end

            test= splt(len-1,:);

            if(str2num(test(1:2))==0)
                test(1:2)='24';
                splt(len-1,:)= test;
            end

            date=strvcat(date,[splt(len-3,:),splt(len-2,:),splt(len-1,:)]);
            startdate=[startdate; datenum(strcat(splt(len-3,:),splt(len-2,:)),'yyyymmddHHMMSS')];
            enddate=[enddate; datenum(strcat(splt(len-3,:),splt(len-1,:)),'yyyymmddHHMMSS')];
            colormodes=[colormodes, {'g'}];

            % Second file format type
        elseif(len==5)

            date=strvcat(date,[splt(len-4,:),splt(len-3,:),splt(len-2,:),splt(len-1,:)]);
            startdate=[startdate; datenum(strcat(splt(len-4,:),splt(len-3,:)),'yyyymmddHHMMSS')];
            enddate=[enddate; datenum(strcat(splt(len-2,:),splt(len-1,:)),'yyyymmddHHMMSS')];

            colormodes=[colormodes, {'m'}];


        end
    end
end

mi=min(startdate);
ma=max(enddate);

plot(linspace(mi,ma,xlen), zeros(1, xlen), '-s','MarkerEdgeColor','r',...
    'MarkerFaceColor','r'), hold on

for m=1:length(startdate)
    y=interp1(linspace(startdate(m), enddate(m),xlen), (m+1)*ones(1,xlen), linspace(mi,ma,xlen));
    if(maxy<max(y))
        maxy=max(y);
    end
    plot(linspace(mi,ma,xlen),y,'-s','MarkerEdgeColor', cell2mat(colormodes(m)),'MarkerFaceColor',cell2mat(colormodes(m)));
    h=text(mi,m+1+0.3,cell2mat(totfiles(m)));
    set(h,'Interpreter','none');
    set(h,'FontSize',8);
end


axis([mi ma 0 round(maxy)+1])
datetick('x',21,'keepticks')
grid minor
