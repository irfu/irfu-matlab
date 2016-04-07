RE = 6371*1e3; % m, Earth radius

% define orbits
r_aps = ([15 25 61]+0); % m, apogee
r_pers = ([3 3 13]+0); % m, perigee
TtotOrb = [1 1 1]*60*60*24*365; % the time to spend in that orbit

dt = 120; % s, timestep, 120s=2min 
t=0; x=[]; y=[];
for orb = [1 2 3]; 
    [tOrb,xOrb,yOrb] = orbit(r_pers(orb),r_aps(orb),round(TtotOrb(orb)),'E','dt',dt);
    % add all orbits in same vectors
    t = [t tOrb+t(end) NaN];
    x = [x xOrb NaN];
    y = [y yOrb NaN];
end
t=t(2:end);
x = x/RE;
y = y/RE;

%% plot the coverage     
% Set up grid    
rmax = ceil(max(abs([x y]))*1.1); % m
binsize = 1; % RE, bin side, 
nbins = 2*rmax/binsize;
edges = linspace(-rmax,rmax,nbins+1);
centers = (edges(1:nbins)+binsize/2);
dt = diff(t(1:2));

% Do the bininng
[X,Y] = meshgrid(edges,edges);
NT = histcn([x(:) y(:)], edges, edges);   
TT = NT*dt;
Ttot = sum(sum(TT));
ndays = sum(sum(TT))/60/60/24;

if 0 % Illustrate gridding    
    hca=subplot(1,3,1);
    hold(hca,'on');    
    mesh(hca,X,Y,X*0-1);
    plot(hca,x,y,'b')
    axis(hca,'equal')
    set(hca,'xlim',rmax*[-1 1],'ylim',rmax*[-1 1])
    title(hca,['Orbit during T_{tot} = ' num2str(ndays,'%.0f') ' days'])
    xlabel(hca,'R_E')
    ylabel(hca,'R_E')
    
    hca=subplot(1,3,2);
    hold(hca,'on');    
    mesh(hca,X,Y,X*0-1);
    plot(hca,x,y,'b',x,y,'bo')
    axis(hca,'equal')
    set(hca,'xlim',[10 12],'ylim',[10 12])
    title(hca,['time between two consecutive ''o'', dt = ' num2str(dt,'%.0f') ' s'])
    xlabel(hca,'R_E')
    ylabel(hca,'R_E')  
    
    hca=subplot(1,3,3);
    ax_pos=get(hca,'position');   
    %[cX,cY] = meshgrid(centers,centers);
    %surf(hca,cX/RE,cY/RE,TT/60/60)
    surf(hca,Y,X,X*0,TT/60)
    view(hca,[0 0 1])
    axis(hca,'equal')
    set(hca,'xlim',rmax*[-1 1],'ylim',rmax*[-1 1])%,'clim',dt*[0 20])
    title(hca,['Orbit coverage, T_{tot} = ' num2str(ndays,'%.0f') ' days'])
    xlabel(hca,'R_E')
    ylabel(hca,'R_E') 
    shading(hca,'flat')
    
    cb = colorbar('peer',hca);
    ylabel(cb,['total time spent in bin [hours]'])
    set(hca,'position',ax_pos);
        
    % grayscale colormap
    cmap=colormap('gray');    
    if 1 % make background white        
        cmap=flipdim(cmap,1);
        cmap(2:3,:)=[]; % make difference between white and first gray step more distinct
    else
        cmap(2:10,:)=[];
    end
    colormap(cmap);
end
if 1 % figure with bowshock + magnetopause            
    hca=axes; hold(hca,'on')              
    surf(hca,Y,X,X*0,TT/60/60)
    view(hca,[0 0 1])
    set(hca,'xlim',rmax*[-1 1],'ylim',rmax*[-1 1])%,'clim',dt*[0 20])
    %title(hca,['Orbit coverage, T_{tot} = ' num2str(ndays,'%.0f') ' days'])
    xlabel(hca,'x_{GSE} [R_E]')
    ylabel(hca,'y_{GSE} [R_E]') 
    shading(hca,'flat')
    
    % add colorbar
    cb = colorbar('peer',hca);
    ylabel(cb,['total time spent in bin [hours]'])            
    
    % grayscale colormap
    cmap=colormap('gray');    
    if 1 % make background white        
        cmap=flipdim(cmap,1);
        cmap(2:3,:)=[]; % make difference between white and first gray step more distinct
    else
        cmap(2:10,:)=[];
    end
    colormap(cmap);
    
    % add bowshock and magnetopause
    % load stand off distances (from Svein)
    filepath = '/Users/Cecilia/M4/thor/orbit_coverage/';
    fid = fopen([filepath 'bowshock_magnetopause_pos_solarwind.txt']);
    format ='%*s%*s%f%f%*f';
    A = textscan(fid,format);
    R_BS = sort(A{1}); % RE
    R_MP = sort(A{2}); % RE
    N = numel(R_BS);
    indR = round([1 0.25*N 0.5*N 0.75*N N]); % percentiles
    %linestyle = {'--','--','-','--','--'};
    linewidth = [1.5 1.5 3 1.5 1.5];
    
    cMP = [0 0.0 1];
    cBS = [1 0.0 0];    
    for kk = 1:numel(indR)%numel(indR):-1:1%1:numel(indR)
        [xmp,ymp] = boundary(R_MP(indR(kk)),'mp');
        plot(hca,-xmp,ymp,'linewidth',linewidth(kk),'color',cMP,'linestyle','--')
        [xbs,ybs] = boundary(R_BS(indR(kk)),'bs');
        plot(hca,-xbs,ybs,'linewidth',linewidth(kk),'color',cBS,'linestyle','-')
    end   
       
    % other stuff
    box(hca,'on')
    set(hca,'ylim',[-20 20],'xlim',[-30 10])
    axis(hca,'square')
    
    % flip x and y labels   
    xticks = get(hca,'xtick'); nticks = numel(xticks);
    xlabels = repmat(' ',nticks,3);
    for pp = 1:numel(xticks)
        xstr = num2str(-xticks(pp)); 
        xlabels(pp,1:numel(xstr)) = xstr;
    end
    yticks = get(hca,'ytick'); nticks = numel(yticks);
    ylabels = repmat(' ',nticks,3);
    for pp = 1:numel(yticks)
        ystr = num2str(-yticks(pp)); 
        ylabels(pp,1:numel(ystr)) = ystr;
    end
    set(hca,'xticklabel',xlabels,'yticklabel',ylabels)
    
    % add Earth
    patch(cos(-pi/2:0.01:pi/2),sin(-pi/2:0.01:pi/2),'k')    
    plot(cos(pi/2:0.01:pi*4/2),sin(pi/2:0.01:pi*4/2),'k')
    
    % set font size
    set(gcf,'defaultAxesFontSize',16);
    set(gcf,'defaultTextFontSize',16);
    set(gcf,'defaultAxesFontUnits','pixels');
    set(gcf,'defaultTextFontUnits','pixels');
end
