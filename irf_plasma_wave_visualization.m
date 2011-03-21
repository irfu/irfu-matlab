function [A,im,map]=irf_plasma_wave_visualization(flag,dB,dVi,dVe,N,f,T,kx,kz,X,Z,init_type)
% Usage:
% [h]=irf_plasma_wave_visualization(flag,dB,dVi,dVe,N,f,T,kx,kz,X,Z,init_type)
%
% visualize plasma waves (ion,electron locations and field lines)
%
% Input:
% 1.    flag  - 'Alfven','whistler', if not recognized use wave polarization
% 2.    dB  - wave polarization (Z is along field, Bo = 1)
% 3.    dVi - ion velocity
% 4.    dVe - electron velocity
% 5.    N  - number of particles
% 6.    f - wave frequency (can be complex), 20 frames per s
% 7.    T - length of simulation
% 8.    kz - wave vector in z, default is 2pi
% 9.    kx - wave vector in x, default is 0
% 10.   X - x length, default is 1
% 11.   Z - Z length, default is 1
% 12.   init_type % 0 - fixed square 1 - extended square 2- zoomed square
%
% Output: 	A - movie matrix
%           im - image matrix
%           map - color map
%
% Examples:
%   irf_plasma_wave_visualization('alfven')
%   irf_plasma_wave_visualization('',[-1i*.05*2*pi .0],[.05 1i*0],[.05 1i*0],1500,1,1,0,2*pi/1,1,1,1)
%   irf_plasma_wave_visualization('',[.1 0],[.01 1i*0],[-.01 1i*0],1500,[.3+1i*.05],8,2*pi/.5,2*pi/1)
%   surface wave
%   irf_plasma_wave_visualization('',[.1 0],[0.03 1i*.03],[0.03 1i*.03],2000,1,5,2*pi/.5,1i*1.5*pi,1,1,1)
%
% to create animated gif
% imwrite(im,map,'anim.gif','DelayTime',0,'LoopCount',inf)

global Xplot_lim Yplot_lim
flag_generate_output_files=0;
if nargout>=1, flag_generate_output_files=1;end
if 1, % check input
    if nargin < 1 ,
        help irf_plasma_wave_visualization;
        return
    end
    
    if nargin < 12, init_type=1; end
    if nargin < 11, Z=1; end
    if nargin < 10, X=1; end
    if nargin < 9, kz=2*pi; end
    if nargin < 8, kx=0; end
    if nargin < 7, if nargin==6, T=2/real(f); else T=2; end; end % default 2 wave periods
    if nargin < 6, f=1; end
    if nargin < 5, N=100; end
    if nargin < 4, dVe=[0,0]; end
    if nargin < 3, dVi=[.1,0]; end
    if nargin < 2, dB=.1; end
    if nargin == 1 && ischar(flag), % example demos
        switch lower(flag)
            case 'alfven'
                Bo=1;                       % background magnetic field is amplitude 1
                dVix=.2;dViz=.2;
                dVex=dVix;dVez=dViz;
                N=1500;
                f=1+1i*.0;w=2*pi/f;T=1;
                Lx=1;kx=2*pi/Lx;
                Lz=1;kz=2*pi/Lz;
                dBx=dVix*kz*Bo/w;
                dBz=-dVix*kx*Bo/w;
                dB=[dBx dBz];
                dVi=[dVix dViz];
                dVe=[dVex dVez];
                init_type=2;
            case 'whistler'
                Bo=1;                       % background magnetic field is amplitude 1
                dVix=0;dViz=0;
                dVex=.2;dVez=.2;
                N=1500;                     % number of particles
                f=1+1i*.0;w=2*pi/f;T=1;
                Lx=1;kx=2*pi/Lx;            % wave period in X
                Lz=1;kz=2*pi/Lz;            % wave period in Z
                dBx=dVex*kz*Bo/w;
                dBz=-dVex*kx*Bo/w;
                dB=[dBx dBz];
                dVi=[dVix dViz];
                dVe=[dVex dVez];
                init_type=2;
        end
    end
end
if 1, % initialize figure
    set(0,'defaultLineLineWidth', 1.5);
    fn=figure(101);clf;
    set(fn,'Renderer','painters')
    set(fn,'color','white')
    set(gcf,'PaperUnits','centimeters')
    xSize = 12; ySize = 12;
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[10 10 xSize*50 ySize*50])
    h=axes('position',[0.1 0.1 0.8 0.8]); % [x y dx dy]
    cla;
    set(h,'visible','off')
end
if 1, % initialize location of particles
    if init_type==0,
        r_i=rand(N,2); % random ions
        r_i(:,1)=r_i(:,1)*X;
        r_i(:,2)=r_i(:,2)*Z;
        r_e=rand(N,2); % random electrons
        r_e(:,1)=r_e(:,1)*X;
        r_e(:,2)=r_e(:,2)*Z;
%        dVxlim=max(abs(dVi(1)),abs(dVi(1)));
%        dVylim=max(abs(dVi(2)),abs(dVi(2)));
        Xplot_lim=[0 X];
        Yplot_lim=[0 Z];
    elseif init_type==1,
        r_i=rand(N,2); % random ions
        r_i(:,1)=r_i(:,1)*X;
        r_i(:,2)=r_i(:,2)*Z;
        r_e=rand(N,2); % random electrons
        r_e(:,1)=r_e(:,1)*X;
        r_e(:,2)=r_e(:,2)*Z;
        dVxlim=max(abs(dVi(1)),abs(dVi(1)));
        dVylim=max(abs(dVi(2)),abs(dVi(2)));
        Xplot_lim=[-dVxlim/f X+dVxlim/f];
        Yplot_lim=[-dVylim/f Z+dVylim/f];
    else
        r_i=rand(N,2); % random ions
        r_i(:,1)=r_i(:,1)*(X+2*abs(dVi(1))/real(f))-abs(dVi(1))/real(f);
        r_i(:,2)=r_i(:,2)*(Z+2*abs(dVi(2))/real(f))-abs(dVi(2))/real(f);
        r_e=rand(N,2); % random electrons
        r_e(:,1)=r_e(:,1)*(X+2*abs(dVe(1))/real(f))-abs(dVe(1))/real(f);
        r_e(:,2)=r_e(:,2)*(Z+2*abs(dVe(2))/real(f))-abs(dVe(2))/real(f);
        Xplot_lim=[0 X];
        Yplot_lim=[0 Z];
    end
    plot_particles(h,r_i,r_e);
end
if 1, % initialize field lines
    Nfield=20; % number of field lines
    Npoints=20; % number of field line points
    field_lines_init=zeros(Nfield,Npoints,2); % x and y coordinates of field lines
    xfield=0:X/Nfield:X;
    for j=1:Nfield,
        field_lines_init(j,:,1)=xfield(j);
        field_lines_init(j,:,2)=0:Z/(Npoints-1):Z;
    end
end
if 1, % step through wave
    dt=1/20; % 20 frames per second
    time_steps=0:dt:T;
    r_ions=zeros([size(r_i) length(time_steps)]);   % initialize matrix
    r_electrons=r_ions;                             % intialize matrix
    for j=1:length(time_steps), % calculate particle positions
        t=time_steps(j);
        dri=1i*dVi/(2*pi*f);
        dre=1i*dVe/(2*pi*f);
        dr_ions(:,1)=real(dri(1).*exp(-1i*2*pi*f*t+1i*kx*r_i(:,1)+1i*kz*r_i(:,2)));
        dr_ions(:,2)=real(dri(2).*exp(-1i*2*pi*f*t+1i*kx*r_i(:,1)+1i*kz*r_i(:,2)));
        dr_electrons(:,1)=real(dre(1).*exp(-1i*2*pi*f*t+1i*kx*r_e(:,1)+1i*kz*r_e(:,2)));
        dr_electrons(:,2)=real(dre(2).*exp(-1i*2*pi*f*t+1i*kx*r_e(:,1)+1i*kz*r_e(:,2)));
        r_ions(:,:,j)=r_i+dr_ions;
        r_electrons(:,:,j)=r_e+dr_electrons;
    end
    if 1, % calculate field line position
        field_lines=zeros(length(time_steps),Nfield,Npoints,2); % field line time evolution
        if kx==0,E=-dB(1)*2*pi*f/kz; else E=dB(2)*2*pi*f/kx;end
        dro=-1i*E/(2*pi*f);
        for j=1:length(time_steps), % calculate field lines position
            t=time_steps(j);
            dr_field=real(dro.*exp(-1i*2*pi*f*t+1i*kx*field_lines_init(:,:,1)+1i*kz*field_lines_init(:,:,2)));
            field_lines(j,:,:,1)=field_lines_init(:,:,1)+dr_field;
            field_lines(j,:,:,2)=field_lines_init(:,:,2);
        end
    end
    if 1, % plot wave animation
        A=moviein(length(time_steps)); % create matlab movie matrix
        for jj=1:length(time_steps), % plot wave animation
            pause(.05);
            plot_particles(h,r_ions(:,:,jj),r_electrons(:,:,jj),squeeze(field_lines(jj,:,:,:)));
            if flag_generate_output_files,
                f=getframe(h);
                A(:,jj)=f;
                if jj==1, % initialize animated gif matrix
                    [im,map] = rgb2ind(f.cdata,256,'nodither');
                    im(1,1,1,20) = 0;
                else
                    im(:,:,1,jj) = rgb2ind(f.cdata,map,'nodither');
                end
            end
        end
    end
end

    function  plot_particles(h,r_ions,r_electrons,field_lines,mark_ions,mark_electrons)
        %        global X Z
        % mark_electrons={marker_size,marker_type,marker_color}
        if nargin<6,
            mark_electrons.markersize=15;
            mark_electrons.marker='.';
            mark_electrons.color='blue';
        end
        if nargin<5,
            mark_ions.markersize=20;
            mark_ions.marker='.';
            mark_ions.color='red';
        end
        if nargin<4, field_lines=[]; end
        if nargin<3, r_electrons=[];end
        if nargin<2, r_ions=[];end
        if nargin<1, h=gca;end
        cla(h);set(h,'xlim',Xplot_lim,'ylim',Yplot_lim);
        if r_ions, % plot ions
            line(r_ions(:,1),r_ions(:,2),'linestyle','none',...
                'markersize',mark_ions.markersize,...
                'marker',mark_ions.marker,...
                'color',mark_ions.color)
        end
        if r_electrons, % plot ions
            line(r_electrons(:,1),r_electrons(:,2),'linestyle','none',...
                'markersize',mark_electrons.markersize,...
                'marker',mark_electrons.marker,...
                'color',mark_electrons.color)
        end
        if ~isempty(field_lines), % plot field lines
            for jjj=1:size(field_lines,1),
                line(field_lines(jjj,:,1),field_lines(jjj,:,2),'linestyle','-',...
                    'color','black')
            end
        end
    end
end
