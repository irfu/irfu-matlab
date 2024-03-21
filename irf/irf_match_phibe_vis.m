function out = irf_match_phibe_vis(type,varargin)
% IRF_MATCH_PHIBE_VIS Visualizes the match from other match functions.
%   Used together with irf_match_phibe_dir.m and/or irf_match_phibe_v.m to
%   illustrate the matching that is made.
%
%   output = IRF_MATCH_PHIBE_VIS(type,req_in1,req_in2,...,op_in1,...)
%
%   Examples:
%       % Direction
%       gif_stuff_dir = IRF_MATCH_PHIBE_VIS('direction',x,y,z,corr_dir,intEdt,Bz,Ek,En);
%       imwrite(gif_stuff_dir.im,gif_stuff_dir.map,'mygif_dir.gif','DelayTime',0.01,'LoopCount',inf);
%
%       % Velocity
%       i_n=50; % if more than one densitiy, choose one by specifying index
%       gif_stuff_v = IRF_MATCH_PHIBE_VIS('velocity',phi_E,phi_B(:,[1 i_n]),v,n(i_n));
%       imwrite(gif_stuff_v.im,gif_stuff_v.map,'mygif_v.gif','DelayTime',0.01,'LoopCount',inf);
%
%       % Density vs velocity
%       figure; h=axes;
%       axis_handle = IRF_MATCH_PHIBE_VIS('velocity/density',h,n,v,corr_v);
%
%   See also IRF_MATCH_PHIBE_DIR, IRF_MATCH_PHIBE_V

switch type
  case 'direction' % gif with normalized phi_B and intEdt for different propagation directions
    % Read input
    x=varargin{1};
    y=varargin{2};
    z=varargin{3};
    corr_dir=varargin{4};
    A1=varargin{5}; % intEdt
    A2=varargin{6}; % Bz
    B=varargin{7}; % could be Ek,dEk
    C=varargin{8}; % could be En,dEn

    if ~exist('str_title','var')
      str_title=[];
    end
    % Initialize figure
    fig=figure;
    set(gcf,'color','white'); % white background for figures (default is grey)
    set(gcf,'defaultAxesFontSize',14);
    set(gcf,'defaultTextFontSize',14);
    set(gcf,'defaultAxesFontUnits','pixels');
    set(gcf,'defaultTextFontUnits','pixels');
    set(gcf,'paperpositionmode','auto') % to get the same printing as on screen

    n_frames=max(size(corr_dir));
    n_ims=100; % make hundred images

    set(fig,'color','white');
    set(fig,'position',[560   531   886   395])
    h(1)=axes('position',[0.070    0.640    0.6750    0.270]);
    h(3)=axes('position',[0.070    0.370    0.6750    0.270]);
    h(4)=axes('position',[0.070    0.100    0.6750    0.270]);
    h(2)=axes('position',[0.800    0.370    0.1250    0.2150]); % small direction plot

    tint=[A2(1,1) A2(end,1)];
    index0=find(corr_dir(:,1)==max(corr_dir(:,1))); % mark highest correlation with red indicator
    ind=0;

    for k=fix(linspace(1,n_frames,n_ims))
      ind=ind+1;

      % normalized potential match plot
      irf_plot(h(1),{[A1(:,1) A1(:,k+1)./repmat(max(max(abs(A1(:,k+1)))),size(A1,1),1)],...
        [A2(:,1) A2(:,2)./repmat(max(abs(A2(:,2))),size(A2,1),1)]},'comp');
      set(h(1),'ylim',[-1.1 1.1]);
      ylabel(h(1),'Normalized potential')
      irf_legend(h(1),{'\phi_E','\phi_B'},[0.02 0.9]);

      title(h(1),str_title)
      grid(h(1),'off'); hold(h(1),'off')

      % electric field
      irf_plot(h(3),B(:,[1 k+1])); ylabel(h(3),'E_k'); hold(h(3),'off')
      ylimk=[min(min(B(:,2:end))) max(max(B(:,2:end)))]; set(h(3),'ylim',ylimk);
      irf_plot(h(4),C(:,[1 k+1])); ylabel(h(4),'E_n'); hold(h(4),'off')
      ylimn=[min(min(C(:,2:end))) max(max(C(:,2:end)))]; set(h(4),'ylim',ylimn);
      irf_zoom(h([1 3:4]),'x',tint);
      grid(h(3),'off');grid(h(4),'off')

      % direction plot
      quiver3(h(2),0,0,0,x(k,1),x(k,2),x(k,3));
      hold(h(2),'on')
      quiver3(h(2),0,0,0,x(index0,1),x(index0,2),x(index0,3),'r')
      %quiver3(h(2),0,0,0,n_hat(1),n_hat(2),n_hat(3),'g')
      plot3(h(2),x(:,1),x(:,2),x(:,3),'b')
      axis(h(2),'equal')
      set(h(2),'xlim',1.1*[-1 1],'ylim',1.1*[-1 1],'zlim',1.1*[-1 1])
      hold(h(2),'off')
      view(h(2),irf_cross(x(1,:),y(1,:)))
      corr_str=['corr = ',num2str(corr_dir(k),'%.003f')];
      b_str=['\delta B_{max}=',num2str(max(abs(A2(:,2))),'%.2f'),' nT'];
      vec_str_x=['k=[',num2str(x(k,1),'%0.1f'),' ',...
        num2str(x(k,2),'%0.1f'),' ',...
        num2str(x(k,3),'%0.1f'),']'];
      vec_str_y=['n=[',num2str(y(k,1),'%0.1f'),' ',...
        num2str(y(k,2),'%0.1f'),' ',...
        num2str(y(k,3),'%0.1f'),']'];
      vec_str_z=['B=[',num2str(z(k,1),'%0.1f'),' ',...
        num2str(z(k,2),'%0.1f'),' ',...
        num2str(z(k,3),'%0.1f'),']'];
      title_right_str={b_str,vec_str_x,vec_str_y,vec_str_z,corr_str};
      %title(h(2),['B_{max}=',num2str(max(abs(Bz(:,2))),'%.2f'),'\newline', corr_str]) %'\phi_{max}=',num2str(max(max(abs(phiE(:,2:end)))),'%.2f'),
      title(h(2),title_right_str)
      grid(h(2),'off')

      % collect frames
      f=getframe(fig);
      A(:,ind)=f;
      if k==1 % initialize animated gif matrix
        [im,map] = rgb2ind(f.cdata,256,'nodither');
        im(1,1,1,n_ims) = 0;
      else
        im(:,:,1,ind) = rgb2ind(f.cdata,map,'nodither');
      end
    end
    vis.A=A;
    vis.im=im;
    vis.map=map;
  case 'velocity' % gif with phi_B and phi_E for different v
    % Required input: phi_E,phi_B,v
    phi_E=varargin{1};
    phi_B=varargin{2};
    v=varargin{3};
    if ~isempty(varargin{4}); str_n=num2str(varargin{4},'%.2f'); else, str_n='?'; end
    nv=length(v);

    % Adjust title_str and add density and B0
    if ~exist('title_str','var')
      title_str=[];
    end
    %title_str=[title_str,', ',num2str(B0,'%.f'),' nT, ',num2str(n,'%.2f'), 'cc'];
    ylims=[floor(min(phi_B(:,end))/100) ceil(max(phi_B(:,end))/100)]*110;
    ylims=[floor(min(phi_B(:,end))) ceil(max(phi_B(:,end)))];

    fig=figure('name','Velocity match','position',[560 560 1000 400]);
    set(gcf,'color','white'); % white background for figures (default is grey)
    set(gcf,'defaultAxesFontSize',14);
    set(gcf,'defaultTextFontSize',14);
    ind=0;
    for k=1:nv
      h=irf_plot({phi_B,phi_E(:,[1 1+k])},'comp');
      str_legend{1}=['\phi_B(n=',str_n,'cc)'];
      str_legend{2}=['\phi_E(v=',num2str(v(k),'%0.f'),'km/s)'];
      irf_legend(h,str_legend,[0.02 0.95]);
      irf_zoom(h,'x',[phi_E(1,1) phi_E(end,1)]);
      set(h,'ylim',ylims); grid off;
      ylabel(h,'Potential [V]')
      title(h,title_str)
      hold off;

      % collect frames
      f=getframe(fig);
      A(:,k)=f;
      if k==1 % initialize animated gif matrix
        [im,map] = rgb2ind(f.cdata,256,'nodither');
        im(1,1,1,k) = 0;
      else
        im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
      end
    end
    pause(1)
    vis.A=A;
    vis.im=im;
    vis.map=map;
  case 'velocity/density' % 2D plot of correlation=f(v,n)
    [ax,args,nargs] = axescheck(varargin{:});
    if isempty(ax); ax=axes; end
    n=args{1};
    v=args{2};
    corr_v=args{3};
    pcolor(ax,v,n,log10(abs(corr_v)));
    shading(ax,'flat')
    title(ax,'Correlation')
    xlabel(ax,'Velocity')
    ylabel(ax,'Density')
    axc=colorbar('peer',ax);
    ylabel(axc,'log_{10}|sum(log_{10}|\phi_E.\phi_B|)|')
    vis.ax=ax;
    vis.axc=axc;
end

out=vis;