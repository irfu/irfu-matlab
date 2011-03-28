function hout=c_pl_sc_orient(spacecraft,time,phase_time_series,magnetic_field,velocity,action)
%C_PL_SC_ORIENT   Plots the orientation of the EFW probes
%
%   h = C_PL_SC_ORIENT;
%   h = C_PL_SC_ORIENT(ic);
%   h = C_PL_SC_ORIENT(ic,t);
%   h = C_PL_SC_ORIENT(ic,t,a);
%   h = C_PL_SC_ORIENT(ic,t,a,b);
%   h = C_PL_SC_ORIENT(ic,t,a,b,v);
%   ic - spacecraft number
%   t  - time in isdat epoch
%   a  - time vector of the satellite phase in degrees
%   b  - magnetic field in GSE reference frame
%   v  - velocity vector [vx vy vz] in GSE which will be marked in the
%        plots, e.g. magnetopause velocity
%
% $Id$

persistent t a b h phase v ic phaseHndl timeHndl figNumber ...
    vec1Hndl vec2Hndl vec1flag vec2flag ...
    flag_v1 flag_v2 v1 v2;

if nargin==1 && ischar(spacecraft)
    action=spacecraft;
    irf_log('proc',['action=' action]);
elseif   (nargin < 6)
    action='initialize';
end
if nargin==0, time=[2005 12 31 01 01 01];spacecraft=1;end
if isempty(ic), ic=spacecraft; end

switch lower(action)
    case 'initialize'
        if length(time)==1, % define time of interest when initializing
            start_time=fromepoch(time);
            t=time;
        elseif length(time)==6,
            start_time=time;
            t=toepoch(time);
        elseif exist('t','var'),
            start_time=fromepoch(t);
        else
            irf_log('fcal','Check time format');return;
        end
        if nargin<6, flag_v=1; end
        if nargin<5, flag_v=0; end
        if nargin>4, % use mangetic field that is given as input
            b=c_coord_trans('GSE','ISR2',magnetic_field,'cl_id',ic);
        end
        if nargin>2, % use phase given as input
            a=phase_time_series;
        end
        if flag_v == 1, v=velocity; end
        % See if spacecraft orientation figures is open
        ch = get(0,'ch');indx=[];
        if ~isempty(ch),
            chTags = get(ch,'Tag');
            indx = find(strcmp(chTags,'cplscor'));
        end
        if isempty(indx),
            figNumber=figure( ...
                'Name',['Cluster s/c' num2str(ic) ' orientation'], ...
                'Tag','cplscor');
            set(figNumber,'Position',[10 10 600 600]);
            delete(findall(gcf,'Type','axes'))
            h=[];
        else
            figure(ch(indx));        
            delete(findall(gcf,'Type','axes'))
            h=[];
            figNumber=gcf;
        end
        set(figNumber,'defaultAxesFontSize',12);
        set(figNumber,'defaultTextFontSize',10);
        h(1)=subplot(2,2,1);axis equal;axis([-50 50 -50 50]);axis manual;title('Spin plane');
        h(2)=subplot(2,2,2);axis equal;axis([-50 50 -50 50]);axis manual;title('View along B');
        h(3)=subplot(2,2,3);axis equal;axis([-50 50 -50 50]);axis manual;title('View towards sun');
        h(4)=subplot(2,2,4);axis off;
        figuserdata=get(figNumber,'userdata');
        figuserdata.h=[];
        figuserdata.h=h;
        set(figNumber,'userdata',figuserdata);

        %====================================
        % The vector 1 entering
        labelStr='0';
        callbackStr='c_pl_sc_orient(''plot'')';
        vec1flag=uicontrol('style','checkbox','units','normalized','Position',[0.5 0.2 .2 .05],'string','vector 1 [GSE]','Callback',callbackStr);
        vec1Hndl=uicontrol( ...
            'Style','edit', ...
            'Units','normalized', ...
            'Position',[0.7 0.2 .2 .05], ...
            'String',labelStr, ...
            'Callback',callbackStr);
        %====================================
        % The vector 2 entering
        labelStr='0';
        callbackStr='c_pl_sc_orient(''plot'')';
        vec2flag=uicontrol('style','checkbox','units','normalized','Position',[0.5 0.25 .2 .05],'string','vector 2 [GSE]','Callback',callbackStr);
        vec2Hndl=uicontrol( ...
            'Style','edit', ...
            'Units','normalized', ...
            'Position',[0.7 0.25 .2 .05], ...
            'String',labelStr, ...
            'Callback',callbackStr);
        %====================================
        % The phase entering
        labelStr='0';
        callbackStr='c_pl_sc_orient(''phase'')';
        uicontrol('style','text','units','normalized','Position',[0.5 0.15 .2 .05],'string','phase [deg]')
        phaseHndl=uicontrol( ...
            'Style','edit', ...
            'Units','normalized', ...
            'Position',[0.7 0.15 .2 .05], ...
            'String',labelStr, ...
            'Callback',callbackStr);
        %====================================
        % The time entering
        labelStr=mat2str(fromepoch(t));
        callbackStr='c_pl_sc_orient(''time'')';
        timeHndl=uicontrol( ...
            'Style','edit', ...
            'Units','normalized', ...
            'Position',[0.5 0.1 .3 .05], ...
            'String',labelStr, ...
            'Callback',callbackStr);
        %====================================
        % The CLOSE button
        labelStr='Close';
        callbackStr='close(gcf)';
        closeHndl=uicontrol( ...
            'Style','pushbutton', ...
            'Units','normalized', ...
            'Position',[0.5 0 .1 .05], ...
            'String',labelStr, ...
            'Callback',callbackStr);
        menus
        c_pl_sc_orient('read_phase_and_b');
    case 'read_phase_and_b'
        if ~isempty(b), % check if getting B data is necessary
            if t>=b(1,1) && t<=b(end,1), % time within interval of B
                flag_get_b_data=0;
            else % get B data
                flag_get_b_data=1;
            end
        else
            flag_get_b_data=1;
        end
        if flag_get_b_data % get B data
            [ok,b]=c_load('diB?',ic);
            if ~any(ok) % no B loaded
                flag_read_isdat=1;
            else
                b(:,3)=-b(:,3);b(:,4)=-b(:,4); % go to DS reference frame instead of DSI
                if t>=b(1,1) && t<=b(end,1), % time within interval of B
                    flag_read_isdat=0;
                else
                    flag_read_isdat=1;
                end
            end
            if flag_read_isdat, % read B in di ref frame from isdat (use CSDD PP data)
                DATABASE=c_ctl(0,'isdat_db');
                data = getData(ClusterDB(DATABASE,c_ctl(0,'data_path')),t-5,5,ic,'b','nosave');
                if isempty(data),
                    irf_log('load','Could not read B field, using B=[0 0 1] nT in DS ref frame'); % first col is time
                    b=[1 0 0 NaN];
                else
                    b=data{3};
                    b=c_coord_trans('GSE','DSC',b,'cl_id',ic);
                end
            end
        end
        if ~isempty(a), % check if getting phasew data is necessary
            if t>=a(1,1) && t<=a(end,1), % time within interval of B
                flag_get_phase_data=0;
            else % get phase data
                flag_get_phase_data=1;
            end
        else
            flag_get_phase_data=1;
        end
        if flag_get_phase_data % get phase data
            [ok,phase_time_series]=c_load('A?',ic);
            if ~any(ok), % check if phase info is on disk
                flag_read_isdat=1;
            else
              a=phase_time_series;
              if t>=a(1,1) && t<=a(end,1), % time within phase interval of B
                flag_read_isdat=0;
              else
                flag_read_isdat=1;
              end
            end
            if flag_read_isdat % read phase form isdat
                [tt,phase_data] = caa_is_get('db.irfu.se:0',t-5,10,ic,'ephemeris','phase_2');
                phase_time_series=[double(tt) double(phase_data)];
                if isempty(phase_time_series),
                    phase_time_series=[1 0]; % default using 0 phase
                end
            end
            a=phase_time_series;
        end
        phase=c_phase(t,a);phase(1)=[];
        set(phaseHndl,'string',num2str(phase,'%3.1f'));
        c_pl_sc_orient('plot');
    case 'time'
        t=toepoch(eval(get(timeHndl, 'string')));
        c_pl_sc_orient('read_phase_and_b');       
    case 'phase'
        phase=str2double(get(phaseHndl, 'string'));
        c_pl_sc_orient('plot');
    case 'plot'
        figuserdata=get(figNumber,'userdata');
        h=figuserdata.h;
        
        flag_v1=get(vec1flag, 'value');
        if flag_v1==1, 
            v1=eval(['[' get(vec1Hndl,'string') ']']);
            if length(v1)==1, flag_v1=0;end;
        end
        flag_v2=get(vec2flag, 'value');
        if flag_v2==1, 
            v2=eval(['[' get(vec2Hndl,'string') ']']);
            if length(v2)==1, flag_v2=0;end;
        end
        
        phase_p1=phase/180*pi + 3*pi/4 ;
        phase_p3=phase_p1     - pi/2   ;
        phase_p2=phase_p1     + pi     ;
        phase_p4=phase_p1     + pi/2 ;
        phase_heea=phase/180*pi-(30)/180*pi;
        phase_leea=phase_heea+pi;
        phase_rapid=phase/180*pi + 60.167/180*pi; % rapid phase
        phase_sunsensor=phase/180*pi + 26.367/180*pi; % the location o fsun sensor
        
        rp1=[44*cos(phase_p1) 44*sin(phase_p1) 0]; % in DS reference frame
        rp2=[44*cos(phase_p2) 44*sin(phase_p2) 0];
        rp3=[44*cos(phase_p3) 44*sin(phase_p3) 0];
        rp4=[44*cos(phase_p4) 44*sin(phase_p4) 0];
        dphi=5/180*pi; % the half size of heea leea rapid azimuthal sectors that are plotted
        sec_length=15; % the length of plotted sectors (the length of efw booms is 44)
        rheea=sec_length*[cos(phase_heea) sin(phase_heea);cos(phase_heea-dphi) sin(phase_heea-dphi);cos(phase_heea+dphi) sin(phase_heea+dphi)];
        rleea=sec_length*[cos(phase_leea) sin(phase_leea);cos(phase_leea-dphi) sin(phase_leea-dphi);cos(phase_leea+dphi) sin(phase_leea+dphi)];
        rrapid=sec_length*[cos(phase_rapid) sin(phase_rapid);cos(phase_rapid-dphi) sin(phase_rapid-dphi);cos(phase_rapid+dphi) sin(phase_rapid+dphi)];
        rsunsensor=sec_length*[cos(phase_sunsensor) sin(phase_sunsensor)];
        
        for ip=1:4,c_eval('rp?_gse=c_coord_trans(''DSC'',''GSE'',[t rp?],''cl_id'',ic);rp?_gse(1)=[];',ip),end
        bfield=irf_resamp(b,t);
        bxs=irf_norm(irf_cross(bfield,[0 0 0 1]));
        bxsxb=irf_norm(irf_cross(bxs,bfield)); % (bxs)xb
        bn=irf_norm(bfield);
        bn_gse=c_coord_trans('DSC','GSE',bn,'cl_id',ic);
        b_elevation=-asin(bn(4))*180/pi;
        angle_deg_p34_vs_b=acos(bn(2)*cos(phase_p4)+bn(3)*sin(phase_p4))*180/pi; % acos(bx*rx+by*ry)
        angle_deg_p12_vs_b=acos(bn(2)*cos(phase_p2)+bn(3)*sin(phase_p2))*180/pi;
        
        if flag_v1==1,
            vn1_gse=[bn(1,1) irf_norm(v1)];
            vn1_ds=c_coord_trans('GSE','DSC',vn1_gse,'cl_id',ic);
            vn1_elevation=-asin(vn1_ds(4))*180/pi;
        end
        if flag_v2==1,
            vn2_gse=[bn(1,1) irf_norm(v2)];
            vn2_ds=c_coord_trans('GSE','DSC',vn2_gse,'cl_id',ic);
            vn2_elevation=-asin(vn2_ds(4))*180/pi;
        end
        
        for ip=1:4,c_eval('rp?_b=[irf_dot(rp?,bxs,1) irf_dot(rp?,bxsxb,1) irf_dot(rp?,bn,1)];',ip),end
        
        aa=0:.1:2*pi;x_circle=cos(aa);y_circle=sin(aa);
        
        axes(h(1));cla
        text(0,50,'sun','verticalalignment','top','horizontalalignment','center','fontweight','demi');
        text(50,0,'dawn','rotation',90,'verticalalignment','bottom','horizontalalignment','center','fontweight','demi');
        patch(x_circle*1.5,y_circle*1.5,x_circle*0+1);hold on; % plot spacecraft
        patch(x_circle*1.5,y_circle*1.5,x_circle*0-1);         % plot spacecraft
        patch([0 rheea(2,2) rheea(3,2)], [0 rheea(2,1) rheea(3,1)],'b'); % heea
        patch([0 rleea(2,2) rleea(3,2)], [0 rleea(2,1) rleea(3,1)],'g'); % leea
        patch([0 rrapid(2,2) rrapid(3,2)], [0 rrapid(2,1) rrapid(3,1)],'k'); % rapid
        line([0 rsunsensor(2)],[0 rsunsensor(1)],[0 0],'Color','y','LineWidth',2); % sun sensor direction
        text(50,50,'Sun sensor','verticalalignment','top','horizontalalignment','right','fontweight','demi','color','y');
        text(50,45,'LEEA','verticalalignment','top','horizontalalignment','right','fontweight','demi','color','g');
        text(50,40,'HEEA','verticalalignment','top','horizontalalignment','right','fontweight','demi','color','b');
        text(50,35,'Rapid','verticalalignment','top','horizontalalignment','right','fontweight','demi','color','k');
        
        bnproj=[0 bn(2)/norm(bn(2:3)) bn(3)/norm(bn(2:3))];
        hl=line([0 bnproj(3)*25],[0 bnproj(2)*25]);set(hl,'color','red','linewidth',.4);       % B direction
        hl=line([0 bn(3)*25],[0 bn(2)*25]);set(hl,'color','red','linewidth',2);       % B direction
        text(30*bnproj(3),30*bnproj(2),'B');           % label
        text(-49,-48,['B elevation=' num2str(b_elevation,2) ' deg'],'fontsize',8);           % label
        if flag_v1==1, % plot v1 vector
            vn1proj=[0 vn1_ds(2)/norm(vn1_ds(2:3)) vn1_ds(3)/norm(vn1_ds(2:3))];
            hl=line([0 vn1proj(3)*25],[0 vn1proj(2)*25]);set(hl,'color','k','linewidth',.4);       % V direction
            hl=line([0 vn1_ds(3)*25],[0 vn1_ds(2)*25]);set(hl,'color','k','linewidth',2);       % V direction
            text(30*vn1proj(3),30*vn1proj(2),'V');           % label
            text(-49,44,['V1 elevation=' num2str(vn1_elevation,2) ' deg'],'fontsize',8);           % label
        end
        if flag_v2==1, % plot v1 vector
            vn2proj=[0 vn2_ds(2)/norm(vn2_ds(2:3)) vn2_ds(3)/norm(vn2_ds(2:3))];
            hl=line([0 vn2proj(3)*25],[0 vn2proj(2)*25]);set(hl,'color','b','linewidth',.4);       % V direction
            hl=line([0 vn2_ds(3)*25],[0 vn2_ds(2)*25]);set(hl,'color','b','linewidth',2);       % V direction
            text(30*vn2proj(3),30*vn2proj(2),'V');           % label
            text(-49,38,['V2 elevation=' num2str(vn2_elevation,2) ' deg'],'fontsize',8);           % label
        end
        
        for aa=0:pi/12:2*pi, % plot grid
            hl=line([0 100*cos(aa)],[0 100*sin(aa)]);
            set(hl,'linestyle',':','color','green','linewidth',.2);
        end
        for ip=1:4; % plot probes
            c_eval('line([0 rp?(2)],[0 rp?(1)]);',ip);
            c_eval('patch(rp?(2)+x_circle*0.4,rp?(1)+y_circle*0.4,x_circle*0+1,''facecolor'',''black'',''edgecolor'',''none'');',ip);
            c_eval('text(rp?(2)*.9,rp?(1)*.9,num2str(?),''fontweight'',''bold'');',ip);
        end
        
        axes(h(3));cla
        text(0,50,'Z_{GSE}','verticalalignment','top','horizontalalignment','center','fontweight','demi');
        text(50,0,'-Y_{GSE}','rotation',90,'verticalalignment','bottom','horizontalalignment','center','fontweight','demi');
        patch(x_circle*1.5,y_circle*1.5,x_circle*0+1);hold on; % plot spacecraft
        patch(x_circle*1.5,y_circle*1.5,x_circle*0-1);         % plot spacecraft
        bnproj=[0 0 bn_gse(3)/norm(bn_gse(3:4)) bn_gse(4)/norm(bn_gse(3:4))];
        hl=line([0 -bnproj(3)*25],[0 bnproj(4)*25]);set(hl,'color','red','linewidth',.4);       % B direction
        hl=line([0 -bn_gse(3)*25],[0 bn_gse(4)*25]);set(hl,'color','red','linewidth',2);       % B direction
        text(-30*bnproj(3),30*bnproj(4),'B');           % label
        if flag_v1==1, % plot v vector
            vn1proj=[0 0 vn1_gse(3)/norm(vn1_gse(3:4)) vn1_gse(4)/norm(vn1_gse(3:4))];
            hl=line([0 -vn1proj(3)*25],[0 vn1proj(4)*25]);set(hl,'color','k','linewidth',.4);       % V direction
            hl=line([0 -vn1_gse(3)*25],[0 vn1_gse(4)*25]);set(hl,'color','k','linewidth',2);       % V direction
            text(-30*vn1proj(3),30*vn1proj(4),'V');           % label
        end
        if flag_v2==1, % plot v vector
            vn2proj=[0 0 vn2_gse(3)/norm(vn2_gse(3:4)) vn2_gse(4)/norm(vn2_gse(3:4))];
            hl=line([0 -vn2proj(3)*25],[0 vn2proj(4)*25]);set(hl,'color','b','linewidth',.4);       % V direction
            hl=line([0 -vn2_gse(3)*25],[0 vn2_gse(4)*25]);set(hl,'color','b','linewidth',2);       % V direction
            text(-30*vn2proj(3),30*vn2proj(4),'V');           % label
        end
        
        for aa=0:pi/12:pi/2,
            hl=line(x_circle*44*sin(aa),y_circle*44*sin(aa));
            set(hl,'linestyle',':','color','green','linewidth',.2);
        end
        for ip=1:4;
            c_eval('line([0 -rp?_gse(2)],[0 rp?_gse(3)]);',ip);
            c_eval('patch(-rp?_gse(2)+x_circle*0.4,rp?_gse(3)+y_circle*0.4,x_circle*0+1,''facecolor'',''black'',''edgecolor'',''none'');',ip);
            c_eval('text(-rp?_gse(2)*.8,rp?_gse(3)*.8,num2str(?));',ip);
        end
        
        axes(h(2));cla
        text(0,50,'(BxS)xS','verticalalignment','top','horizontalalignment','center','fontweight','demi');
        text(50,0,'BxS','rotation',90,'verticalalignment','bottom','horizontalalignment','center','fontweight','demi');
        patch(x_circle*1.5,y_circle*1.5,x_circle*0+1);hold on; % plot spacecraft
        patch(x_circle*1.5,y_circle*1.5,x_circle*0-1);         % plot spacecraft
        for aa=0:pi/12:pi/2,
            hl=line(x_circle*44*sin(aa),y_circle*44*sin(aa));
            set(hl,'linestyle',':','color','green','linewidth',.2);
        end
        for ip=1:4;
            c_eval('line([0 rp?_b(1)],[0 rp?_b(2)]);',ip);
            c_eval('patch(rp?_b(1)+x_circle*0.4,rp?_b(2)+y_circle*0.4,x_circle*0+1,''facecolor'',''black'',''edgecolor'',''none'');',ip);
            c_eval('text(rp?_b(1)*.8,rp?_b(2)*.8,num2str(?));',ip);
        end
        
        % add text 
        cla(h(4))
        irf_legend(h(4),['c_pl_sc_orient() ' datestr(now)],[0,1],'fontsize',8,'interpreter','none','color',[0.5 0.5 0.5]); 
        irf_legend(h(4),['Cluster spacecraft C' num2str(ic)],[0,.9],'fontsize',10);
        irf_legend(h(4),['angle beteen p34 (Ey WBD) and B ' num2str(angle_deg_p34_vs_b,'%3.1f') 'deg'],[0,.8],'interpreter','none','fontsize',10);
        irf_legend(h(4),['angle beteen p12 (Ez WBD) and B ' num2str(angle_deg_p12_vs_b,'%3.1f') 'deg'],[0,.7],'interpreter','none','fontsize',10);
        xp=0;yp=.9;dyp=-0.1;
        yp=yp+dyp;
    case 'c1'
        ic=1;
        c_pl_sc_orient('read_phase_and_b');
    case 'c2'
        if ic~=2
            a=[];b=[];ic=2;
            c_pl_sc_orient('read_phase_and_b');
        end
    case 'c3'
        if ic~=3
            a=[];b=[];ic=3;
            c_pl_sc_orient('read_phase_and_b');
        end
    case 'c4'
        if ic~=4
            a=[];b=[];ic=4;
            c_pl_sc_orient('read_phase_and_b');
        end
end
if nargout,hout=h;end


function menus
% generate menus
if isempty(findobj(gcf,'type','uimenu','label','&Options'))
    hcoordfigmenu=uimenu('Label','&Options');
    uimenu(hcoordfigmenu,'Label','C1','Callback','c_pl_sc_orient(''C1'')')
    uimenu(hcoordfigmenu,'Label','C2','Callback','c_pl_sc_orient(''C2'')')
    uimenu(hcoordfigmenu,'Label','C3','Callback','c_pl_sc_orient(''C3'')')
    uimenu(hcoordfigmenu,'Label','C4','Callback','c_pl_sc_orient(''C4'')')
    user_data = get(gcf,'userdata');
    user_data.coordfigmenu=1;
    set(gcf,'userdata',user_data);
end

