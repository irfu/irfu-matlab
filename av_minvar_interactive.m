function status = av_minvar_interacitve(x,column)
% AV_MINVAR_INTERACTIVE interactively do the minimum variance analysis
% AV_MINVAR_INTERACTIVE(X,COLUMN)
%  X - vector to use, [x(:,column(1)) x(:,column(2)) x(:,column(3))]
%  COLUMN - which columns to use, if not given use 2,3,4
%
% You can access the results through variable 'ud' that is defined as global
% ud.l - eigenvalues  ud.l(1), ud.l(2),ud.l(3)
% ud.v - eigenvectors (ud.v(1,:), ..), also ud.v1, ud.v2. ud.v3
%  ud.Xminvar - data in minimum variance coordinates
%

evalin('caller','clear ud; global ud;');
if nargin < 1, help av_minvar_interactive;return; end
if size(x,2)<3, disp('Vector has too few components');return;end
if nargin < 2,
 if size(x,2)==3, column=[1 2 3];end
 if size(x,2)>3, column=[2 3 4];end
end

% X is used for minimum variance estimates

if min(column)==1, time_vector=1:size(x,1);
elseif min(column)>1, time_vector=x(:,1);
end

X=[time_vector x(:,column)];X=av_abs(X);
dgud={}; % structure to pass all information to manager function
dgud.X=X;
dgud.from = 1; % first click with mouse is 'from', second is 'to'
dgud.cancel = 0;
tlim = [min(X(:,1)) max(X(:,1))];

dgh=figure;clf;av_figmenu;
h(1)=subplot(4,1,1);
av_tplot(X);
set(h(1),    'buttondownfcn', 'av_minvar_interactive_manager(''ax'')');zoom off;
av_pl_info(['av\_minvar\_interactive() ' datestr(now)]); % add information to the plot
set(h(1),'layer','top');
ax=axis;grid on;
legend('x','y','z','abs');
dgud.mvar_intervals=patch([tlim(1) tlim(2) tlim(2) tlim(1)],[ax(3) ax(3) ax(4) ax(4)],[-1 -1 -1 -1],'y','buttondownfcn', 'av_minvar_interactive_manager(''ax'')');

h(2)=subplot(4,1,2);

h(3)=subplot(4,2,5);

h(4)=subplot(4,2,6);

dgud.h=h;

xp=0.2;yp=0.2;
dgud.fromtext=uicontrol('style', 'text', 'string', 'From:', 'position', [xp yp 0.1 0.03],'units','normalized','backgroundcolor','red');
dgud.fromh = uicontrol('style', 'edit', ...
      'string', strrep(datestr(datenum(fromepoch(tlim(1))), 0),' ','_'), ...
    'callback', 'av_minvar_interactive_manager(''from'')', ...
    'backgroundcolor','white','position', [xp+0.11 yp 0.2 0.03],'units','normalized');

yp=0.15;
dgud.totext=uicontrol('style', 'text', 'string', 'To:', 'position', [xp yp 0.1 0.03],'units','normalized','backgroundcolor','white');
dgud.toh=uicontrol('style', 'edit', ...
    'string', strrep(datestr(datenum(fromepoch(tlim(2))), 0),' ','_'), ...
    'callback', 'av_minvar_interactive_manager(''to'')','backgroundcolor','white', 'position', [xp+0.11 yp 0.2 0.03],'units','normalized');


xp=0.15;yp=0.1;
uch1 = uicontrol('style', 'text', 'string', 'Low pass filter f/Fs = ','position', [xp yp 0.15 0.03],'units','normalized','backgroundcolor','white');
dgud.filter = uicontrol('style', 'edit', ...
      'string', '1', ...
    'callback', 'c_4_v_update(''dt'')', ...
    'backgroundcolor','white','position', [xp+0.15 yp 0.1 0.03],'units','normalized');

subplot(4,2,8);axis off;
dgud.result_text=text(0,0.8,'result');

set(dgh, 'userdata', dgud);

gcbf=gca;
av_minvar_interactive_manager('from');

return


