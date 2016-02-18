%% Test data
T   = EpochTT('2002-03-04T09:30:00Z'):.2...
     :EpochTT('2002-03-04T10:30:00Z');      % define time line as EpochTT object
t   = T.tts - T.tts(1);                     % define relative time in s from start
x   = exp(0.001*(t)).*sin(2*pi*t/180);      % define function x(t)=exp(0.001(t-to))*sin(t-to)
TS1 = irf.ts_scalar(T,x);                   % define scalar TSeries object
y   = exp(0.001*(t)).*cos(2*pi*t/180);	      % y(t)=exp(0.001(t-to))*cos(t-to)
TS2 = irf.ts_vec_xy(T,[x y]);
TS3 = irf.ts_vec_xy(T,[.8*x 1.1*y]);

%%
h = irf_plot(2);
hca = irf_panel('panel A');
irf_plot(hca,TS1)
hca = irf_panel('panel B');
irf_plot(hca,TS2)

%%
%XXX TODO: must be using array H as input
%h = irf_plot(2);
%irf_plot(h,{TS1,TS2})
irf_plot({TS1,TS2})

%%
irf_plot({TS2,TS3},'comp')

%%
h = irf_plot(2);
hca = irf_panel('panel A');
irf_plot(hca,{TS2.x,TS3.x},'comp')
hca = irf_panel('panel B');
irf_plot(hca,{TS2.y,TS3.y},'comp')
