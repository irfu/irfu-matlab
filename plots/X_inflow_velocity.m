% plot inflow velocity as a funciton of distance to current sheet
%
% model 2D reconnection. y- out of plane direction, x-normal direction,
% z-along the outflow
%
% Bx= constant, Bz - changes from 0 in the center of current sheet to 1 on
% edge
% Ey=1 - tangential electric field is constant
%
% v_inflow=Ey*Bz/B^2=Ey*Bz/(Bz^2+Bx^2)

Bx=[0.05 0.1 0.2 0.5 1.]; % different cases
Bz=0:0.01:1; Bz=Bz(:);
clear vin vout;
legend_str='legend(';
for j=1:length(Bx)
  vin(:,j)=Bz./(Bz.*Bz+Bx(j)*Bx(j));
  vout(:,j)=Bx(j)./(Bz.*Bz+Bx(j)*Bx(j));
  legend_str=[legend_str '''Bx=' num2str(Bx(j),2) ''','];
end
legend_str(end)=[]; % remove last comma
legend_str=[legend_str ');']; %

%%%%%%%%%%%%%%  Figure 1 %%%%%%%%%%%%%%

figure;
subplot(2,1,1);
plot(Bz,vin)
eval(legend_str);
title('2D reconnection geometry. x-normal,z-outflow.')
grid on
xlabel('Bz')
ylabel('V inflow')
ht=irf_pl_info([mfilename '  ' char(datetime("now","Format","dd-MMM-uuuu HH:mm:ss"))],gca,[0,1 ]);

subplot(2,1,2);
plot(Bz,vout)
eval(legend_str);
grid on
xlabel('Bz')
ylabel('V outflow')


