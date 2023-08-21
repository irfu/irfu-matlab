% whamp.m - used to produce whamp plots
% run whamp using whamp_run.bat and whamp_run.pl in the directory
% during the run will be saved file 'wh' on which works this routine
%
k=-1;
while(k~=0)  % ==== MAIN LOOP ====
  disp('0 quit')
  disp('1 load wh file, 11 - update wh matrice')
  disp('2 new figure')
  disp('3 plot frequency contour plot')
  disp('4 plot other contour plot')
  disp('5 colormaps')
  disp('7 change z,p to lin/log scale')
  disp('8 print')
  disp('9 hold')
  disp('10 matlab session')
  k=input('');
  if k==1 % ========================== k=1 ===
    clear q;
    disp('Loading wh file into variable wh');
    q=input('File name (default ./wh):','s');
    if size(q), eval(['load ' q '; wh=' q ';']);
    else, load .\wh;
    end
    clear q;
    disp('preparing p,z vectors and f,fim matrices')
    [p,z,f,fim] = whamp.m2xyz(wh(:,1:4));
    plin=p;zlin=z;plog=log10(p);zlog=log10(z);
    p=plog;z=zlog;
    zlabel='log_{10} k_{par}';plabel='log_{10} k_{perp}';
  end
  if k==11 % ========================== k=11 ===
    wh_temp = wh;
    load .\wh
    for i1=1:length(wh(:,1))
      for i2=1:length(wh_temp(:,1))
        if ((wh(i1,1)==wh_temp(i2,1)) && (wh(i1,2)==wh_temp(i2,2)))
          wh_temp(i2,:)=wh(i1,:);
        end
      end
    end
    wh=wh_temp;
    clear wh_temp;
    disp('preparing p,z vectors and f,fim matrices')
    [p,z,f,fim] = whamp.m2xyz(wh(:,1:4));
    plin=p;zlin=z;plog=log10(p);zlog=log10(z);p=plog;z=zlog;
  end
  if k==2 % ========================== k=2 ===
    figure
  end
  if k==3 % ========================== k=3 ===
    irf_whamp_contour(p,z,f);
    xlabel(plabel);ylabel(zlabel);
  end
  if k==4 % ========================== k=4 ===
    disp('What to plot?')
    disp('1 Column')
    disp('2 log10(Column)')
    disp('3 Expression, eg. wh(:,3)./wh(:,2)')
    disp('4 (w/z)')
    disp('5 log10(w/z) in eV')
    disp('6 gamma')
    q1 = input('');
    clear xx;
    if q1 == 1
      q2 = input('column number=');
      [d1,d2,d3,xx]=whamp.m2xyz(wh(:,[1 2 3 q2]));
    end
    if q1 == 2
      q2 = input('column number=');
      s = size(wh);
      wh(:,s(2)+1)=log10(wh(:,q2));
      [d1,d2,d3,xx]=whamp.m2xyz(wh(:,[1 2 3 s(2)+1]));
    end
    if q1 == 3
      q2 = input('Expression=','s');
      clear d1;
      d1 = eval(q2);
      s = size(wh);
      wh(:,s(2)+1)=d1;
      [d1,d2,d3,xx]=whamp.m2xyz(wh(:,[1 2 3 s(2)+1]));
    end
    if q1 == 4
      s = size(wh);
      wh(:,s(2)+1)=wh(:,3)./wh(:,2);
      [d1,d2,d3,xx]=whamp.m2xyz(wh(:,[1 2 3 s(2)+1]));
    end
    if q1 == 5
      s = size(wh);
      q3 = 13841^2*5.686e-12/2;
      q4=input('e,H,O (default-e) >','s');
      if strcmpi(q4, 'H'), q3=q3*1836.2; end
      if strcmpi(q4, 'O'), q3=q3*1836.2*4; end
      wh(:,s(2)+1)=log10((wh(:,3)./wh(:,2)).^2*q3);
      [d1,d2,d3,xx]=whamp.m2xyz(wh(:,[1 2 3 s(2)+1]));
    end
    if q1 == 6
      disp('Value to plot xx is calculated as')
      disp('xx = sign(fim).*(MIN+max(-MIN,log10(abs(fim))))');
      disp('gamma=10^(xx-MIN)')
      q2=input('MIN specifies the lowest power for gamma');
      s = size(wh);
      wh(:,s(2)+1) = sign(wh(:,4)).*(q2+max(-q2,log10(abs(wh(:,4)))));
      [d1,d2,d3,xx]=whamp.m2xyz(wh(:,[1 2 3 s(2)+1]));
    end
    clear d1 d2 d3 s;
    irf_whamp_contour(p,z,xx);
    xlabel(plabel);ylabel(zlabel);
  end
  if k==5 % ========================== k=5 ===
    disp('1 light gray colormap (good for printing)')
    disp('2 gamma type colormap')
    q2 = input('');
    if q2 == 1
      gg = flipud(gray(100));
      gr = gg(1:70,:); clear gg
      colormap(gr);
      clear gr;
      colorbar;
    end
    if q2 == 2
      cc = caxis;	% use default caxis
      q3=input('For min/max input 0/1');
      if q3==0, cc(1)=-min(abs(cc));cc(2)=min(abs(cc));end
      if q3==1, cc(1)=-max(abs(cc));cc(2)=max(abs(cc));end
      caxis(cc);

      cm = cc(1):(cc(2)-cc(1))/100:cc(2);
      cm = rot90(cm,-1);
      xcm = ones(length(cm),3);	% colormap matrice
      q1=input('gray scale y/n? ','s');
      q3=input('show positive values y/n','s');
      if q3=='y'
        xcm(:,3)=xcm(:,3)-(sign(cm)+1).*cm/2./cc(2);
        xcm(:,2)=xcm(:,3);
        if q1 == 'y', xcm(:,1)=xcm(:,3);end
      end
      q3=input('show negative values y/n','s');
      if q3=='y'
        xcm(:,2)=xcm(:,2)+(sign(cm)-1).*cm/2./cc(1);
        xcm(:,1)=xcm(:,1)+(sign(cm)-1).*cm/2./cc(1);
        if q1 =='y',xcm(:,3)=xcm(:,3)+(sign(cm)-1).*cm/2./cc(1);end
      end
      if q1=='y',xcm=(xcm+ones(size(xcm)))/2;end % to make it lighter for printer
      colormap(xcm);
      colorbar;
    end
  end
  if k==7 % ========================== k=7 ===
    zp_scale=irf_ask('lin/log scale [%]','zp_scale','lin');
    if strcmp(zp_scale, 'lin'), z=zlin;p=plin;zlabel='k_{par}';plabel='k_{perp}';
    elseif strcmp(zp_scale, 'log'), z=zlog;p=plog;zlabel='log_{10} k_{par}';plabel='log_{10} k_{perp}';
    else, disp(['scale not changed. scale: ' zp_scale]);
    end
  end
  if k==8 % ========================== k=8 ===
    prfig;
  end
  if k==9 % ========================== k=9 ===
    if (eval('ishold') == 0)
      hold on; disp(' HOLD ON ')
    else
      hold off; disp('HOLD OFF')
    end
  end
  if k==10 % ========================== k=10 ===
    q2 = 'disp(''matlab session, q-exit'')';
    while(q2~='q')
      eval(q2);
      q2=input('matlab>','s');
    end
  end
end % ==== MAIN LOOP ====

hold off
disp('bye!')

