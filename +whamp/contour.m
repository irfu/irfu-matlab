function contour(p,z,f)
% WHAMP.CONTOUR(p,z,f) 
%
% plots contours with some options that are often used
% 
q2 = input('contour or surface, c/s','s');
if q2=='c'
	disp('Levels')
	disp('1 Automaticaly linear')
	disp('2 Automaticaly log')
	disp('3 [0 0.1 0.2 ... 1.0]')
	disp('4 [0 0.01 0.02 0.03 0.04 0.05 0.06]')
	disp('5 [1e-6 1e-5 1e-4 1e-3 1e-2]')
	disp('9 Specify levels')
	q1 = input('');
	if ((q1 == 1) || (q1==2))
		q2=input('Number of levels?');
		if q1==1,cs=whamp.contour(p,z,f,q2);
		else, cs=whamp.contour(p,z,log10(f),q2);
		end
	end
	if q1 >2 
		if q1==3, v = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];end
		if q1==4, v = [0 0.01 0.02 0.03 0.04 0.05 0.06];end
		if q1==5, v=[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1];end
		if q1==9
			q2=input('level1 level2 ...=','s');
			eval(['v=[' q2 ']']);
		end
		cs = whamp.contour(p,z,f,v);
	end
	if (input('Label whamp.contours? y/n','s')=='y'), clabel(cs);end
	xlabel('k_{perp}');ylabel('k_{par}');
else
	pcolor(p,z,f);
	colorbar;
	clear q2
	q2 = input('Color axis (if nothing automatically) cmin cmax =','s');
	if (~isempty(q2))
		caxis(eval(['[' q2 ']']));
	else
		caxis('auto')
	end
	colorbar;
	q2 = irf_ask('shading (flat,faceted,interp) [%] >','q2','flat');
	eval(['shading ' q2]);
end
