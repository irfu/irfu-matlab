%% THOR key science region definition
%% Bow shock
thetaRef = 0:5:30;
disp(thetaRef);
out={};
out{1}='======';
out{2}='theta ';
out{3}='======';
for ii = 1:numel(thetaRef)
	out{ii+3} = [num2str(thetaRef(ii),'%02.0f') '    '] ;
end
for R0 = [13 14 15]
	xBS=R0:-0.01:0;
	yBS=sqrt(0.04*(xBS-R0).^2-45.3*(xBS-R0)); % original F/G model adds rstandoff^2=645
	rBS = sqrt(xBS.^2 + yBS.^2);
	thetaBS = atan2d(yBS,xBS);
	rRef = interp1(thetaBS,rBS,thetaRef);
	out{1}=[out{1} '======'];
	out{2}=[out{2} 'R=' num2str(R0,'%4.0f') '  '];
	out{3}=[out{3} '======'];
	for ii = 1:numel(thetaRef)
		out{ii+3} = [out{ii+3} num2str(rRef(ii),'%04.2f') ' '];
	end
end
disp('Bow shock, Rin=13, Rout=15');
disp(sprintf('%s\n',out{:}));