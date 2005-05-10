function plan = caa_pl_coverage(onoff,cover,tint)
%CAA_PL_COVERAGE  pl on/off and coverage information
%
% caa_pl_coverage(onoff,cover)
%
% $Id$

% Copyright 2005 Yuri Khotyaintsev

DT=20;

clf
axes;
set(gca,'YLim',[0 5]);
hold on
for cl_id=1:4
	oo = onoff{cl_id};
	co = cover{cl_id};
	io = find(oo(:,1)>=tint(1) & oo(:,1)<=tint(2));
	if isempty(io)
		io = find(oo(:,1)<tint(1));
		if isempty(io)
			io = find(oo(:,1)>tint(2));
			io = io(1);
		else, io = io(end);
		end
	end
	if oo(io(1),2)==0
		if io(1)>1, io=[io(1)-1; io]; end
	end
	if oo(io(end),2)==1
		if io(end)<length(oo(:,1)), io=[io; io(end)+1]; end
	end
	
	ic = find(co(:,2)>=tint(1) & co(:,1)<=tint(2));

	for j=1:2:length(io)-1
		irf_plot([[oo(io(j),1) oo(io(j+1),1)]; (5-cl_id)*[1 1]]','g-x')
	end
	for j=1:length(ic)
		if co(ic(j),3)==1
			irf_plot([[co(ic(j),1) co(ic(j),2)]; (5-cl_id)*[1 1]-.2]','r-x')
		else
			irf_plot([[co(ic(j),1) co(ic(j),2)]; (5-cl_id)*[1 1]-.2]')
		end
	end
	
end
hold off
set(gca,'YTick',[1 2 3 4],'YTickLabel',['C4'; 'C3'; 'C2'; 'C1'])
irf_zoom(tint,'x',gca)

if nargout<1, return, end
% Planning intervals
for cl_id=1:4
	co = cover{cl_id};
	ii = find(co(:,1)>=tint(1) & co(:,2)<=tint(2));
	co = co(ii,:);
	
	% Remove all NM intervals which are shorter than 10 min
	ii = find(co(:,2)-co(:,1)<600 & co(:,3)==0);
	if ~isempty(ii)
		disp(['throwing away ' num2str(length(ii)) ' short intervals'])
		co(ii,:) = [];
	end
	cover_tmp(cl_id) = {co};
end
plan = [];
% We first find 4-7 min intervals of BM1 which are indicative of starting 
% long and good intervals of NM.
cl_id = 1;
while cl_id<=4
	co = cover_tmp{cl_id};
	dt = co(:,2)-co(:,1);
	ii = find(dt<7*60 & dt>4*60 & co(:,3)==1);
	if ~isempty(ii)
		ii = fliplr(torow(ii));
		for j=ii
			if find_close(co(j,1),1,cover_tmp,cl_id)
				if j<length(co(:,1))
					% We are really interested in the interval which follows this one
					ts = co(j+1,1);
					te = co(j+1,2);
					sc_list = 1:4;
					sc_list(find(sc_list==cl_id)) = [];
					for cli=sc_list
						co_tmp = cover_tmp{cli};
						ii_tmp = find(co_tmp(:,3)==0);
						co_tmp = co_tmp(ii_tmp,1:2);
						ii_tmp = find(co_tmp(:,1)>ts-DT & co_tmp(:,1)<ts+DT);
						if ~isempty(ii_tmp)
							if co_tmp(ii_tmp,1)>ts, ts = co_tmp(ii_tmp,1); end
							if co_tmp(ii_tmp,2)>te, te = co_tmp(ii_tmp,2); end
							disp(['including C' num2str(cli) ' ' epoch2iso(co_tmp(ii_tmp,1)) ' - ' epoch2iso(co_tmp(ii_tmp,2))])
							% clear the interval
							ii_tmp = find(cover_tmp{cli}(:,1)==co_tmp(ii_tmp,1));
							cover_tmp{cli}(ii_tmp-1:ii_tmp,:) = [];
						end
					end
					plan = [plan; [ts te]];
					disp(['saving ' epoch2iso(ts) ' - ' epoch2iso(te)])
					disp(['1 clearing C' num2str(cl_id) ' ' epoch2iso(co(j,1)) ' - ' epoch2iso(co(j,2))])
					% clear the intervals
					co(j:j+1,:) = [];
				else
					% this was the last interval, throw it away anyway
					disp(['2 throwing away C' num2str(cl_id) ' ' epoch2iso(co(j,1)) ' - ' epoch2iso(co(j,2))])
					cover_tmp = throw_away(co(j,1),1,cover_tmp,cl_id);
					co(j,:) = [];
				end
			else
				disp(['3 throwing away C' num2str(cl_id) ' ' epoch2iso(co(j,1)) ' - ' epoch2iso(co(j,2))])
				co(j,:) = [];
			end
		end
	end
	cover_tmp(cl_id) = {co};
	cl_id = cl_id + 1;
end

if ~isempty(plan)
	hold on
	for j=1:length(plan(:,1))
		irf_plot([[plan(j,1) plan(j,2)]; 4*[1 1]+.2]','k-x')
	end
	hold off
end

function res=find_close(ts,tmmode,cover,scn)
	DT = 20; % 20 seconds
	sc_list = 1:4;
	sc_list(find(sc_list==scn)) = [];
	res = 0;
	for cli=sc_list
		co = cover{cli};
		ii = find(co(:,3)==tmmode);
		co = co(ii,1);
		if ~isempty(find(co>ts-DT & co<ts+DT)), res = 1; break, end
	end
	
function cover_new = throw_away(ts,tmmode,cover,scn)
	DT = 20; % 20 seconds
	sc_list = 1:4;
	sc_list(find(sc_list==scn)) = [];
	res = 0;
	cover_new = cover;
	for cli=sc_list
		co = cover{cli};
		ii = find(co(:,1)>ts-DT & co(:,1)<ts+DT & co(:,3)==tmmode);
		if ~isempty(ii)
			disp(['4 throwing away C' num2str(cli) ' ' epoch2iso(co(ii,1)) ' - ' epoch2iso(co(ii,2))])
			co(ii,:) = [];
		end
		cover_new(cli) = {co};
	end
