function out=c_pl_tx(varargin)
%C_PL_TX    Plot data from all four Cluster spacecraft in the same plot with right Cluster colors
%   c=c_pl_tx(x1,x2,x3,x4,[column,t_unit_in_original_units,t_origo_in_original_units])
%   c=c_pl_tx(x1,x2,x3,x4,[column,1,[dt1 dt2 dt3 dt4]])
%   c=c_pl_tx(x1,x2,x3,x4,[column,t_unit_in_original_units,t_origo_in_original_units])
%   c=c_pl_tx('x?',column,1,[dt1 dt2 dt3 dt4])
%   plot all 4 cluster values, time is in the first column, <column> gives which column to plot
%   if column is vector then create subplots so that in each subplot is corresponding column
%
%   Example:
%      c_pl_tx('av_abs(B?)') - will plot 3 components + magnitude of B1:4.
%
%   $Id$

error(nargchk(1,8,nargin))

args = varargin;

if isstr(args{1})
	% We have 4 arguments
	for cl_id=1:4
		ttt = evalin('caller',av_ssub(args{1},cl_id)); 
		eval(av_ssub('x? =ttt;',cl_id)); clear ttt
	end
	if length(args) > 1, args = args(2:end); 
	else, args = ''; end
else
	% We have 8 arguments
	c_eval('x? = args{?};');
	if length(args) > 4, args = args(5:end); 
	else, args = ''; end
end

if length(args) < 1
	% try to guess the size of the matrix
	column = size(x1,2);
	if column > 2, column = 2:column; end
else
	column = args{1};
end

if length(args) < 2, t_unit_in_original_units=1;
else, t_unit_in_original_units = args{2};
end
if length(args) < 3, t_origo_in_original_units=0;
else, t_origo_in_original_units = args{3};
end

if (length(t_origo_in_original_units) == 4)
 tt=t_origo_in_original_units;ts1=tt(1);ts2=tt(2);ts3=tt(3);ts4=tt(4);
else
 ts1=t_origo_in_original_units(1);ts2=ts1;ts3=ts1;ts4=ts1;
end
tu=t_unit_in_original_units;
ts=t_origo_in_original_units;

if length(column) == 1,
  h=plot((x1(:,1)-ts1)/tu,x1(:,column),'k',(x2(:,1)-ts2)/tu,x2(:,column),'r',(x3(:,1)-ts3)/tu,x3(:,column),'g',(x4(:,1)-ts4)/tu,x4(:,column),'b');
  c=get(h(1),'parent');
  grid on
else
  clf;
  for j=1:length(column),
    c(j)=av_subplot(length(column),1,-j);
    plot((x1(:,1)-ts1)/tu,x1(:,column(j)),'k',(x2(:,1)-ts2)/tu,x2(:,column(j)),'r',(x3(:,1)-ts3)/tu,x3(:,column(j)),'g',(x4(:,1)-ts4)/tu,x4(:,column(j)),'b');
    grid on
  end
end

if x1(1,1) > 9e8,  add_timeaxis(c); end

if nargout > 0, out = c; end


