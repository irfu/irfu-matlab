function [ax,args,nargs] = axescheck(varargin)
%AXESCHECK Process Axes objects from input list
%   [AX,ARGS,NARGS] = AXESCHECK(ARG1,ARG2,...) looks for Axes provided in
%   the input arguments. It first checks if ARG1 is an Axes. If so, it is
%   removed from the list in ARGS and the count in NARGS.
args = varargin;
nargs = nargin;
ax=[];
if (nargs > 0) && all(all(ishghandle(args{1},'axes')))
  ax = args{1};
  args = args(2:end);
  nargs = nargs-1;
end
end