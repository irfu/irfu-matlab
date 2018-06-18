function [out]=irf_newxyz(inp,x,y,z)
%IRF_NEWXYZ   Rotate vector into a new coordinate system
%
% [out]=irf_newxyz(inp,x,y,z)
% inp,out - column vector if more than 3 columns assume that first is time and 2nd-4th is X Y Z
% x,y,z - vectors (x=[xx xy xz], y= [yx yy yz]; )
%          that give new coordinates, if some is 0 then calculate it from other
% if x=0 -> x=yxz,z=xxy
% if y=0 -> y=zxx,x=yxz
% if z=0 -> z=xxy,y=zxx

if x==0
 x=cross(y,z);
 z=cross(x,y);
end
if y==0
 y=cross(z,x);
 x=cross(y,z); % to make sure x and z are orthogonal
end
if z==0
 z=cross(x,y);
 y=cross(z,x);
end

x=x/norm(x);
y=y/norm(y);
z=z/norm(z);

out=inp;
if isa(inp,'TSeries')
  % XXX: this must be moved to TSeries
  % and here we should give an error
  if inp.tensorOrder==0, error('Does not work on scalars'),
  elseif inp.tensorOrder>1, error('Not implemented yet')
  end
  if strcmp(inp.tensorBasis(1:3),'xyz')
    out.data(:,1)=inp.data(:,1:3)*x';
    out.data(:,2)=inp.data(:,1:3)*y';
    out.data(:,3)=inp.data(:,1:3)*z';
    out.coordinateSystem = 'new_xyz';
  else, error('Not implemented yet')
  end
  return
end
  
if size(out,2)==3
  out(:,1)=inp(:,1:3)*x';
  out(:,2)=inp(:,1:3)*y';
  out(:,3)=inp(:,1:3)*z';
elseif size(out,2)>3
  out(:,2)=inp(:,2:4)*x';
  out(:,3)=inp(:,2:4)*y';
  out(:,4)=inp(:,2:4)*z';
else
	errStr = 'Coordinate transformation not possible when input has less than 3 columns';
	irf.log('critical',errStr);
	error('irf_newxyz:too_few_components',errStr);
end
