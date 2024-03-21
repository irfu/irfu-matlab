function TsOut = mms_dsl2gse(TsIn, defatt, direction)
%MMS_DSL2GSE  transform TS from DSL to GSE
%
% TsOut = mms_dsl2gse(TsIn, Defatt, [direction])
% TsOut = mms_dsl2gse(TsIn, spin_axis_gse, [direction])
%
% Transfrorms a vector from DSL to GSE (of GSE->DSL for direction=-1)
%
% Example usage with DEFATT:
%  tint = irf.tint('2015-05-09T14:00:00Z/2015-05-09T17:59:59Z');
%  defatt = mms.db_get_variable('mms2_ancillary_defatt','zra',tint);
%  defatt.zdec = mms.db_get_variable('mms2_ancillary_defatt','zdec',tint).zdec;
%  gseB = mms_dsl2gse(B_dmpa,defatt)
%
% Example usage with SAX:
%  sax = [0 0.7071 0.7071];
%  gseB = mms_dsl2gse(B_dmpa,sax);

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if nargin<3, direction = 1; end
if abs(direction)~=1
  direction = 1;
  irf.log('warning','using GSE->DSL')
end

defatt = mms_removerepeatpnts(defatt);

if isa(defatt,'TSeries')
  errS = 'Not implemented yet';
  irf.log('critical',errS), error(errS)
elseif isstruct(defatt) && isfield(defatt,'time') && ...
    isfield(defatt,'zra')  && isfield(defatt,'zdec')
  ttDefatt = EpochTT(defatt.time);
  [x,y,z] = sph2cart(defatt.zra*pi/180,defatt.zdec*pi/180,1);
  saxTmp = irf.geocentric_coordinate_transformation(...
    [ttDefatt.epochUnix x y z],'gei>gse');
  spin_axis = saxTmp(:,2:4);

  %XXX: FIXME add quaternion interpolation
  tData = ttDefatt - ttDefatt(1); data = double(spin_axis);
  newData = irf_resamp([tData data],TsIn.time-ttDefatt(1));
  spin_axis = newData(:,2:end);
elseif isvector(defatt) && length(defatt)==3
  spin_axis = defatt;
else
  errS = 'unrecognized DEFATT/SAX input';
  irf.log('critical',errS), error(errS)
end

Rx = spin_axis(:,1); Ry = spin_axis(:,2); Rz = spin_axis(:,3);
a = 1./sqrt(Ry.^2+Rz.^2);
M = zeros(length(a),3,3);
M(:,1,:) = [a.*(Ry.^2+Rz.^2) -a.*Rx.*Ry -a.*Rx.*Rz];
M(:,2,:) = [0*a a.*Rz	-a.*Ry];
M(:,3,:) = [Rx	Ry Rz];

TsOut = TsIn.transform('xyz');
inp = TsOut.data;

if direction == -1
  out = mult_mat(M,inp);
elseif direction==1
  out = mult_mat(transpose_mat(M),inp);
else
  disp('No coordinate transformation done!')
end

TsOut.data = out;

end

function out=transpose_mat(inp)
out = inp;
if numel(size(inp))==2 || (size(inp,2)~=size(inp,3))
  error('not impemented');
end
for ii=1:size(inp,2)
  for jj=ii+1:size(inp,2)
    out(:,ii,jj) = inp(:,jj,ii);
    out(:,jj,ii) = inp(:,ii,jj);
  end
end
end

function out = mult_mat(inp1,inp2)
dimInp1 = numel(size(inp1));
dimInp2 = numel(size(inp2));
if (dimInp1==dimInp2)
  numOfMult=size(inp1,3);
  T = zeros(size(inp1,1),size(inp1,2),size(inp2,3));
  for ii=1:size(inp1,2)
    for jj=1:size(inp2,3)
      for kk=1:numOfMult
        T(:,ii,jj)=T(:,ii,jj)+inp1(:,ii,kk).*inp2(:,kk,jj);
      end
    end
  end
elseif (dimInp1==3) && (dimInp2==2)
  numOfOutp=size(inp1,2);
  numOfInp=size(inp2,2);
  T = inp2(:,1)*zeros(1,numOfOutp);
  for ii=1:numOfOutp
    for jj=1:numOfInp
      T(:,ii)=T(:,ii)+inp1(:,ii,jj).*inp2(:,jj);
    end
  end
end
out = T;
end