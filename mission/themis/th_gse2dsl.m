function out = th_gse2dsl(inp,thId)
%TH_GSE2DSL  Trasform THEMIS data from GSE to DSL
%
%  vec_DSL = th_gse2dsl(vec_GSE,thId)
%
% NOTE: requires l1/state files in /data/themis

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if ~any(thId=='abcde'), error('THEMIS ID must be a, b, c, d, e'), end


t = inp(:,1);

spinras=th_read_l2(['th' thId '_spinras'],t([1 end])');
spindec=th_read_l2(['th' thId '_spindec'],t([1 end])');

[xspin,yspin,zspin] = sph2cart(spinras(:,2)*pi/180,spindec(:,2)*pi/180,1);
geiSAX = [spinras(:,1) xspin yspin zspin];

sax = irf.geocentric_coordinate_transformation(geiSAX,'gei>gse');
sax = irf_resamp(sax,t);
Rx = sax(:,2);
Ry = sax(:,3);
Rz = sax(:,4);
a = 1./sqrt(Ry.^2+Rz.^2);
M = zeros(length(a),3,3);
M(:,1,1) = a.*(Ry.^2+Rz.^2);
M(:,1,2) = -a.*Rx.*Ry;
M(:,1,3) = -a.*Rx.*Rz;
M(:,2,1) = 0;
M(:,2,2) = a.*Rz;
M(:,2,3) = -a.*Ry;
M(:,3,1) = Rx;
M(:,3,2) = Ry;
M(:,3,3) = Rz;
out = inp;
out(:,2:4) = mult_mat(M,inp(:,2:4));
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