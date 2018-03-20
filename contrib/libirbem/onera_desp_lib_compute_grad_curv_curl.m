function [grad_par,grad_perp,grad_drift,curvature,Rcurv,curv_drift,curlB,divB] = onera_desp_lib_compute_grad_curv_curl(Bgeo,B,gradBmag,diffB)
%***************************************************************************************************
% Copyright 2011, T.P. O'Brien
%
% This file is part of IRBEM-LIB.
%
%    IRBEM-LIB is free software: you can redistribute it and/or modify
%    it under the terms of the GNU Lesser General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    IRBEM-LIB is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public License
%    along with IRBEM-LIB.  If not, see <http://www.gnu.org/licenses/>.
%
%***************************************************************************************************
%
% [grad_par,grad_perp,grad_drift,curvature,Rcurv,curv_drift,curlB,divB] = onera_desp_lib_compute_grad_curv_curl(Bgeo,B,gradBmag,diffB)
% returns gradient/curvature force factors and div/curl B
% Bgeo is N x 3, B is N x 1
% gradBmag is N x 3
% diffB is N x 3 x 3 diffB(:,i,j) = dB(:,i)/dx(j)
% all arguments and derivatives are in Cartesian GEO
% grad_par Nx1 gradient of Bmag along B nT/RE
% grad_perp Nx3 ! gradient of Bmag perpendicular to B nT/RE
% grad_drift Nx3 ! (bhat x grad_perp)/B, 1/RE (part of gradient drift velocity)
% curvature Nx3! (bhat dot grad)bhat, 1/RE (part of curvature force)
% Rcurv Nx1 ! 1/|curvature| RE (radius of curvature)
% curv_drift Nx3 ! (bhat x curvature), 1/RE (part of curvature drift)
% curlB Nx3! curl of B (nT/RE) (part of electrostatic current term)
% divB Nx1! divergence of B (nT/RE) (should be zero!)

onera_desp_lib_load;

ntime = numel(B);

Nmax = onera_desp_lib_ntime_max; % maximum array size in fortran library
if ntime > Nmax % break up into multiple calls
    grad_par = nan(ntime,1);
    grad_perp = nan(ntime,3);
    grad_drift = nan(ntime,3);
    curvature = nan(ntime,3);
    Rcurv = nan(ntime,1);
    curv_drift = nan(ntime,3);
    curlB = nan(ntime,3);
    divB = nan(ntime,1);
    for i = 1:Nmax:ntime
        ii = i:min(i+Nmax-1,ntime);
        [grad_par(ii),grad_perp(ii,:),grad_drift(ii,:),curvature(ii,:),Rcurv(ii),curv_drift(ii,:),curlB(ii,:),divB(ii)] = ...
            onera_desp_lib_comput_grad_curv_curl(Bgeo(ii,:),B(ii),gradBmag(ii,:),diffB(ii,:,:));
    end
    return
end

grad_parPtr = libpointer('doublePtr',nan(Nmax,1));
grad_perpPtr = libpointer('doublePtr',nan(3,Nmax));
grad_driftPtr = libpointer('doublePtr',nan(3,Nmax));
curvaturePtr = libpointer('doublePtr',nan(3,Nmax));
RcurvPtr = libpointer('doublePtr',nan(Nmax,1));
curv_driftPtr = libpointer('doublePtr',nan(3,Nmax));
curlBPtr = libpointer('doublePtr',nan(3,Nmax));
divBPtr = libpointer('doublePtr',nan(Nmax,1));
Bgeo = cat(1,Bgeo,nan(Nmax-ntime,3));
B = cat(1,B,nan(Nmax-ntime,1));
gradBmag = cat(1,gradBmag,nan(Nmax-ntime,3));
diffBperm = cat(1,diffB,nan([Nmax-ntime,3,3])); % pad
diffBperm = reshape(diffBperm,[Nmax 9])';

calllib('onera_desp_lib','compute_grad_curv_curl_',ntime,Bgeo',B,gradBmag',diffBperm,...
    grad_parPtr,grad_perpPtr,grad_driftPtr,curvaturePtr,RcurvPtr,curv_driftPtr,curlBPtr,divBPtr);
% have to do this next bit because Ptr's aren't really pointers
grad_par = get(grad_parPtr,'value');
grad_par = grad_par(1:ntime);
grad_perp = get(grad_perpPtr,'value')';
grad_perp = grad_perp(1:ntime,:);
grad_drift = get(grad_driftPtr,'value')';
grad_drift = grad_drift(1:ntime,:);
curvature = get(curvaturePtr,'value')';
curvature = curvature(1:ntime,:);
Rcurv = get(RcurvPtr,'value');
Rcurv = Rcurv(1:ntime);
curv_drift = get(curv_driftPtr,'value')';
curv_drift = curv_drift(1:ntime,:);
curlB = get(curlBPtr,'value')';
curlB = curlB(1:ntime,:);
divB = get(divBPtr,'value');
divB = divB(1:ntime);

% the flag value is actually -1d31
grad_par(grad_par<-1e30) = nan;
grad_perp(grad_perp<-1e30) = nan;
grad_drift(grad_drift<-1e30) = nan;
curvature(curvature<-1e30) = nan;
Rcurv(Rcurv<-1e30) = nan;
curv_drift(curv_drift<-1e30) = nan;
curlB(curlB<-1e30) = nan;
divB(divB<-1e30) = nan;
