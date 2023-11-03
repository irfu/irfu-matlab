function [t1_d,t2_d,t1_dd,t2_dd]=irf_corr_deriv(x1,x2,fla)
%IRF_CORR_DERIV  Correlate the derivatives of two time series
%
% [t1_d,t2_d,t1_dd,t2_dd]=irf_corr_deriv(x1,x2,[flag]);
%      Find time instants of common steapest gradients [default] or
%      zero crossings [flag=1] and maxima and minima.
%
% x1,x2 - signals to correlate, first column is time, 2nd is signal
% t1_d,t2_d - positions of maxima and minima
% t1_dd,t2_dd - time instants of common highest
%               first derivatives (steepest gradients) [default]
%             - times of zero crossings [flag=1]
%

if nargin<3, fla = 0; end

k=1:2;

%%%%%%%%%%%%%%% 1st derivative %%%%%%%%%%%

c_eval('dx?=[x?(1:end-1,1)+0.5*diff(x?(:,1)) diff(x?(:,2))];',k)
c_eval('ind_zeros?=find(sign(dx?(1:end-1,2).*dx?(2:end,2))<0);',k)
c_eval('if ind_zeros?(1)==1, ind_zeros?(1)=[];end',k)
c_eval('ind_zeros?_plus=find(dx?(ind_zeros?-1,2)-dx?(ind_zeros?,2)>0);',k)
c_eval('ind_zeros?_minu=find(dx?(ind_zeros?-1,2)-dx?(ind_zeros?,2)<0);',k)
c_eval('ind=ind_zeros?(ind_zeros?_plus);xx=dx?;t_zeros?_plus=xx(ind,1)+(xx(ind+1,1)-xx(ind,1)).*1./(1+abs(xx(ind+1,2))./abs(xx(ind,2)));',k)
c_eval('ind=ind_zeros?(ind_zeros?_minu);xx=dx?;t_zeros?_minu=xx(ind,1)+(xx(ind+1,1)-xx(ind,1)).*1./(1+abs(xx(ind+1,2))./abs(xx(ind,2)));',k)

%Remove repeating points
c_eval('t_zeros?_plus(find(diff(t_zeros?_plus)==0)) = [];',k)

% Define identical pairs of two time axis
[t1_d_plus,t2_d_plus]=irf_find_closest(t_zeros1_plus,t_zeros2_plus); %#ok<ASGLU>
[t1_d_minu,t2_d_minu]=irf_find_closest(t_zeros1_minu,t_zeros2_minu); %#ok<ASGLU>

c_eval('t?_d = sortrows([t?_d_plus;t?_d_minu]);',k)

if fla
  %%%%%%%%%%%%%%% zero crossings %%%%%%%%%%%
  c_eval('ind_zeros?=find(sign(x?(1:end-1,2).*x?(2:end,2))<0);ind_zeros?(ind_zeros?==1)=[];',k)
  c_eval('ind_zeros?_plus=find(x?(ind_zeros?-1,2)-x?(ind_zeros?,2)>0);',k)
  c_eval('ind_zeros?_minu=find(x?(ind_zeros?-1,2)-x?(ind_zeros?,2)<0);',k)
  c_eval('ind=ind_zeros?(ind_zeros?_plus);xx=x?;t_zeros?_plus=xx(ind,1)+(xx(ind+1,1)-xx(ind,1)).*1./(1+abs(xx(ind+1,2))./abs(xx(ind,2)));',k)
  c_eval('ind=ind_zeros?(ind_zeros?_minu);xx=x?;t_zeros?_minu=xx(ind,1)+(xx(ind+1,1)-xx(ind,1)).*1./(1+abs(xx(ind+1,2))./abs(xx(ind,2)));',k)
else
  %%%%%%%%%%%%%%% 2nd derivative %%%%%%%%%%%

  c_eval('ddx?=[dx?(1:end-1,1)+0.5*diff(dx?(:,1)) diff(dx?(:,2))];',k)
  c_eval('ind_zeros?=find(sign(ddx?(1:end-1,2).*ddx?(2:end,2))<0);ind_zeros?(ind_zeros?==1)=[];',k)
  c_eval('ind_zeros?_plus=find(ddx?(ind_zeros?-1,2)-ddx?(ind_zeros?,2)>0);',k)
  c_eval('ind_zeros?_minu=find(ddx?(ind_zeros?-1,2)-ddx?(ind_zeros?,2)<0);',k)
  c_eval('ind=ind_zeros?(ind_zeros?_plus);xx=ddx?;t_zeros?_plus=xx(ind,1)+(xx(ind+1,1)-xx(ind,1)).*1./(1+abs(xx(ind+1,2))./abs(xx(ind,2)));',k)
  c_eval('ind=ind_zeros?(ind_zeros?_minu);xx=ddx?;t_zeros?_minu=xx(ind,1)+(xx(ind+1,1)-xx(ind,1)).*1./(1+abs(xx(ind+1,2))./abs(xx(ind,2)));',k)
end

% Define identical pairs of two time axis
[t1_dd_plus,t2_dd_plus]=irf_find_closest(t_zeros1_plus,t_zeros2_plus); %#ok<ASGLU>
[t1_dd_minu,t2_dd_minu]=irf_find_closest(t_zeros1_minu,t_zeros2_minu); %#ok<ASGLU>

c_eval('t?_dd = sortrows([t?_dd_plus;t?_dd_minu]);',k)
