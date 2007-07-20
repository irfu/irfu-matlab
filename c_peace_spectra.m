function [spec,pa] = c_peace_spectra(peace)
%C_PEACE_SPECTRA  construct different spectra (parallel, perp, antiparallel, ..)
%
% spec = c_peace_spectra(peace)
%
% Read PEACE structure obtained from C_PEACE_READ_QJAS_CDF
%
% Input: 
%     peace: PEACE structure obtained from C_PEACE_READ_QJAS_CDF
%    option:
%      'parallel' - get pitch angles closest to parallel
%  'antiparallel' - get pitch angles closest to parallel
%          'perp' - get pitch angles closest to parallel
%
% Output: 
%      spec: Data structure with the following fields
%              t: Time stamps (epoch) 
%              f: energy levels 
%        f_label: energy label
%              p: cell array with spectra (1-parallel, 2-perp, 3-antipar}
%        p_label: cell arrray with spectra label
%          
%        pa: cell array with pitch angles at which measurements were taken (can be
%        important for e.g. 'parallel' option)
%
%    See also C_PEACE_READ_QJAS_CDF
%
% $Id$

% construct the zero result
spec.t=peace.t;
spec.f=peace.level;
spec.f_unit=peace.level_unit;
spec.f_label=['Energy [' spec.f_unit ']'];
spec.p={zeros(size(peace.psd,1),size(peace.psd,3))};spec.p{2}=spec.p{1};spec.p{3}=spec.p{1};
spec.p_unit=peace.psd_unit;
spec.p_label={['psd par [' spec.p_unit ']'], ['psd perp [' spec.p_unit ']'],['psd antipar [' spec.p_unit ']']};

% constract measured pitch angle matrix (pitch angle for every time and energy)
theta=peace.theta;
theta_all=peace.psd;
for j=1:length(theta),
    theta_all(:,j,:)=peace.theta(j);
end
theta_all(find(~peace.psd))=NaN; % fast(dirty) solution assuming zero counts means pitch angle no measured
[theta_par,theta_par_index]=min(theta_all,[],2); 
[theta_antipar,theta_antipar_index]=max(theta_all,[],2);

peace_options={'parallel','perp','antiparallel'};
for jj=1:3;
    switch lower(peace_options{jj})
        case 'parallel'
            for jt=1:length(spec.t),
                for jf=1:length(spec.f),
                    spec.p{jj}(jt,jf)=peace.psd(jt,theta_par_index(jt,1,jf),jf);
                end
            end
            pa{jj}=squeeze(theta_par);
        case 'antiparallel'
            for jt=1:length(spec.t),
                for jf=1:length(spec.f),
                    spec.p{jj}(jt,jf)=peace.psd(jt,theta_antipar_index(jt,1,jf),jf);
                end
            end
            pa{jj}=squeeze(theta_antipar);
        case 'perp' % read only 90 degree values
            spec.p{jj}=squeeze(peace.psd(:,7,:));
            pa{jj}=zeros(size(spec.p{jj}));pa{jj}=pa{jj}+90;
    end
end
