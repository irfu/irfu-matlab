function [spec_pa,spec_en] = c_peace_spectra(peace)
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
%      spec_pa: Data structure with the following fields
%              t: Time stamps (epoch) 
%              f: energy levels 
%        f_label: energy label
%              p: cell array with spectra (1-parallel, 2-perp, 3-antipar}
%        p_label: cell arrray with spectra label
%             pa: cell array with pitch angles at which measurements were taken (can be
%        important for e.g. 'parallel' option)
%
%    See also C_PEACE_READ_QJAS_CDF
%
% $Id$

% construct the zero result
<<<<<<< c_peace_spectra.m

% constructing spectrograms for given pitch angle and all energies
spec_pa.t=peace.t;
spec_pa.dt=peace.dt;
spec_pa.f=peace.level;
spec_pa.f_unit=peace.level_unit;
spec_pa.f_label=['Energy [' spec_pa.f_unit ']'];
spec_pa.p={zeros(size(peace.psd,1),size(peace.psd,3))};spec_pa.p{2}=spec_pa.p{1};spec_pa.p{3}=spec_pa.p{1};
spec_pa.p_unit=peace.psd_unit;
spec_pa.p_label={['psd par \newline [' spec_pa.p_unit ']'], ['psd perp \newline [' spec_pa.p_unit ']'],['psd antipar \newline [' spec_pa.p_unit ']']};
=======
spec.t=peace.t;
spec.dt=peace.dt;
spec.f=peace.level;
spec.f_unit=peace.level_unit;
spec.f_label=['Energy [' spec.f_unit ']'];
spec.p={zeros(size(peace.psd,1),size(peace.psd,3))};spec.p{2}=spec.p{1};spec.p{3}=spec.p{1};
spec.p_unit=peace.psd_unit;
spec.p_label={['psd par \newline [' spec.p_unit ']'], ['psd perp \newline [' spec.p_unit ']'],['psd antipar \newline [' spec.p_unit ']']};
>>>>>>> 1.4

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
            for jt=1:length(spec_pa.t),
                for jf=1:length(spec_pa.f),
                    spec_pa.p{jj}(jt,jf)=peace.psd(jt,theta_par_index(jt,1,jf),jf);
                end
            end
            spec_pa.pa{jj}=squeeze(theta_par);
        case 'antiparallel'
            for jt=1:length(spec_pa.t),
                for jf=1:length(spec_pa.f),
                    spec_pa.p{jj}(jt,jf)=peace.psd(jt,theta_antipar_index(jt,1,jf),jf);
                end
            end
            spec_pa.pa{jj}=squeeze(theta_antipar);
        case 'perp' % read only 90 degree values
            spec_pa.p{jj}=squeeze(peace.psd(:,7,:));
            spec_pa.pa{jj}=zeros(size(spec_pa.p{jj}));spec_pa.pa{jj}=spec_pa.pa{jj}+90;
    end
end


% constructing spectrograms for given energy and all pitch angles
spec_en.t=peace.t;
spec_en.dt=peace.dt;
spec_en.f=peace.theta;
spec_en.df=peace.theta_delta;
spec_en.f_unit='[deg]';
spec_en.f_label=['Pitch angle [deg]'];
spec_en.p={zeros(size(peace.psd,1),size(peace.psd,2))};
spec_en.p_unit=peace.psd_unit;
spec_en.p_label={['psd par \newline [' spec_en.p_unit ']'], ['psd perp \newline [' spec_en.p_unit ']'],['psd antipar \newline [' spec_en.p_unit ']']};

% constract measured pitch angle matrix (pitch angle for every time and energy)
for jen=1:length(peace.level),
    spec_en.p{jen}=spec_en.p{1};
    spec_en.p_label{jen}=['E=' num2str(peace.level(jen),2) 'eV \newline ' spec_en.p_unit ];
    spec_en.p{jen}=peace.psd(:,:,jen);
end
