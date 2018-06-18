function fout=pitch_angle_average(fin,theta,thetaRange,deltaTheta,thetaDimension)
% IRF.PITCH_ANGLE_AVERAGE average data over pitch angle range
%
%	fout = IRF.PITCH_ANGLE_AVERAGE(fin,theta,thetaRange,[deltaTheta],[thetaDimension])
%		f - matrix to average
%		theta - pitch angle vector in degrees
%		thetaRange - theta range to average [theta_min theta_max] (default all)
%		deltaTheta - step in pitch angle (if not given, calculate)
%		thetaDimension - which dimension of matrix f correspond to pitch angle
%
%	fout matrix has one dimension less than fin, the pitch angle dimension is removed.
%

if nargin < 3 || isempty(thetaRange)
	thetaRange = [min(theta(:)) max(theta(:))];
end
if nargin < 4 || isempty(deltaTheta) % define dtheta
	deltaTheta = (theta(2)-theta(1))/2;
end
if nargin < 5 % define theta dimension
	szind = (size(fin) == numel(theta));
	if ~any(szind)
		errStr = 'Pitch angle dimension cannot be identified!';
		irf.log('critical',errStr);
		error(['irf.pitch_angle_average: ' errStr]);
	elseif sum(szind) > 1 % more than 1 dimension corresponds theta vector size
		irf.log('warning','WARNING!!! More than 1 dimension is of the theta vector size, assuming last one is pitch angle');
		thetaDimension = find(szind,1,'last');
	else
		thetaDimension = find(szind);
	end
end

indTheta = (theta >= thetaRange(1) & theta <= thetaRange(2));
dataExist = ~isnan(fin);
thetaSteradian = bsxfun(@plus,zeros(size(fin)),...
	indTheta.*(cosd(theta-deltaTheta)-cosd(theta+deltaTheta)));
thetaSteradian(~dataExist) = NaN;

totalSteradian = irf.nansum(thetaSteradian,thetaDimension);

ftemp = bsxfun(@times,fin,thetaSteradian);
fout = irf.nansum(ftemp, thetaDimension)./totalSteradian;

