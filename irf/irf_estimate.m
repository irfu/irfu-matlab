function out = irf_estimate(what_to_estimate,varargin)
%IRF_ESTIMATE  estimate values for some everyday stuff 
%   out=IRF_ESTIMATE(what_to_estimate,parameters); % 
%
%  Example:
%    out=irf_estimate('capacitance_sphere',r); % estimate capacitance of sphere,
%                       given radius. Everything in SI units.
%
% To see all possibilities execute
%   irf_estimate 

if nargin==0
    disp(' ');
    disp('out=irf_estimate(what_to_estimate,flag);');
    disp('Possible ''what_to_estimate'' values:');
    disp(' ');
    fid=fopen(which('irf_estimate'));
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        if ~isempty(regexp(tline,'case\s''.*''.*$','match'))
            disp(tline(10:end))
        end
    end
    fclose(fid);
    disp(' ');
    return
end
Units=irf_units;
switch what_to_estimate
    case 'capacitance_disk'     % irf_estimate('capacitance_disk',radius)
        radius=varargin{1};
        out=8*Units.eps0*radius;
    case 'capacitance_sphere'   % irf_estimate('capacitance_sphere',radius)
        radius=varargin{1};
        out=4*pi*Units.eps0*radius;
    case 'capacitance_wire'     % irf_estimate('capacitance_wire',radius,length)
        radius=varargin{1};
        length=varargin{2};
				if isempty(radius) || radius == 0 || isempty(length)
					out = [];
				elseif length < 10*radius
					error('capacitance_wire requires length at least 10 times the radius!');
				else
					L=log(length/radius);
					out=2*pi*Units.eps0*length/L*(1+1/L*(1-log(2)));
				end
    case 'capacitance_cylinder' % irf_estimate('capacitance_cylinder',radius,length)
        % Verolino (1995) Electr. Eng. 78, 201-207 Eq. 21-22
        a=varargin{1};   % radius
        h=varargin{2}/2; % half length
        coef=4*pi^2*a*Units.eps0;
        if h/a>.5 && h/a < 4
            C1=pi*h/2./a*1./(log(16*h./a).^2+pi^2/12);
            out=coef*C1;
        elseif h/a >= 4
            Om=2*(log(4*h./a)-1);
            C1=2*h/pi./a.*(1./Om+(4-pi^2)./Om.^3);
            out=coef*C1;
		else
			disp('length less than diameter, do not have formula yet');
            out=NaN;
%        C1=1./(log(16*h./a)-1/2^4*(log(16*h./a)-2).*(h./a).^2); % Eq 20, does not seem ok
        end
    otherwise
        disp(['!!! irf_estimate: unknown estimate ''' lower(what_to_estimate) ''', not estimating.'])
        disp('CONSIDER ADDING THIS ESTIMATE;)!');
        out=in;
end
