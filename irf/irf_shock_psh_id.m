function PSH = irf_shock_psh_id(shock_data,filter_params,method_params)
%
% IRF_SHOCK_PSH_ID - Identify Phase Space Holes (PSHs) in the reduced ion
% distribution.
%
%   PSH = IRF_SHOCK_PSH_ID(shock_data,filter_params,method_params) Find PSHs in
%   the reduced ion distribution which is contained in the structure
%   shock_data. Input filter_params and method_params define methods for
%   finding PSHs.
%
%   Input parameters:
%
%   shock_data      - Structure containing various data from the shock
%                   crossing
%     > f1Dn      - PDist object with reduced distribution
%     > others?   - (not used)
%  	filter_params   - Filter parameters in a structure
%     > method    - Filter method, either 'gaussian'
%                   (default) or 'runav' (running average)
%     > fw        - Filter values [wt, wv], default is [1,2]
%  	method_params   - Method parameters in a structure
%     > method    - 'contour', 'boundary' , or 'both'
%     > cl        - Number of contour levels (or actual levels in array?)
%     > bl        - Criteria limits set on PSHs in the
%                   boundary method [FminRel,FmeanRel,FmaxRel,Fmin],
%                   default is [1.25,1.5,0,10]
%
%   Output fields in structure:
%    	(C denotes contour, B denotes boundary)
%    	> CN/BN     -   Number of holes detected by each method
%    	> Cxy/Bxy   -	Cell array with time-velocity data for all PSHs
%    	> CTS/BTS   - 	TSeries objects of the holes for easy plotting
%                     	(see example)
%    	> Cparams   -   INSERT EXPLANATION
%    	> filter_params
%    	> method_params
%
%   	DERIVED PROPERTIES LIKE WIDTH AND HEIGHT ARE NOT IMPLEMENTED
%
%
%   Examples:
%    	%%% get and plot PHSs for MMS1 for (Johlander et al. 2018) %%%
%    	tint = irf.tint('2016-01-06T00:32:20/2016-01-06T00:33:30');
%    	% read FPI-DIS data and remove 1-count level
%    	iPDist = mms.get_data('PDi_fpi_brst_l2',tint,1);
%   	iPDistErr = mms.get_data('PDERRi_fpi_brst_l2',tint,1);
%   	iPDist.data(iPDist.data<1.1*iPDistErr.data) = 0;
%     % reduce distribution
%     nvec = [0.95, -0.30, -0.05]; nvec = nvec/norm(nvec);
%     f1Dn = iPDist.reduce('1D',nvec,'vg',linspace(-800,400,100),'nMC',100);
%     % set input
%     shock_data = struct('f1Dn',f1Dn);
%     filter_params = struct('method','runav','fw',[2,4]);
%     method_params = struct('method','both','B_criteria_values',[1.2,1.5,0,10]);
%     % find PSHs (takes some time)
%     PSH = irf_shock_psh_id(shock_data,filter_params,method_params);
%     % plot
%     hca = irf_plot(1,'newfigure');
%     irf_spectrogram(hca,f1Dn.specrec)
%     hold on % stop crying your heart out
%     % plot first PSHs so legends behave
%     irf_plot(hca,PSH.CTS{1},'k','linewidth',2); irf_plot(hca,PSH.BTS{1},'r','linewidth',2)
%     % plot the rest
%     for ii = 2:PSH.CN; irf_plot(hca,PSH.CTS{ii},'k','linewidth',2); end
%     for ii = 2:PSH.BN; irf_plot(hca,PSH.BTS{ii},'r','linewidth',2); end
%     legend(hca,'PSD','Contour','Boundary')

% TODO:
%   - Clean up code and comment a bit
%   - Should this function calculate width, height, etc.?
%   - Pre-allocate arrays (!)
%   - Print out progress at some points in the program
%   - optimization?

% Written by A. P. Dimmock & A. Johlander


%% Handle input

%%%% filter input
if isfield(filter_params,'method')
  filterMethod = filter_params.method;
else
  filterMethod = 'gaussian'; % default value
  % save
  filter_params.method = filterMethod;
end

% filtering width
switch filterMethod
  case 'gaussian'
    if isfield(filter_params,'fw')
      fw = filter_params.fw;
    else
      fw = [1,2]; % default value
      % save
      filter_params.fw = fw;
    end
  case 'runav'
    if isfield(filter_params,'fw')
      fw = filter_params.fw;
    else
      fw = [1,2];% default value
      % save
      filter_params.fw = fw;
    end
  otherwise
    error('Filter method not implemented')
end


%%%% method input
if isfield(method_params,'method')
  pshMethod = method_params.method;
else
  pshMethod = 'both'; % default value
  method_params.method = pshMethod;
end

% other input parameters
if strcmpi(pshMethod,'contour') || strcmpi(pshMethod,'both') % contour
  % set contour levels
  if isfield(method_params,'cl')
    cl = method_params.cl;
  else
    cl = 50; % default value
    method_params.cl = cl;
  end
end
if strcmpi(pshMethod,'boundary') || strcmpi(pshMethod,'both') % boundary
  % many criteria
  if isfield(method_params,'bl')
    BoundaryCriteria = method_params.bl;
  else
    BoundaryCriteria = [1.25,1.5,0,10];% default value
    method_params.bl = BoundaryCriteria;
  end
end

% should throw some error if unknown PSH method


%% Filter

% Gaussian filter
switch filterMethod
  case 'gaussian'
    % filtered distribution
    f1Dn = filt_mat_gaussian(shock_data.f1Dn,fw(1),fw(2));
  case 'runav'
    f1Dn = filt_mat_runav(shock_data.f1Dn,fw(1),fw(2));
end


t = f1Dn.time.epochUnix;  % time in epoch unix
Vn = f1Dn.depend{1,1}(1,:); % Velocity along the shock normal
F = f1Dn.data'; % filtered phase-space density matrix (transposed for some reason)


%% Contour method

%%%% Contour-based method %%%%
if strcmpi(pshMethod,'contour') || strcmpi(pshMethod,'both')

  %%% arrange data to [t Vn M]
  data = ones(length(t)*length(Vn),3)*nan;    % pre-allocate
  c=0;
  for k=1:length(t)
    for kk=1:length(Vn)
      c = c+1;
      data(c,:) = [t(k),Vn(kk),F(kk,k)];
    end
  end
  clear c

  %%% Make data to grid format - filtered data (inefficient)
  %   xv    = linspace(min(data(:,1)), max(data(:,1)), 2000);     % t values
  %   yv    = linspace(min(data(:,2)), max(data(:,2)), 2000);     % Vn values
  %   [X,Y] = meshgrid(xv, yv);                                   % grid of t-Vn
  %   Z     = griddata(data(:,1),data(:,2),data(:,3),X,Y);        % find M-grid

  % Optimized
  [X,Y] = meshgrid(t, Vn);
  Z = F;

  %%% Map the contours of the filtered distribution
  [~, cm] = contourf(X, Y, Z,cl,'edgecolor','k');         % make a contour plot
  CA      = irf_extract_contour_lines(cm);                % extract the contour lines


  for k=1:length(CA) % run for each contour line

    CSS = [CA{k}(1,1),CA{k}(1,2);CA{k}(end,1),CA{k}(end,2)];        % start and end point of the contour line

    if range(CA{1,k}(:,1))<30 % only do this if contour is <30 seconds wide, saves time
      pinc = inpolygon(data(:,1),data(:,2),CA{k}(:,1),CA{k}(:,2));    % index of points inside the contour from raw data

      clevel = mean(interp2(X,Y,Z,(CA{k}(:,1)),(CA{k}(:,2))));        % the level/edge of the contour

      if length(find(pinc))<2 % if contour is small and not enough points inside, threshold set at 2 points
        min_val = double(interp2(X,Y,Z,mean(CA{k}(:,1)),mean(CA{k}(:,2)))); % value at contour center
        min_yes = double(min_val<clevel); % value at contour center is less than the contour level (for small contours)
      else
        min_yes =  double(100*length(find(data(pinc,3)<clevel))/length(find(pinc))>60); % over 60% of points are below the contour level
      end

      CSinfo(k,:) = [...                                  % information about the specific contour
        sum(diff(CSS,1)),...                            % if 0 then the contour is closed
        range(CA{1,k}(:,1)),...                         % width of the contour
        range(CA{1,k}(:,2)),...                         % height of the contour
        mean(CA{1,k}(:,1)),mean(CA{1,k}(:,2)),...       % center point of the contour
        min_yes,...                                     % if 1, then we classify the contour as a minima
        ];

    end

  end


  %%% find how many center points like within a contour
  rr = ones(length(CA),1);
  for k = 1:length(CA)
    rr(k,1) = length(find(inpolygon(CSinfo(:,4),CSinfo(:,5),CA{k}(:,1),CA{k}(:,2))));
  end
  CSinfo = [CSinfo,rr]; % add this information to the CSinfo matrix


  %%%% now we can start to remove some contours
  CA2 = CA; %%% we will gradually eliminate contours from CA2


  %%%% remove the open contours
  CA2(CSinfo(:,1)~=0)       = [];
  CSinfo(CSinfo(:,1)~=0,:)  = [];
  %%%% remove the contours that are maximias
  CA2(CSinfo(:,6)~=1 )      = [];
  CSinfo(CSinfo(:,6)~=1,:)  = [];
  %%%% remove the contours that are too narrow
  %%%%     - based on outliers of >100 test events
  CA2(CSinfo(:,2)<0.2)      = [];
  CSinfo(CSinfo(:,2)<0.2,:) = [];
  %%%% remove the contours that are too wide
  %%%%     - based on outliers of >100 test events
  CA2(CSinfo(:,2)>10)      = [];
  CSinfo(CSinfo(:,2)>10,:) = [];

  % find if centre points are not within the contour
  for k=1:length(CA2)
    outcen(k) = ~inpolygon(CSinfo(k,4),CSinfo(k,5),CA2{k}(:,1),CA2{k}(:,2));
  end
  if ~isempty(CA2)
    CA2(outcen)          = [];
    CSinfo(outcen,:)     = [];
  end
  % %%%% remove the contours that have <2 points inside (this is the rr we added)
  CA2(CSinfo(:,7)<=2)      = [];
  CSinfo(CSinfo(:,7)<=2,:) = [];




  if ~isempty(CA2) % we only continue if we have contours left

    % here we select the larest contour around the center clusters
    for k=1:length(CSinfo(:,1))
      for kk=1:length(CSinfo(:,1))
        % determine what centre points are within each contour
        outcen2(kk) = inpolygon(CSinfo(k,4),CSinfo(k,5),CA2{kk}(:,1),CA2{kk}(:,2));
      end
      ci      = find(outcen2);                % how many center points are in each contour
      [~, bb] = max(CSinfo(outcen2,3));       % we take the contours that have the largest number of centre points
      cc(k)   = ci(bb);
    end

    %%%% Now we do the final selection, which corresponds to the
    %%%% phase-space holes
    CSinfo_final = CSinfo(unique(cc),:);
    CA_final     = CA2(unique(cc));

  else

    CSinfo_final = [];
    CA_final     = [];


  end


  %%%% If we have found phase-space holes, then we compute some parameters
  if ~isempty(CSinfo_final)

    for ll=1:length(CA_final) % run for each hole contour

      %%%% indices of points in each hole contour
      pinc   = inpolygon(data(:,1),data(:,2),CA_final{1,ll}(:,1),CA_final{1,ll}(:,2));
      %%%% level of the hole contour / edge value
      clevel = mean(interp2(X,Y,Z,(CA_final{1,ll}(:,1)),(CA_final{1,ll}(:,2))));

      %%%% compute the ratio between the edge and the minimum.
      if length(find(pinc))<5
        min_val = double(interp2(X,Y,Z,mean(CA_final{ll}(:,1)),mean(CA_final{ll}(:,2))));
      else
        min_val = min(data(pinc,3));
      end

      hole_stats(ll,:) = [
        mean(double(CA_final{ll}(:,1))),double(mean(CA_final{ll}(:,2))),...   % xy center (based on mean)
        clevel/min_val,...                                                    % level/minima ratio
        range(CA_final{ll}(:,1)),...                                          % width [sec]
        range(CA_final{ll}(:,2))];                                            % height [km/s];
    end
  else

    hole_stats = [];

  end
  close(gcf)

  %%%% Save the data into a structure PSH
  PSH.CN          = numel(CA_final); % number of holes
  PSH.Cxy         = CA_final;
  PSH.Cparam      = double(hole_stats);

  % return filtered dist
  PSH.filt_dist = f1Dn;
end




%% Boundary method

%%%% Boundary-based method %%%%
if strcmp(pshMethod,'boundary') || strcmp(pshMethod,'both')

  % First thing to do is to un-transpose F (funny quib here)
  F = F';

  % then find all minima points (can be lots of 'em)
  % find the minimas in x and y
  % check other directions later
  [idmT,idmV] = findmin_mat(F);

  % next find the minimas that are phase space holes
  [rB,FB,phi,~] = boundary_get_psholes(F,idmT,idmV,36); % 36 gives 10 deg intervals

  % finally select only those PSHs that meet the criteria
  [idht,idhv,xB,yB,nholes] = select_psholes(...
    rB,FB,phi,F,idmT,idmV,BoundaryCriteria);

  PSH.BN = nholes;
  PSH.Bxy = cell(1,nholes);

  for ii = 1:nholes
    % convert to unix time and km/s
    tB = interp1(1:length(f1Dn),f1Dn.time.epochUnix,xB(ii,:));
    vB = interp1(1:length(Vn),Vn,yB(ii,:));

    PSH.Bxy{ii} = [tB',vB'];
  end
  % PSH.Bparam      = [];
  % PSH.CTint       = [];

  % probably there are some other parameters to include such as the
  % ratios you use to find the boundaries

end

PSH.filter_params = filter_params;
PSH.method_params = method_params;


%% make time series objects for the holes

if isfield(PSH,'Cxy')
  PSH.CTS = cell(1,PSH.CN);
  for ii = 1:PSH.CN
    PSH.CTS{ii} = irf.ts_scalar(irf_time(PSH.Cxy{ii}(:,1),...
      'epoch>epochtt'),PSH.Cxy{ii}(:,2));
  end
end

if isfield(PSH,'Bxy')
  PSH.BTS = cell(1,PSH.BN);
  for ii = 1:PSH.BN
    PSH.BTS{ii} = irf.ts_scalar(irf_time(PSH.Bxy{ii}(:,1),...
      'epoch>epochtt'),PSH.Bxy{ii}(:,2));
  end
end

end

%% Filtering auxillary functions

function [f1Dfilt] = filt_mat_gaussian(f1D, wt, wv)
% Internal gaussian filter function

F = f1D.data'; % PSD matrix (not flux!) transposed (?)

[ny,nx]=size(F);
filter_gauss_x = fspecial('gaussian', [1,nx], wt);
filter_gauss_y = fspecial('gaussian', [ny,1], wv);

Ff = zeros(size(F));
for i=1:nx
  Ff(:,i)=filter2(filter_gauss_y, F(:,i));
end

for i=1:ny
  Ff(i,:)=filter2(filter_gauss_x, Ff(i,:));
end

f1Dfilt = f1D;
f1Dfilt.data = Ff';

end


function [f1Dfilt] = filt_mat_runav(f1D, w1, w2)
% Internal running average filter function

F = f1D.data; % PSD matrix (not transposed!)

% simple running average with given widths in x and y

[n1,n2] = size(F);
Ff = zeros(n1,n2);

for ii = 1:n1
  id1 = ii-w1:ii+w1;
  if id1(1)<1; id1(id1<1) = []; end
  if id1(end)>n1; id1(id1>n1) = []; end

  for jj = 1:n2
    id2 = jj-w2:jj+w2;
    if id2(1)<1; id2(id2<1) = []; end
    if id2(end)>n2; id2(id2>n2) = []; end

    Ff(ii,jj) = mean(mean(F(id1,id2)));
  end
end

f1Dfilt = f1D;
f1Dfilt.data = Ff;

end

%% Contour method auxillary functions
function CL = irf_extract_contour_lines(C)
M = C.ContourMatrix;
L = C.LevelList;
tic
cc=0;
for k=1:length(L)
  Lin = find(M(1,:)==L(k));
  Llen = M(2,Lin);
  for kk=1:length(Lin)
    cc=cc+1;
    CLx = M(1,Lin(kk)+1:Lin(kk)+Llen(kk));
    CLy = M(2,Lin(kk)+1:Lin(kk)+Llen(kk));
    CL{cc} = [CLx', CLy'];
  end
end
end

%% Boundary method auxilliary functions

function [idMin1,idMin2] = findmin_mat(F)
%FINDMIN_MAT Summary of this function goes here
%   Detailed explanation goes here

[n1,~] = size(F);

% should preallocate idMin1 and 1dMin2 (TODO)
count = 1;
for ii = 1:n1

  % findpeaks is a great function that finds local maxima
  [~,id1] = findpeaks(-F(ii,:)); % in "x"

  for jj = 1:length(id1)
    [~,id2] = findpeaks(-F(:,id1(jj))); % in "y"

    if ismember(ii,id2) % check if mimima in both x and y
      idMin1(count) = ii;
      idMin2(count) = id1(jj);
      count = count+1;
    end
  end
end

end


function [rB,FB,phi,r] = boundary_get_psholes(F,idm1,idm2,nphi)
%GET_PSHOLES Summary of this function goes here
%   Detailed explanation goes here


% finds boundaries for the minima
% will return a lot of junk that has to be cleaned

% directions to look for peaks
phi = linspace(0,360,nphi+1); % must go all the way around

% radius to look in (unit is index), (can be made anisotropic?)
rlim = 60; % probably safe to use large number
% points along radius (should not matter as long as its big)
nr = 120;

r = linspace(0,rlim,nr);

[n1,n2] = size(F);

% mesh for interpolation
[X1,X2] = meshgrid(1:n1,1:n2);

nmin = length(idm1);

% radius and F-value for the peaks around the minima
rB = zeros(nmin,nphi);
FB = zeros(nmin,nphi);

for ii = 1:nmin

  for jj = 1:length(phi)
    x1 = idm1(ii)+r*cosd(phi(jj));
    x2 = idm2(ii)+r*sind(phi(jj));

    % interpolate
    ft = interp2(X1,X2,F',x1,x2);

    % set the boundary as the first peak (somtimes finds "false" peaks)
    [~,idpeak] = findpeaks(ft);
    if ~isempty(idpeak)
      rB(ii,jj) = r(idpeak(1));
      FB(ii,jj) = interp2(X1,X2,F',x1(idpeak(1)),x2(idpeak(1)));
    end
  end
end
end


function [idh1,idh2,xB,yB,nholes] = select_psholes(rB,FB,phi,F,idm1,idm2,boundaryCriteria)
%SELECT_PSHOLES Summary of this function goes here
%   Detailed explanation goes here


% select holes with given parameters
% holes with zero radius and overlaping holes will always be removed

% x and y of holes
xpeak = idm1'+rB.*cosd(phi);
ypeak = idm2'+rB.*sind(phi);

nmin = length(idm1);

% get rid of holes with zero radius (no input)

% if there are zeros trow
for ii = 1:size(rB,1)
  if ismember(0,rB(ii,:))
    rB(ii,:) = nan;
  end
end


% some simple condition relating to relative depth of minima

% min
for ii = 1:nmin
  if min(FB(ii,:))<boundaryCriteria(1)*F(idm1(ii),idm2(ii))
    rB(ii,:) = nan;
  end
end

% mean
for ii = 1:nmin
  if mean(FB(ii,:))<boundaryCriteria(2)*F(idm1(ii),idm2(ii))
    rB(ii,:) = nan;
  end
end

% max
for ii = 1:nmin
  if max(FB(ii,:))<boundaryCriteria(3)*F(idm1(ii),idm2(ii))
    rB(ii,:) = nan;
  end
end

% some simple condition on absolute height of peak (not great condition)
for ii = 1:nmin
  if min(FB(ii,:))<boundaryCriteria(4)
    rB(ii,:) = nan;
  end
end


% get rid of holes belonging to a minima inside another hole (must be last)
for ii = 1:nmin
  for jj = 1:nmin

    if ii == jj
      continue;
    end

    % magic function
    inpol = inpolygon(idm1(jj),idm2(jj),xpeak(ii,:),ypeak(ii,:));

    if inpol && ~isnan(rB(jj,1)) % conflict between ii and jj
      % check which hole is deeper and keep that
      hole_ii_depth = mean(FB(ii,:))/F(idm1(ii),idm2(ii));
      hole_jj_depth = mean(FB(jj,:))/F(idm1(jj),idm2(jj));
      if hole_ii_depth <= hole_jj_depth
        rB(ii,:) = nan;
      else
        rB(jj,:) = nan;
      end
    end

  end
end

% set output
idSelected = find(~isnan(rB(:,1)));

idh1 = idm1(idSelected);
idh2 = idm2(idSelected);
nholes = length(idSelected);

xB = xpeak(idSelected,:);
yB = ypeak(idSelected,:);

end
