function currentIntervals=strong_current_search_brst(sc,currentLim,intervalStart,intervalStop)
%MMS.STRONG_CURRENT_SEARCH_BRST searches through avaiable brst data for time intervals where abs(j)>currentLim
%
%MMS.STRONG_CURRENT_SEARCH_BRST searches through brst (l2pre) data from
%intervalStart until Interval stop looking for time intervals with strong currents.
% IMPORTANT: This example loads the tetrahedron quality from /data/mms/ancillary so
% these files are needed if a different harddrive is used than /data/mms.
% It also uses the irfu_index.mat from l2pre dfg to get all brst data time
% intervals available. Currently it will not work because the
% irfu_index.mat is removed. Will be updated once the database of all
% datafiles is finished
%
%   currentIntervals=MMS.STRONG_CURRENT_SEARCH_BRST(sc,currentLim,intervalStart, intervalStop)
%
%
%   INPUT
%     sc    - Looks through the time intervals of brst data from spacecraft
%     sc. Default is sc=1;
%     currentLim    - Looks for current larger than currentLim. Default is
%     zero.
%     intervalStart - gives the time to start the search in the brst data. Default is no specific start.
%     The time must be given in epoch e.g. intervalStart=irf_time([2015 11 01 00 00 00])
%     intervalStop - gives the time to stop the search in the brst data. Default is no specific stop.
%     The time must be given in epoch e.g. intervalStop=irf_time([2015 11 01 00 00 00])
%
%   OUTPUT
%
%   currentIntervals - [tstart tstop] in epochuniX starting from the
%   strongest current to the lowest found.
%
%% Check variables
if nargin ==0
  sc=1;
  currentLim=0;
  Start=false;
  Stop=false;
elseif nargin <2
  currentLim=0;
  Start=false;
  Stop=false;
elseif nargin <3
  Start=false;
  Stop=false;
elseif nargin<4
  Start=true;
  Stop=false;
elseif nargin==4
  Start=true;
  Stop=true;
elseif nargin > 4
  error('Too many input values. See usage: help mms.strong_current_search_brst')
end
%% Go into correct folder and download index
mms.db_init('local_file_db','/data/mms');
switch sc
  case 1
    oldFolder=cd;
    newFolder= '/data/mms/mms1/dfg/';
    cd(newFolder)
    indexMat=load('irfu_index.mat');
    cd(oldFolder)
  case 2
    oldFolder=cd;
    newFolder= '/data/mms/mms2/dfg/';
    cd(newFolder)
    indexMat=load('irfu_index.mat');
    cd(oldFolder)
  case 3
    oldFolder=cd;
    newFolder= '/data/mms/mms3/dfg/';
    cd(newFolder)
    indexMat=load('irfu_index.mat');
    cd(oldFolder)
  case 4
    oldFolder=cd;
    newFolder= '/data/mms/mms4/dfg/';
    cd(newFolder)
    indexMat=load('irfu_index.mat');
    cd(oldFolder)
end
%% Find the time intervals with burst data
%%MMS1
i=1;
% Looks for time intervals with brst data from Nov, 1 2015 and presently
if Start
  if Stop
    for ii=1:length(indexMat.index)
      %strcmp
      if ~isempty(strfind(indexMat.index(ii).filename,'brst/l2pre')) && (EpochTT(indexMat.index(ii).tstart).epochUnix > intervalStart) && (EpochTT(indexMat.index(ii).tstart).epochUnix < intervalStop)
        timeIntervalMMS(i,:)=[indexMat.index(ii).tstart indexMat.index(ii).tstop];
        i=i+1;
      end
    end
  else
    for ii=1:length(indexMat.index)
      %strcmp
      if ~isempty(strfind(indexMat.index(ii).filename,'brst/l2pre')) && (EpochTT(indexMat.index(ii).tstart).epochUnix > intervalStart)
        timeIntervalMMS(i,:)=[indexMat.index(ii).tstart indexMat.index(ii).tstop];
        i=i+1;
      end
    end
  end
else
  for ii=1:length(indexMat.index)
    %strcmp
    if ~isempty(strfind(indexMat.index(ii).filename,'brst/l2pre'))
      timeIntervalMMS(i,:)=[indexMat.index(ii).tstart indexMat.index(ii).tstop];
      i=i+1;
    end
  end
end
test=1;
amountStrongCurrent=0;
for int=1:length(timeIntervalMMS(:,1))
  try
    Tint=irf.tint(timeIntervalMMS(int,1), timeIntervalMMS(int,2));
    disp(['This is interval nr:',num2str(int), ' and ', num2str(amountStrongCurrent), ' intervals with strong current have been found'])
    % Magnetic Field
    disp('Loading Magnetic fields');
    c_eval('B?=mms.db_get_ts(''mms?_dfg_brst_l2pre'',''mms?_dfg_brst_l2pre_gse'',Tint);',1:4);
    %Resamples according to the time line where all spacecraft has data
    timeStart=max([B1.time.start.epochUnix B2.time.start.epochUnix B3.time.start.epochUnix B4.time.start.epochUnix],[],2);
    timeStart=EpochUnix(timeStart);
    timeStop=min([B1.time.stop.epochUnix B2.time.stop.epochUnix B3.time.stop.epochUnix B4.time.stop.epochUnix],[],2);
    timeStop=EpochUnix(timeStop);
    timeLog=B1.time >= timeStart & B1.time <=timeStop;
    newTime=B1.time(timeLog,:); %#ok<NASGU>
    c_eval('B? = B?.resample(newTime);',1:4);
    % Spacecraft Position
    disp('Loading Spacecraft Position');
    R  = mms.get_data('R_gse',Tint); %Cailbrated position
    if   length(R.gseR1(1,:))==4 && length(R.gseR2(1,:))==4 && length(R.gseR3(1,:))==4 && length(R.gseR4(1,:))==4
      % Assume the first column is time
      c_eval('R? =irf.ts_vec_xyz(R.time, R.gseR?(:,2:4));',1:4);
      clear R
      c_eval('R? = R?.resample(B1);',1:4);
    else
      c_eval('R? =irf.ts_vec_xyz(R.time, R.gseR?);',1:4);
      clear R
      c_eval('R? = R?.resample(B1);',1:4);
    end

    %Quality data comes 2 days late
    % Load quality of tetrahedron
    disp('Checks tetrahedron quality');
    quality=mms.db_get_variable('mms_ancillary_defq','quality',Tint);
    if isempty(quality.quality)
      disp('No tetrahedron quality available right now');
      continue
    else
      quality=irf.ts_scalar(EpochTT(quality.time),quality.quality);
      quality=quality.resample(B1);
      tetrahedronGood= quality.data > 0.7; %#ok<NASGU>
      % Removes all time steps with bad tetrahedron quality
      c_eval('R.C? = R?(tetrahedronGood);',1:4);
      c_eval('B.C? = B?(tetrahedronGood);',1:4);
    end

    if isempty(B1) || isempty(B2) || isempty(B3) || isempty(B4)
      error('Tetrahedron quality is not good enough to use curlometer method');
    else
      disp('Looking for strong currents in brst data');
      curlB = c_4_grad(R,B,'curl');
      j     = curlB.data./1.0e3.*1e-9./(4*pi*1e-7); %A/m^2 if B in nT and R in km
      jabs=irf_abs([curlB.time.epochUnix j],1);
      strongCurrents=jabs>currentLim;
      if sum(strongCurrents) >= 1
        IntervalStrongCurrent=[curlB.time(strongCurrents).epochUnix jabs(strongCurrents,1)];
      else
        disp('No strong currents in this interval')
        continue
      end
    end

    if isempty(IntervalStrongCurrent)
      continue
    else
      if test==1

        currentIntervals=IntervalStrongCurrent;
        amountStrongCurrent=amountStrongCurrent+1;
      else
        disp('Saves the data in structure');

        currentIntervals(length(currentIntervals(:,1))+1:length(IntervalStrongCurrent(:,1))+length(currentIntervals(:,1)),:)=IntervalStrongCurrent;

        amountStrongCurrent=amountStrongCurrent+1;

      end
      test=test+1;

    end
  catch
    warning('Problem using function.  Move to next interval');
    continue
  end
end


%Sorts the currents found from highest to lowest size
[~,IndexSize] = sort(currentIntervals(:,2),'descend');
currentIntervals=currentIntervals(IndexSize,:); %[currentTime currentAbs]

%TODO - still some overlapping between the intervals
jTemp = [currentIntervals ones(size(currentIntervals,1),1)];
for ii = 1:size(jTemp,1)
  if jTemp(ii,3)
    jTemp(jTemp(:,1)<jTemp(ii,1)+30 & jTemp(:,1)>jTemp(ii,1)-30,3)=0;
    jTemp(ii,3)=1;
  end
end

currentIntervals=currentIntervals(logical(jTemp(:,3)),:);
currentIntervals=[currentIntervals(:,1)-30 currentIntervals(:,1)+30 currentIntervals(:,2)];
end
