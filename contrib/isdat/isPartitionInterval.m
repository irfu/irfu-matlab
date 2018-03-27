function [partIntrvlStart, partIntrvlDur]=isPartitionInterval(intrvlStart,intrvldur,partLen,cntLstStart,cntLstDur)
% isPartitionInterval - Partition time intervals.
%
%   [PS,PD]=isPartitionInterval(S,D,L) partitions the time interval, given by
%   start time S and duration D, into intervals of length L. The returned
%   intervals are given as an array of start times PS and interval durations
%   PD. All interval start times are given as numerical row vectors in the
%   format [yyyy mm dd HH MM SS.s] and all interval durations or lengths
%   are given in seconds.
%
%   [PS,PD]=isPartitionInterval(S,D,L,CS,CD) partitions the time interval
%   into portions of maximally D in length which coincides with the
%   additional intervals specified by CS, DS. In other words the function
%   will return data intervals of maximal length L which are both within
%   the main interval S,D AND in CS,CD. This is useful in conjunction with
%   the function isGetContentLite() which returns available dataset
%   time intervals.
%
%   Example:
%     >> [CS,CD]=isGetContentLite(DB,'Test','*','wave','sine','wf');
%     >> [PS,PD]=isPartitionInterval([1993 01 01 00 00 00],60*60,60,CS,CD);

if nargin == 3
  cntLstStart=intrvlStart;
  cntLstDur=intrvldur;
end
cntStartS=toepoch(cntLstStart);
cntStopS=cntStartS+cntLstDur;
intrvlStartS=toepoch(intrvlStart);
intrvlStopS=intrvlStartS+intrvldur;

%Discard all times which are outside of the requested interval
indBeforeStart=find(cntStartS < intrvlStartS);
if ~isempty(indBeforeStart) 
  if cntStopS(indBeforeStart(end)) > intrvlStartS
    cntStartS(indBeforeStart(end))=intrvlStartS;
    indBeforeStart(end)=[];
  end
  cntStartS(indBeforeStart)=[];
  cntStopS (indBeforeStart)=[];
end
indAfterStop=find(cntStopS > intrvlStopS);
if ~isempty(indAfterStop) 
  if cntStartS(indAfterStop(1)) < intrvlStopS
    cntStopS(indAfterStop(1))=intrvlStopS;
    indAfterStop(1)=[];
  end
  cntStartS(indAfterStop)=[];
  cntStopS (indAfterStop)=[];
end
cntDurS=cntStopS-cntStartS;
partIntrvlStartS = [];
partIntrvlDur = [];
for segInd=1:length(cntStartS)
  partSeg=cntStartS(segInd):partLen:cntStopS(segInd);
  partIntrvlStartS = [partIntrvlStartS; partSeg'];
  partIntrvlDur = [partIntrvlDur; partLen*ones(length(partSeg),1)];
  partIntrvlDur(end)=mod(cntDurS(segInd),partLen);
  if partIntrvlDur(end)==0.0
     partIntrvlDur(end)=[];
  end
end
partIntrvlStart=fromepoch(partIntrvlStartS');
