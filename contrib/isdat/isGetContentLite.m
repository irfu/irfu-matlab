function [start,dur]=isGetContentLite(db,proj,mem,inst,sig,sen,chan,param)
% isGetContentLite - Get list of dataset contents on an ISDAT server.
%
%   [START,DUR]=isGetContentLite(DB,PROJ,MEM,INST,SIG,SEN,CHAN,PARAM)
%   This ISDAT function returns an array of intervals for which data
%   files exist for the specificied dataset on database server. DB is
%   an open database connection handle, typically opened by 
%   Mat_DbOpen(). The stings PROJ, MEM, INST, SIG, SEN, CHAN, PARAM
%   specify the dataset. The arrays START and DUR are the start time
%   and duration of the intervals for which data may exist. The start
%   time is a numerical row vector in the format [yyyy mm dd HH MM SS]
%   and the duration is given in seconds.
%
%   Example:
%     [s,d]=isGetContentLite(DB,'Test', '*', 'wave', 'sine', 'wf', ' ', ' ')
%       returns s=[1993 8 1 15 10 0; 1993 8 2 18 30 0] and d=[100; 100].
%
%   See also: Mat_DbOpen

[dataSetId, err]=Mat_DbName2Spec(db,proj,mem,inst,sig,sen,chan,param);

% Prepare a content request for the Test/*/wave/sine dataset.
desc=zeros(9,1);
desc(1:7) = dataSetId;
desc(8) = -1;
desc(9) = -1;
[periods,spec_ret,message,dcont,err] = Mat_DbGetContent(db,desc);

%Sort them into ascending date order (Isdat doesn't garantee this)
[ignore,sortIndex]=sort(periods(1,:),2);
periods=periods(:,sortIndex);

nsINs=1e9;
start=fromepoch(periods(1,:)+periods(2,:)/nsINs);
dur=(periods(3,:)+periods(4,:)/nsINs)'; 
