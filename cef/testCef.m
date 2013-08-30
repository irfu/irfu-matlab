% Test script for testing reading and writing CEF files
%
%
clear all
close all

file='./testfiles/C4_CP_EFW_L1_P12__20010706_060000_064341_V01.cef'
header_efw='./testfiles/efw'
header_fgm='./testfiles/fgm'
header_glob='./testfiles'
file2='./testfiles/C4_CP_FGM_SPIN__20010706_211607_20010709_062406_V01.cef'




%[header, var]=cefReadMetaData(file2,header_fgm, header_glob);
%data=cefReadData(file2,[1 inf],[],[],header_fgm, header_glob);

[data,var,header]=cefRead(file2, [],[], {header_fgm, header_glob});

%[header2, var2]=cefReadMetaData(file,header_efw, header_glob);
%data2=cefReadData(file,[1 inf],[],[],header_efw, header_glob);

[data2, var2, header2]=cefRead(file,[],[], {header_efw,header_glob});

% example of plotting
%plot(cell2mat(data.B_mag__C4_CP_FGM_SPIN));


% example of writing
cefWrite('./testfiles/fil',data, var, header);  
%cefWriteMetaData(header,var,[],'./testfiles/fil','w')
%cefWriteData(data,var,'./testfiles/fil')

include_tags={'CL_CH_MISSION.ceh','C1_CH_OBS.ceh','CL_CH_EFW_EXP.ceh'};

cefWrite('./testfiles/fil2',data2, var2, header2, include_tags);
%cefWriteMetaData(header2,var2,[],'./testfiles/fil2','w')
%cefWriteData(data2,var2,'./testfiles/fil2')


%plot the data 
cefPlot(data);
cefPlot(data2);

