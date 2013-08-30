%
% Ceftoolbox 0.4 beta 2 release, 2007-04-24
% By Josef Hook, removethis-jhook (a t) rssd.esa.int
%
%  
%  General:
%                        cefRead              Main read function
%                        cefWrite             Main write function
%
% Low level functions:  
%                        cefReadData          Low level read function
%                        cefReadMetaData      Low level read function
%                        cefWriteMetaData     Low level write function   
%                        cefWriteData         Low level write function   
%  Date conversion:
%                        cefTimeToMjs
%                        cefTimeToDatenum
%                        cefMjsToTime
%                        cefTimeSpanToMjs
%                        cefMjd2kToDatenum
%                        cefMjsToDatenum
%                        cefDatenumToMjs
%                        cefMjsToMjd2k
%                        cefMjd2kToMjs
%                        cefTimeToDatenum
%                        cefTimeTojd2k
%                        cefMjsToDate           Split MJS to date parts
%
%
%    CEF options: 
%                        cefGet                 Return default options
%
%    Helper functions:
%                        cefPlot
%                        cefTransform           General transformation
%                                               function for transformations 
%                                               from different coordinatesystems 
%                                               like GSE, ECI 
%                        cefPackageData         Package raw data into a
%                                               format used by cefWriteData()
%                        cefMeanInterp          Mean window interpolator
%                        cefUnion               Union for logicals
%                        cefSymbolColormap
%                        cefPlotTimeline 
%                        cefSplit               Perl like split function
%                        cefChunk               Extract a chunk of cefdata
%                                               given start and end date.
%                        cefRemoveNan           Remove/replace NAN
%                                               numbers
%                
%                        cefRemoveFillval       Removes/replaces fillvalues
%                        cefVar2cdfVar          Convert CEF variable metadata to cdf formated metadata.  
%     Small applications:
%                        test
%
%

