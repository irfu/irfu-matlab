function [output_string]=av_ssub(input_string,num1,num2)
%AV_SSUB change all appearences of '?' and/or '!' in string to some number
% [output_string]=av_ssub(input_string,num) change all appearence of '?' in input_string to num
%  num is converted to string using NUM2STR function
% 
% [output_string]=av_ssub2(input_string,num1,num2)
%  change all appearences of '?' to num1 and '!' to num2
%
% Ex: for ic=1:4,eval(av_ssub('R?=r?;C?=R?.^2;',ic)),end
%     is the same as R1=r1;C1=R1.^2;R2=r2;C2=R2.^2;...

if nargin == 2,
 output_string=strrep(input_string,'?',num2str(num1));
elseif nargin == 3,
 string_a=strrep(input_string,'?',num2str(num1));
 output_string=strrep(string_a,'!',num2str(num2));
end
