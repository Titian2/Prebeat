function [data] =  dynamic(inarr,outlen,padding) 

%
% [data] =  dynamic(inarr,outlen,padding) 
%
% This function takes a set of data and combines them with a set length of
% NaNs to allow storage of multiple datasets of different lengths. The
% length of data is concatinated with NaN's to reach the desired length 
%
% used for making 'dynamic' arrays
%
% INPUTS
%
% inarr   = input array of data
% outlen  = desired length of output array 
% padding = type of padding value e.g. 'NaN', 'zeros' , 'ones' 
% 
% OUTPUTS 
%  
% data   = concatinated array with nans of length outlen 
%
% Author S.Tait -2019 
% 

if  size(inarr,2)>1
    if strcmpi('NaN',padding)
        data = [inarr; NaN(outlen-length(inarr(:,1)),size(inarr,2))];
    elseif strcmpi('ones', padding)
        data = [ inarr  ;ones(outlen-length(inarr),size(inarr,2))];
    elseif strcmpi('zeros',padding)
        data = [ inarr  ;zeros(outlen-length(inarr),size(inarr,2))];
    elseif strcmpi('NaT',padding)
        data = [ inarr  ;NaT(outlen-length(inarr),size(inarr,2))];
    end
else
    if strcmpi('NaN',padding)
        data = [ inarr  ;NaN(outlen-length(inarr),1)];
    elseif strcmpi('ones', padding)
        data = [ inarr  ;ones(outlen-length(inarr),1)];
    elseif strcmpi('zeros',padding)
        data = [ inarr  ;zeros(outlen-length(inarr),1)];
    elseif strcmpi('NaT',padding)
         data = [ inarr  ;NaT(outlen-length(inarr),1)];
    end
end
end 