function [sorted_files,sorted_folders] = fnames(pattern,directory)
%
% [sorted_files,sorted_folders] = fnames(pattern,directory)
%
% Function generates sorted cell array of files matching input pattern in
% current directory 
%
%
% INPUTS 
% 
% pattern       : Filename pattern, uses same input as 'dir' 
%           see dir 
% directory     : Contains path which can be used to call file names from
%                 different directory. Default is current path 
% 
% OUTPUTS 
% fnames        : Cell array of file names matching pattern. Will be sorted
%               if numbers are present in the name
%
% Requirements 
% natsort.m    
%       see natsort.m
%
% author S.tait 2020 
here = pwd; 

if exist('directory','var') && ~isempty(directory)
    cd(directory) 
end 
    
    
files = dir(pattern); 
folders = {files.folder};
files= {files.name}; 
folders = [files' folders'];
sorted_files= natsort(files)'; 
for i =1:length(sorted_files)
    [~,b] = ismember(sorted_files(1),folders(:,1));
    sorted_folders(i,:) = folders(b,2);
    clear b 
end 
cd(here)
    


