function rawText = tex_import(fileName) 
%
%
%
% Function imports a text file and outputs as a cell array 
%
%
%
% Author S.Tait 2021
%

fid = fopen(fileName);
try
    % load each line as a single string
    tmp = textscan(fid, '%s', 'delimiter', '\n');
    fclose(fid);
catch
    cprintf('err', '\n ERR: No File with specified name could be found')
    cprintf('err', '\n Please check it has been entered correctly\n')
    return
end
% textscan wraps its results in a cell, remove that wrapping
rawText = tmp{1};

end 
% % check for instances where more than one space is used
% out = regexp(rawText,'\s{2,}','Match');
% lines = find(~cellfun(@isempty,out));