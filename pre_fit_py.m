
function pre_fit_py
%
% function copies fit_ringdowns_degenerate_cmd.ipynb to current
% directory and replaces the definition of the prefix variable in the
% notebook with the path of the current directory 
%
% Author S.Tait 2022

if ismac 
! cp ~/Dropbox/Python/FitBeats/fit_ringdowns_degenerate_cmd.ipynb . 
elseif ispc 
    pyfilepath = extractBefore(which('prebeat.m'),'prebeat.m');
    copyfile(fullfile(pyfilepath,'fit_ringdowns_degenerate_cmd.ipynb'),pwd)
end 

rawText =tex_import(char(fnames('fit_ringdowns*cmd.ipynb')));

out = regexp(rawText,'prefix = ','Match');
lines = find(~cellfun(@isempty,out));


if ismac

oldpath = regexprep(rawText(lines),{'"prefix = \\"','\\"\\n",'},'');

rawText(lines) = strrep(rawText(lines),oldpath,pwd);

end 

   
if ispc
    str = rawText(lines);
    idx = regexp(str,'"') ;
    t =char(str);
    oldpath = t(idx{:}(2):idx{:}(3));

    rawText(lines) = strrep(rawText(lines),oldpath,strcat('r',pwd,'"'));
end 


if ~isempty(find(contains(rawText,'makeplot=False')))
    rawText(find(contains(rawText,'makeplot=False'))) = strrep(rawText(find(contains(rawText,'makeplot=False'))),'makeplot=False','makeplot=True');
end 

fid2 = fopen('fit_ringdowns_degenerate_cmd.ipynb','wt');

for k =1:numel(rawText)
    
fprintf(fid2,'%s',string(rawText{k}));
end
fclose(fid2);

cprintf('hyper','\n\nRunning Python Analysis\n\n')

!papermill fit_ringdowns_degenerate_cmd.ipynb output.ipynb
pause(5)
if ismac
! mv ./fit_ringdowns_degenerate_cmd.ipynb ./deleteme.ipynb
! mv ./output.ipynb ./fit_ringdowns_degenerate_cmd.ipynb
! rm deleteme.ipynb
end 

if ispc

movefile("fit_ringdowns_degenerate_cmd.ipynb","deleteme.ipynb")
movefile("output.ipynb","fit_ringdowns_degenerate_cmd.ipynb")
delete("deleteme.ipynb")

end 
end  