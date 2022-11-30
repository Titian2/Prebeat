
function pre_fit_py
%
% function copies fit_ringdowns_degenerate_cmd.ipynb to current
% directory and replaces the definition of the prefix variable in the
% notebook with the path of the current directory 
%
% Author S.Tait 2021 


! cp ~/Dropbox/Python/FitBeats/fit_ringdowns_degenerate_cmd.ipynb . 

rawText =tex_import(char(fnames('fit_ringdowns*cmd.ipynb')));

out = regexp(rawText,'prefix = ','Match');
lines = find(~cellfun(@isempty,out));

oldpath = regexprep(rawText(lines),{'"prefix = \\"','\\"\\n",'},'');

rawText(lines) = strrep(rawText(lines),oldpath,pwd);

fid2 = fopen('fit_ringdowns_degenerate_cmd.ipynb','wt');

for k =1:numel(rawText)
    
fprintf(fid2,'%s',string(rawText{k}));
end
fclose(fid2);

cprintf('hyper','\n\nRunning Python Analysis\n\n')

!papermill fit_ringdowns_degenerate_cmd.ipynb output.ipynb
pause(5)
! mv ./fit_ringdowns_degenerate_cmd.ipynb ./deleteme.ipynb
! mv ./output.ipynb ./fit_ringdowns_degenerate_cmd.ipynb
! rm deleteme.ipynb
end  