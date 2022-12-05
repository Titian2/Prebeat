
function pre_fit_py
%
% function copies fit_ringdowns_degenerate_cmd.ipynb to current
% directory and replaces the definition of the prefix variable in the
% notebook with the path of the current directory 
%
%
% Currently PRE_ALPHA 
% not tested on windows systems 
% 
% Author S.Tait 2022

currentVersion = 0.2; 

pyfilepath = extractBefore(which('prebeat.m'),'prebeat.m');
copyfile(fullfile(pyfilepath,'fit_ringdowns_degenerate_cmd.ipynb'),pwd)
% 
% % % ! cp ~/Dropbox/Python/FitBeats/fit_ringdowns_degenerate_cmd.ipynb . 
% 
% rawText =tex_import(char(fnames('fit_ringdowns*cmd.ipynb')));
% 
% out = regexp(rawText,'prefix = ','Match');
% lines = find(~cellfun(@isempty,out));
% 
% oldpath = regexprep(rawText(lines),{'"prefix = \\"','\\"\\n",'},'');
% 
% rawText(lines) = strrep(rawText(lines),oldpath,pwd);
% 
% fid2 = fopen('fit_ringdowns_degenerate_cmd.ipynb','wt');
% 
% for k =1:numel(rawText)
%     
% fprintf(fid2,'%s',string(rawText{k}));
% end
% fclose(fid2);

cprintf('hyper','\n\nRunning Python Analysis\n\n')

!papermill fit_ringdowns_degenerate_cmd.ipynb output.ipynb

pause(5)

movefile("fit_ringdowns_degenerate_cmd.ipynb","deleteme.ipynb")
movefile("output.ipynb","fit_ringdowns_degenerate_cmd.ipynb")
delete("deleteme.ipynb")

end  