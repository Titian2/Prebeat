
function pre_fit_py
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

notebookFile = 'fit_ringdowns_degenerate_cmd.ipynb';

modify_ipynb(fullfile(pwd,notebookFile),pwd);


cprintf('hyper','\n\nRunning Python Analysis\n\n')

!papermill fit_ringdowns_degenerate_cmd_modified.ipynb output.ipynb
pause(5)

movefile("fit_ringdowns_degenerate_cmd.ipynb","deleteme.ipynb")
movefile("fit_ringdowns_degenerate_cmd_modified.ipynb","deleteme2.ipynb")
movefile("output.ipynb","fit_ringdowns_degenerate_cmd.ipynb")
delete("deleteme.ipynb")
delete("deleteme2.ipynb")


function modify_ipynb(file_path, new_prefix)
    % Read the file content as a JSON
    fid = fopen(file_path, 'r', 'n', 'UTF-8');
    raw = fread(fid, inf, 'uint8=>char')';
    fclose(fid);
    data = jsondecode(raw);
    
    new_prefix = strcat('"', new_prefix, '"');
    if ispc 
        new_prefix = strcat('r',new_prefix);
    end

    % Modify content
    for i = 1:size(data.cells,1)
        cell = data.cells(i);
        if strcmp(cell.cell_type, 'code')
            for j = 1:length(cell.source)
                line = cell.source{j};
                if startsWith(strtrim(line), 'prefix =')
                    % Split the line at '=' and only modify the part after it
                    data.cells(i).source(j) = join({'prefix =', new_prefix,newline}); 
                    break; % Break out if the line is found and modified
                end
            end
        end
    end

    % Write the modified content back to file
    modified_content = jsonencode(data,'PrettyPrint',true);
    modified_path = strrep(file_path,'.ipynb','_modified.ipynb');

    fid = fopen(modified_path, 'w', 'n', 'UTF-8');
    fprintf(fid,'%s', modified_content);
    fclose(fid);
end




end