function  csvout  = csvstring(inarray,seperate,buffer)
%     csvout  = csvstring(inarray,seperate,buffer)
%
% Function takes in a cell array and outputs csv string of the input cell
% array
%
%
% Author S.Tait 2021

for i=1:numel(inarray)
    if i ==1
        if exist('seperate','var')
            if ~exist('buffer','var')
            csvout = strcat("'",inarray{i},"'");
            else 
            csvout = strcat(string(buffer),inarray{i},string(buffer));
            end 
        else
            csvout = inarray{i};
        end
    else
        if exist('seperate','var')
            csvout  = strcat(csvout,',',"'", inarray{i},"'");
        else
            
            csvout = strcat(csvout,',',inarray{i});
        end
    end
    
end



end