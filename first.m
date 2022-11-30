function a= first(array,plus)
%
% takes in an array returning only the first element
%
%
% Author  S.Tait 2019 

try 
    if isa(array  ,'double')
        
        if exist('plus','var')
            if plus> 0
                try 
                a = array(1,plus);
                catch 
                    a = array(plus);
                end 
            else
                cprintf('err','ERR:\tOffset value must be greater than 0\n')
            end
        else
            a = array(1,1);
        end
    else
        if exist('plus','var')
            if plus> 0
                try 
                a = array(1,plus);
                catch 
                    a = array(plus);
                end 
            else
                cprintf('err','ERR:\tOffset value must be greater than 0\n')
            end
        else
            a = array(1);
        end
    end
catch 
     cprintf('err', 'ERR:\t Array Must be larger than (1x1)\n')
     a = NaN; 
end 




end 