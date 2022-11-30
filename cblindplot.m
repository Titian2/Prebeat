function [cb] = cblindplot(input)

cb = [0,0,1;1,0,0;0.660156250000000,0.660156250000000,0.660156250000000;0,0,0;1,0.644531250000000,0;1,0,1;0,0.500000000000000,0.500000000000000;0,0,0.542968750000000;0,0.390625000000000,0;0,1,1;0.597656250000000,0.195312500000000,0.796875000000000;1,1,0];

if exist('input','var') 
    for i =1:numel(input)
        try
        tmp(i,:) = cb(input(i),:);
        catch 
            cprintf('err','\nInput value exceeeds array dimensions. Please specify values between 1-9\n')
            return 
        end 
    end 
    cb = tmp; 
end 
    
    
    
    

end 