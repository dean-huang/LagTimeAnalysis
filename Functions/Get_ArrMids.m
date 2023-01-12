function [mids] = Get_ArrMids(arr)
% Takes an array (single column or row) and gives you the midpoints.

if min(size(arr)) == 1
    mids = (arr(1:end-1)+arr(2:end))/2;
else
    disp('Please restructure the input so that it is either a single column or single row')
end

end

