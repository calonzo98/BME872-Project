function [threshold] = auto_thresholding(gm, factor)

[row, col] = size(gm);
var = stdfilt(gm);

T = zeros(1, 'like', var);
for i = 1:row
    for j = 1:col
            T = T + var(i,j);
    end
end

T = T/numel(var); %global mean

threshold = T*factor;

end 
