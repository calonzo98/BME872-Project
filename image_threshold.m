function img_out = image_threshold(img_in, T)
img_in = double(img_in);
% img_in = cast(img_in, 'uint8');
[row, col] = size(img_in); %Get the input image size
img_out = zeros(size(img_in));
for i=1:row
    for j=1:col
        if img_in(i,j) <T
            img_out(i,j) = 0;
        else 
            img_out(i,j) = img_in(i,j);
        end 
    end 
end 


end 