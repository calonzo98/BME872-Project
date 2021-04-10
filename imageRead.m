function [img, info] = imageRead(path,imageformat)

if strcmp(imageformat, '.mhd') == 1
    [img, info] = read_mhd(path);
elseif strcmp(imageformat, '.raw') == 1
    [img, info] = read_mhd(path);
elseif strcmp(imageformat, '.dcm') == 1
   img = dicomread(path);
   info = dicominfo(path);
elseif strcmp(imageformat, '.png') == 1
    img = imread(path);
    info = 0;
elseif strcmp(imageformat, '.pgm') ==1
    img = imread(path);
    info = 0;
else 
    disp("Invalid File Type")
    
end 

end 