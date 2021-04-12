function out_img = apply_mask(img, img_mask, type)

if (strcmp(type, 'ImageType')) == 1
    I = img.data;
    Mask = img_mask.data;
    out_img = I.*Mask;
    
elseif (strcmp(type, 'double')) == 1 
    I = img;
    Mask = img_mask;
    out_img = I.*Mask;
end