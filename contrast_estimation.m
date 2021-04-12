function CCIQ = contrast_estimation(img1, img2)

%Step 1 Find local measure names edge-based criterion (ECC)

ECC_CT1 = ECC(img1);
ECC_CT2 = ECC(img2);

ECC_XY = 2/(1+(1-ECC_CT2)/(1-ECC_CT1));

img1_u = cast(img1, 'uint8');
img2_u = cast(img2, 'uint8');
% Step 2 Find Entropy 
entropy1 = entropy(img1_u);
entropy2 = entropy(img2_u);
H_XY = 2/(1+(entropy1/entropy2));

% Step 3 Correlation coeff

imgcorr = corrcoef(img1(:), img2(:));

R_XY = imgcorr(1,2);

% Step 4 Ambsolute Mean Brightness Error

AMBE_XY = (2.*(mean(img1(:))).*(mean(img2(:))))/((mean(img1(:)).^2)+(mean(img2(:)).^2));

%% Step 5 Linear Combination 

CCIQ = (ECC_XY+H_XY+R_XY+AMBE_XY)/4;


end 


