%The function CCIQ = contrast_estimation(img1, img2) takes two images and
%retruns a value corresponding to the quality of the contrast between the
%two images. A higher value relates to a higher contrast
%
%This metric includes a local index named edge-based contrast criterion
%(ECC) and three global measures; entropy, correlation coeff and mean
%intensity. An optimal linear combination of these quantities are achieved
%using the particle swam optimization
%
%The ECC is the first term which is based on the fact that human perception
%system is very sensitive to edges. If the original image and the modified 
%one are illustrated as X and Y respectively, when ECC(Y) is greater than 
%ECC(X), it could be concluded that the contrast is improved. This measure
%is normalized, consequently a value of ECCN (X,Y) greater than unity 
%indicates the contrast improvement of Y 
%
%Entropy is the second term: , a high contrast image has a large entropy 
%Entropy + ECC can approximate image contrast, but can;t evaluate the
%quality in case of over enhanced images. In this situation, the structure
%of the original image is not preserved, it is necessary to employ some
%other quantities that could measure the structural similarity between the
%original image and enhanced image -- correlation coef and mean intensity
%
%Correlation coef computes the degree of linear correlation b/t two images
%
%mean intesntiy - absolute mean brightness error (AMBE) is a measure to
%describe brightness preservation of the contras enhanced image and shows
%the diff between the two images mean intensity. 
%
%this method was vhosen due to the superiod perfromance compared to other
%methods such as RMSE, SROCC and PLCC
%

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

%Step 5 Linear Combination 
% a1 = optimvar('a1');
% a2 = optimvar('a2');
% a3 = optimvar('a3');
% a4 = optimvar('a4');
% prob = optimproblem;
% prob.Objective = a1*ECC_XY + a2*H_XY + a3*R_XY + a4*AMBE_XY;
% prob.Constraints.cons1 = a1+a2+a3+a4==1;
% prob.Constraints.cons2 = a1>=0;
% prob.Constraints.cons3 = a1<=1;
% prob.Constraints.cons4 = a2>=0;
% prob.Constraints.cons5 = a2<=1;
% prob.Constraints.cons6 = a3>=0;
% prob.Constraints.cons7 = a3<=1;
% prob.Constraints.cons8 = a4>=0;
% prob.Constraints.cons9 = a4<=1;
% 
% sol = solve(prob, 'Solver', 'intlinprog');

% CCIQ = (sol.a1*ECC_XY) + (sol.a2*H_XY) + (sol.a3*R_XY) + (sol.a4*AMBE_XY);

a1 = 0.6655;
a2 = 0.3174;
a3 = 0.2183;
a4 = 0.4448;

CCIQ = (a1*ECC_XY) + (a2*H_XY) + (a3*R_XY) + (a4*AMBE_XY);

end 


