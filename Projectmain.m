%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BME872: Biomedical Image Analysis
% Projecct: Automated Image Quality Assessment in Mdedical Images  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Names: Claudia Alonzo & Alexandra Zsivanov
% Student IDs: 500745327 & 500750592
%

%% Housekeeping
close all
clear all

%% Noise Quality in Images 
% Assuming that noise is additive, stationary and has zero mean --> white
% noise 
% The goal is to estimate the standard deviation 
% Average method is the most robust 


% Test out no-reference noisiness metric 4

%% Initialize CT Images 
[training_post_half, infoCThalf] = imageRead('C:\Users\Claudia\Documents\MATLAB\BME872\Lab1\Lab1-LungCT\noise_0.5x_post.mhd', '.mhd');
noiseCT_half = training_post_half.data;
noiseCT_half=noiseCT_half(:,:,143);
[training_post_1, infoCT1] = imageRead('C:\Users\Claudia\Documents\MATLAB\BME872\Lab1\Lab1-LungCT\training_post.mhd', '.mhd');
noiseCT1 = training_post_1.data;
noiseCT1=double(noiseCT1(:,:,143));
[training_post_10, infoCT10] = imageRead('C:\Users\Claudia\Documents\MATLAB\BME872\Lab1\Lab1-LungCT\noise_10x_post.mhd', '.mhd');
noiseCT10 = training_post_10.data;
noiseCT10=noiseCT10(:,:,143);
[mask, info_mask] = imageRead('C:\Users\Claudia\Documents\MATLAB\BME872\Lab1\Lab1-LungCT\training_mask.mhd', '.mhd');

%% sample plotting 
% figure;
% imshow(noiseCT_half_slice, [])
% imdisplayrange;
% colorbar6
% title('CT Image with .5x noise ');
% test = cast(noiseCT_half, 'uint8');
% figure;
% [bins, freq]=intensityHistogram(noiseCT_half, 5,1);

%% Initializing Brain MRI2 Images 
% Initializing mammogram files 
brainMRI1 = load('brainMRI_1.mat'); brainMRI1 = brainMRI1.vol;brainMRI1 = brainMRI1(:,:,90);
brainMRI2 = load('brainMRI_2.mat'); brainMRI2 = brainMRI2.vol;brainMRI2 = brainMRI2(:,:,90);
brainMRI3 = load('brainMRI_3.mat'); brainMRI3 = brainMRI3.vol;brainMRI3 = brainMRI3(:,:,90);
brainMRI4 = load('brainMRI_4.mat'); brainMRI4 = brainMRI4.vol;brainMRI4 = brainMRI4(:,:,90);
brainMRI5 = load('brainMRI_5.mat'); brainMRI5 = brainMRI5.vol;brainMRI5 = brainMRI5(:,:,90);
brainMRI6 = load('brainMRI_6.mat'); brainMRI6 = brainMRI6.vol;brainMRI6 = brainMRI6(:,:,90);

%% initialize image size 

[H, W] = size(noiseCT1);

%% Step 1 -> edge detection 
    
% [BW, threshout,GV, GH] = edge(noiseCT1, 'sobel');
% G = abs(GV)+abs(GH);
% figure;
% imshow(BW, []);
% figure;
% imshow(G, []);
% figure;
% intensityHistogram(G, 500, 1);

%% Step 1 edge detection different way -> using fspecial 
h = fspecial('sobel');
v = rot90(fspecial('sobel'));
gh = imfilter(noise3, h);
gv = imfilter(noise3, v);
[G, gdir] = imgradient(gh,gv);
figure;
subplot(2,2,1)
imshow(G, []);
subplot(2,2,3)
imshow(gh, [])
title ('horizontal');
subplot(2,2,4)
imshow(gv, []);
title('verticle');
figure;
imhist(brainMRI1);

%% Step 2: Threshold for edge map test1
% [counts, bins] = imhist(brainMRI1);
% cdf = cumsum(counts) / sum(counts);
% plot(cdf);

%% Step 2: Threshold for edge map test1

%% Ploting CDF  -- works like a charm!!!
[hist, normhist] = histo_norm(cast(G, 'uint8'));
cdf1 = hist_cum(normhist);
figure;
plot(cdf1);
%% creating a loop to find intensity x when cdf y=0.1 (i.e 10% of entire image intensities)
[r,c] = size(cdf1);
for i = 1:r
    for j = 1:c
        if (cdf1(i,j)>=0.1)&&(cdf1(i,j)<=0.105)
            T = j;
        end
    end
end

%% Try autothreshold 
Ttest = auto_thresholding(G, 1);

%% otsu method 
t = graythresh(noise3);
bw = imbinarize(noise3,t);
imshow(bw);
%% find edge map 
Gth = noiseCT1>t;
Gth = double(Gth);
imshow(Gth,[]);

%%
G_edge = edge(G,'sobel', Ttest);
subplot(1,2,1)
imshow(G_edge, []);
subplot(1,2,2)
imshow(G, []);

%% Step 2 Threshold for edge map 
% figure;
% imhist(G);
% figure;
% histogram(G,'Normalization','cdf');
% x = prctile(G, 0.1);
% y = quantile(G,0.1);
% imshow(x, []);

%otsu -- thresholding 

%% Excluding edges found in edge map from original image 
% %  -- this doesnt make sense because we can't subtract this binary image
% from the original greyscale image -- does not match 
% img = G-Gth;
% subplot(1,2,1)
% imshow(im2double(G), []);
% subplot(1,2,2)
% imshow(img, []);
%% Step 3 Laplacian 
% multiply Gth with laplacian to suppress image stuctures 
kernel = [1 -2 1; -2 4 -1; 1 -2 1];
img_lap = double(imfilter(G_edge, kernel));
imshow(img_lap, []);

% kernel = [1 -2 1; -2 4 -1; 1 -2 1];
% Gth_lap = conv2(Gth, kernel);
% imshow(Gth_lap, []);

%% Step 4 Calculating standard deviation
std_test6 = sqrt(pi/2)*(1/(6*(H-2)*(W-2))).*sum(abs(img_lap), [1 2]);

%std_noiseCT1 = sqrt(pi/2)*(1/(6*(H-2)*(W-2))).*mean(abs(img_lap(:)));

%% Plotting standard deviation 



%% Performance in estimating the standard deviation 
% add noise to CT1 and find estimation ratio
std_CT1 = noise_estimation(noiseCT1); % assume this to be base image (i.e no added noise)
std_CT10 = noise_estimation(noiseCT10); % assume this to be base image (i.e no added noise)
std_CThalf = noise_estimation(noiseCT_half);
noise1= noiseCT1 + (1* randn(512, 512)); % adding noise with std 1
noise2= noiseCT1 + (5* randn(512, 512)); % adding noise with std 2
noise3= noiseCT1 + (10* randn(512, 512)); % adding noise with std 3
q1 = (1* randn(512, 512));
q2 = (5* randn(512, 512));
q3 = (10* randn(512, 512));
std_noise1a = std2(q1); % which is basically 1 
std_noise2a = std2(q2);
std_noise3a = std2(q3);
std_noise1e = noise_estimation(noise1);
std_noise2e = noise_estimation(noise2);
std_noise3e = noise_estimation(noise3);

%estimation ratio -- sigma estimated/simga added
e1 = std_noise1e/std_noise1a;
e2 = std_noise2e/std_noise2a;
e3 = std_noise3e/std_noise3a;


%% brainmri noise metrics 
std_brain1 = noise_estimation(brainMRI1);
std_brain2 = noise_estimation(brainMRI2);
std_brain3 = noise_estimation(brainMRI3);
std_brain4 = noise_estimation(brainMRI4);
std_brain5 = noise_estimation(brainMRI5);
std_brain6 = noise_estimation(brainMRI6);

x = [std_brain1, std_brain2, std_brain3, std_brain4, std_brain5, std_brain6];
y = (1:6);

scatter(y,x);


% % making a noise image of standard deviation of 10 gray level
% % we can range it from 0-50
% noiseOnlyImage = 1* randn(H, W); % we added 1 std 
% % adding noise image to gray scale image 
% noiseAddedImage = double(noiseCT1)+noiseOnlyImage;
% % figure;
% % imshow(noiseAddedImage, []);
% % figure;
% % imshow(noiseOnlyImage, []);
% %compute std of noisy image
% [std_addednoise_CT1, Ttest] = noise_estimation(noiseAddedImage);
% % Ttest = noise_estimation(noiseAddedImage);
% 
% estimation_ratio = std_noiseCT1/std_addednoise_CT1;


%% Contrast Quality Metrics 

%% Step 1 Find local measure names edge-based criterion (ECC)

ECC_CT1 = ECC(brainMRI1);
ECC_CT2 = ECC(brainMRI6);

ECC_XY = 2/(1+(1-ECC_CT2)/(1-ECC_CT1));

%% Step 2 Find Entropy 
% MaskedCT =  apply_mask(training_post_1, mask, 'ImageType');
% mask = MaskedCT(:,:,143);
entropy1 = entropy(brainMRI1);
entropy2 = entropy(brainMRI6);

H_XY = 2/(1+(entropy1/entropy2));


%% Step 3 Correlation coeff

imgcorr = corrcoef(brainMRI1(:), brainMRI6(:));

R_XY = imgcorr(1,2);

%% Step 4 Ambsolute Mean Brightness Error

AMBE_XY = (2.*(mean(brainMRI1(:))).*(mean(brainMRI2(:))))/((mean(brainMRI1(:)).^2)+(mean(brainMRI6(:)).^2));

%% Step 5 Linear Combination 

% Creating a function handle of the linear combination 
fun = @(A) A(1)*ECC_XY + A(2)*H_XY + A(3)*R_XY + A(4)*AMBE_XY;
f = @(A) A*4;
nvars = 4;
x = particleswarm(fun, nvars,0, 1)


CCIQ = (ECC_XY+H_XY+R_XY+AMBE_XY)/4;

%% Testing 
yref = ECC_XY + H_XY + R_XY + AMBE_XY;
M = [ECC_XY(:), H_XY(:), R_XY(:), AMBE_XY(:)];
b = yref(:);
ABCD = M\b;
A = ABCD(1); B = ABCD(2); C = ABCD(3); D = ABCD(4);

