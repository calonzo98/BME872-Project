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
[training_post_1, infoCT1] = imageRead('C:\Users\Claudia\Documents\MATLAB\BME872\Lab1\Lab1-LungCT\noise_0.5x_post.mhd', '.mhd');
noiseCT1 = training_post_1.data;
noiseCT1=double(noiseCT1(:,:,143));
[training_post_10, infoCT10] = imageRead('C:\Users\Claudia\Documents\MATLAB\BME872\Lab1\Lab1-LungCT\noise_0.5x_post.mhd', '.mhd');
noiseCT10 = training_post_10.data;
noiseCT10=noiseCT10(:,:,143);

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
gh = imfilter(noiseCT1, h);
gv = imfilter(noiseCT1, v);
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
imhist(noiseCT1);

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
%% find edge map 
Gth = image_threshold(G,Ttest);
imshow(Gth,[]);

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
img_lap = imfilter(Gth, kernel);
imshow(img_lap, []);

% kernel = [1 -2 1; -2 4 -1; 1 -2 1];
% Gth_lap = conv2(Gth, kernel);
% imshow(Gth_lap, []);

%% Step 4 Calculating standard deviation
std_noiseCT1 = sqrt(pi/2)*(1/(6*(H-2)*(W-2))).*sum(abs(img_lap), [1 2]);


%% test -- adding noise to noiseCT1
% ultimate goal is to determine which image has the most noise 

%% Performance in estimating the standard deviation 
% add noise to CT1 and find estimation ratio
std_CT1 = noise_estimation(noiseCT1); % assume this to be base image (i.e no added noise)
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
