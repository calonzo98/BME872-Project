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
figure;
imshow(noiseCT_half_slice, [])
imdisplayrange;
colorbar6
title('CT Image with .5x noise ');
test = cast(noiseCT_half, 'uint8');
figure;
[bins, freq]=intensityHistogram(noiseCT_half, 5,1);

%% Initializing Brain MRI2 Images 
% Initializing mammogram files 
brainMRI1 = load('brainMRI_1.mat'); brainMRI1 = brainMRI1.vol;brainMRI1 = brainMRI1(:,:,90);
brainMRI2 = load('brainMRI_2.mat'); brainMRI2 = brainMRI2.vol;brainMRI2 = brainMRI2(:,:,90);
brainMRI3 = load('brainMRI_3.mat'); brainMRI3 = brainMRI3.vol;brainMRI3 = brainMRI3(:,:,90);
brainMRI4 = load('brainMRI_4.mat'); brainMRI4 = brainMRI4.vol;brainMRI4 = brainMRI4(:,:,90);
brainMRI5 = load('brainMRI_5.mat'); brainMRI5 = brainMRI5.vol;brainMRI5 = brainMRI5(:,:,90);
brainMRI6 = load('brainMRI_6.mat'); brainMRI6 = brainMRI6.vol;brainMRI6 = brainMRI6(:,:,90);


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
imshow(G, []);
figure;
imhist(noiseCT1);

%% Step 2: Threshold for edge map test1
% [counts, bins] = imhist(brainMRI1);
% cdf = cumsum(counts) / sum(counts);
% plot(cdf);

%% Step 2: Threshold for edge map test1

%% Ploting CDF  -- works like a charm!!!
Gtest = cast(G, 'uint8');
[hist, normhist] = histo_norm(Gtest);
cdf1 = hist_cum(normhist);
figure;
plot(cdf1);
%% creating a loop to find intensity x when cdf y=0.1 (i.e 10% of entire image intensities)
[r,c] = size(cdf1);
for i = 1:r
    for j = 1:c
        if (cdf1(i,j)>=0.1)&&(cdf1(i,j)<=0.101)
            T = j;
        end
    end
end

%% find edge map 

Gth = image_threshold(G,T);
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
%% Step 3 Laplacian 
% multiply Gth with laplacian to suppress image stuctures 
kernel = [1 -2 1; -2 4 -1; 1 -2 1];
Gth_lap = conv2(Gth, kernel);
imshow(Gth_lap, []);

%% Step 4 Calculating standard deviation
std_noiseCT1 = sqrt(pi/2)*(1/(6*(512-2)*(512-2))).*sum(abs(Gth_lap), [1 2]);


%% test -- adding noise to noiseCT1
% making a noise image og standard deviation of 10 gray levels 
noiseOnlyImage = 10 * randn(512, 512);
% adding noise image to gray scale image 
noiseAddedImage = double(noiseCT1)+noiseOnlyImage;

%compute std of noisy image
[std_addednoise_CT1, Ttest] = noise_estimation(noiseAddedImage);


estimation_ratio = std_noiseCT1/std_addednoise_CT1;


