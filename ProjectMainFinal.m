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


%% Initializing Brain MRI2 Images 
% Initializing mammogram files 
brainMRI1 = load('brainMRI_1.mat'); brainMRI1 = brainMRI1.vol;brainMRI1 = brainMRI1(:,:,90);
brainMRI2 = load('brainMRI_2.mat'); brainMRI2 = brainMRI2.vol;brainMRI2 = brainMRI2(:,:,90);
brainMRI3 = load('brainMRI_3.mat'); brainMRI3 = brainMRI3.vol;brainMRI3 = brainMRI3(:,:,90);
brainMRI4 = load('brainMRI_4.mat'); brainMRI4 = brainMRI4.vol;brainMRI4 = brainMRI4(:,:,90);
brainMRI5 = load('brainMRI_5.mat'); brainMRI5 = brainMRI5.vol;brainMRI5 = brainMRI5(:,:,90);
brainMRI6 = load('brainMRI_6.mat'); brainMRI6 = brainMRI6.vol;brainMRI6 = brainMRI6(:,:,90);

%% Image Quality Metric 1: Noise 
 
%% LungCT
std_CT1 = noise_estimation(noiseCT1); 
std_CT10 = noise_estimation(noiseCT10); 
std_CThalf = noise_estimation(noiseCT_half);

xct = 1:3;
yct = [std_CThalf, std_CT1, std_CT10];

figure; 
scatter(xct(1),yct(1));
hold on;
scatter(xct(2),yct(2));
scatter(xct(3),yct(3));
legend('CT 0.5x noise', 'CT 1x noise', 'CT 10x noise');
xlabel('LungCT image');
ylabel('Standard Deviation');
title('Noise Quality Metric of LungCT Images');
axis tight;
hold off;

%% Brain MRI2
std_brain1 = noise_estimation(brainMRI1);
std_brain2 = noise_estimation(brainMRI2);
std_brain3 = noise_estimation(brainMRI3);
std_brain4 = noise_estimation(brainMRI4);
std_brain5 = noise_estimation(brainMRI5);
std_brain6 = noise_estimation(brainMRI6);

x = [std_brain1, std_brain2, std_brain3, std_brain4, std_brain5, std_brain6];
y = (1:6);

figure;
scatter(y(1),x(1));
hold on;
scatter(y(2),x(2));
scatter(y(3),x(3));
scatter(y(4),x(4));
scatter(y(5),x(5));
scatter(y(6),x(6));
legend('brainMRI1', 'brainMRI2', 'brainMRI3', 'brainMRI4', 'brainMRI5', 'brainMRI6');
xlabel('Brain MRI2 image');
ylabel('Standard Deviation');
title('Noise Quality Metric of Brain MRI2 Images');
hold off;

%% LungCT and Brain
x1 = (1:9);
y1 = [std_brain1, std_brain2, std_brain3, std_brain4, std_brain5, std_brain6, std_CThalf, std_CT1, std_CT10];

figure; 
scatter(x1(1),y1(1));
hold on;
scatter(x1(2),y1(2));
scatter(x1(3),y1(3));
scatter(x1(4),y1(4));
scatter(x1(5),y1(5));
scatter(x1(6),y1(6));
scatter(x1(7),y1(7));
scatter(x1(8),y1(8));
scatter(x1(9),y1(9));
legend('brainMRI1', 'brainMRI2', 'brainMRI3', 'brainMRI4', 'brainMRI5', 'brainMRI6','CT 0.5x noise', 'CT 1x noise', 'CT 10x noise');
xlabel('LungCT image');
ylabel('Standard Deviation');
title('Noise Quality Metric of LungCT Images');
axis tight;
hold off;
%% Image Quality Metric 2: Contrast

%% LungCT 
CCIQCT1 = contrast_estimation(noiseCT_half, noiseCT1);
CCIQCT2 = contrast_estimation(noiseCT_half, noiseCT10);

x3 = (1:2);
y3 = [CCIQCT1, CCIQCT2];

figure; 
scatter(x3(1),y3(1));
hold on;
scatter(x3(2),y3(2));
legend('CT 1 & 0.5', 'CT 0.5 & 1', 'CT 0.5 & 10');
xlabel('Lung CT Image');
ylabel('CCIQ');
title('Contrast Quality Metric of Brain MRI Images');
hold off;
%% Brain MRI 2 

CCIQ1 = contrast_estimation(brainMRI1, brainMRI2);
CCIQ2 = contrast_estimation(brainMRI1, brainMRI3); 
CCIQ3 = contrast_estimation(brainMRI1, brainMRI4);
CCIQ5 = contrast_estimation(brainMRI1, brainMRI5);
CCIQ6 = contrast_estimation(brainMRI1, brainMRI6);

montage({cast(brainMRI1, 'uint8'), cast(brainMRI2, 'uint8'), cast(brainMRI3, 'uint8'), cast(brainMRI4, 'uint8'), cast(brainMRI5, 'uint8'), cast(brainMRI6, 'uint8')});

x2 = (1:5);
y2 = [CCIQbrain1, CCIQbrain2, CCIQbrain3, CCIQbrain4, CCIQbrain5];

figure; 
scatter(x2(1),y2(1));
hold on;
scatter(x2(2),y2(2));
scatter(x2(3),y2(3));
scatter(x2(4),y2(4));
scatter(x2(5),y2(5));
legend('brainMRI1&2', 'brainMRI1&3', 'brainMRI1&4', 'brainMRI1&5', 'brainMRI1&6');
xlabel('Brain MRI 2 Image');
ylabel('CCIQ');
title('Contrast Quality Metric of Brain MRI Images');
hold off;


%% Image Quality Metric 3


%%% ORIGINAL %%%

MRI_12 = imageQuality_edge(brainMRI1, brainMRI2);
MRI_13 = imageQuality_edge(brainMRI1, brainMRI3);
MRI_14 = imageQuality_edge(brainMRI1, brainMRI4);
MRI_15 = imageQuality_edge(brainMRI1, brainMRI5);
MRI_16 = imageQuality_edge(brainMRI1, brainMRI6);

x = [MRI_12, MRI_13, MRI_14, MRI_15, MRI_16];
y = (1:5);

figure;
scatter(y(1),x(1));
hold on;
scatter(y(2),x(2));
scatter(y(3),x(3));
scatter(y(4),x(4));
scatter(y(5),x(5));
legend('MRI 1 and 2', 'MRI 1 and 3', 'MRI 1 and 4', 'MRI 1 and 5', 'MRI 1 and 6');
title('MSE Values of Original MRIs');
ylabel('Mean Squared Error');
xlabel('MRI Image #');
hold off;

%%% FILTERED MSE %%%
MRI1_filt = filter2(fspecial('average',3),brainMRI1)/255;
MRI2_filt = filter2(fspecial('average',3),brainMRI2)/255;
MRI3_filt = filter2(fspecial('average',3),brainMRI3)/255;
MRI4_filt = filter2(fspecial('average',3),brainMRI4)/255;
MRI5_filt = filter2(fspecial('average',3),brainMRI5)/255;
MRI6_filt = filter2(fspecial('average',3),brainMRI6)/255;

MRI_12f = imageQuality_edge(MRI1_filt, MRI2_filt);
MRI_13f = imageQuality_edge(MRI1_filt, MRI3_filt);
MRI_14f = imageQuality_edge(MRI1_filt, MRI4_filt);
MRI_15f = imageQuality_edge(MRI1_filt, MRI5_filt);
MRI_16f = imageQuality_edge(MRI1_filt, MRI6_filt);

xf = [MRI_12f, MRI_13f, MRI_14f, MRI_15f, MRI_16f];
yf = (1:5);

figure;
scatter(yf(1),xf(1));
hold on;
scatter(yf(2),xf(2));
scatter(yf(3),xf(3));
scatter(yf(4),xf(4));
scatter(yf(5),xf(5));
legend('MRI 1 and 2', 'MRI 1 and 3', 'MRI 1 and 4', 'MRI 1 and 5', 'MRI 1 and 6');
title('MSE Values of Filtered MRIs');
ylabel('Mean Squared Error');
xlabel('MRI Image #');
hold off;

%%% COMBINED LINE GRAPH FOR COMPARISON %%%

figure;
plot(y, x);
hold on;
plot(yf,xf);
legend;
title('MSE Values of Original and Filtered MRIs');
subtitle('Comparison of Edge Quality');
ylabel('Mean Squared Error');
xlabel('MRI Image #');
hold off;


%% Troubleshoot contrast 


% CCIQbrain12 = contrast_estimation(brainMRI1, brainMRI2);
% CCIQbrain21 = contrast_estimation(brainMRI2, brainMRI1); 
% CCIQbrain13 = contrast_estimation(brainMRI1, brainMRI3);
% CCIQbrain31 = contrast_estimation(brainMRI3, brainMRI1); 
% CCIQbrain23 = contrast_estimation(brainMRI2, brainMRI3); 
% CCIQbrain32 = contrast_estimation(brainMRI3, brainMRI2);
% CCIQbrain24 = contrast_estimation(brainMRI2, brainMRI4);
% CCIQbrain42 = contrast_estimation(brainMRI4, brainMRI2);
% CCIQbrain34 = contrast_estimation(brainMRI3, brainMRI4);
% CCIQbrain43 = contrast_estimation(brainMRI4, brainMRI3);
% CCIQbrain14 = contrast_estimation(brainMRI1, brainMRI4);
% CCIQbrain41 = contrast_estimation(brainMRI4, brainMRI1);
% CCIQbrain25 = contrast_estimation(brainMRI2, brainMRI5);
% CCIQbrain52 = contrast_estimation(brainMRI5, brainMRI2);
% CCIQbrain45 = contrast_estimation(brainMRI4, brainMRI5);
% CCIQbrain54 = contrast_estimation(brainMRI5, brainMRI4);
% CCIQbrain35 = contrast_estimation(brainMRI3, brainMRI5);
% CCIQbrain53 = contrast_estimation(brainMRI5, brainMRI3);
% CCIQbrain15 = contrast_estimation(brainMRI1, brainMRI5);
% CCIQbrain51 = contrast_estimation(brainMRI5, brainMRI1);
% CCIQbrain26 = contrast_estimation(brainMRI2, brainMRI6);
% CCIQbrain62 = contrast_estimation(brainMRI6, brainMRI2);
% CCIQbrain65 = contrast_estimation(brainMRI6, brainMRI5);
% CCIQbrain56 = contrast_estimation(brainMRI5, brainMRI6);

% contrast in order of highest to lowest -> brain mri 2, brain mri 3,
% brainmri1, brainmri4, brainmri5, brainmri6