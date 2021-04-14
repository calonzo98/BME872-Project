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


%% initialize image size 

[H, W] = size(brainMRI1);

%% Noise Quality Metric 

%% Step 1 edge detection -> using fspecial 
h = fspecial('sobel');
v = rot90(fspecial('sobel'));
gh = imfilter(brainMRI1, h);
gv = imfilter(brainMRI1, v);
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
%Ttest = auto_thresholding(G, 1);

%% find edge map -- exclude pixels less than T 
Gth = image_threshold(G,T);
imshow(Gth,[]);


%% Step 3 Laplacian 
% multiply Gth with laplacian to suppress image stuctures 
kernel = [1 -2 1; -2 4 -1; 1 -2 1];
img_lap = double(imfilter(Gth, kernel));
imshow(img_lap, []);

% kernel = [1 -2 1; -2 4 -1; 1 -2 1];
% Gth_lap = conv2(Gth, kernel);
% imshow(Gth_lap, []);

%% Step 4 Calculating standard deviation
std_test = sqrt(pi/2)*(1/(6*(H-2)*(W-2))).*sum(abs(img_lap),'all');

%std_noiseCT1 = sqrt(pi/2)*(1/(6*(H-2)*(W-2))).*mean(abs(img_lap(:)));

 
%% Plotting STD 
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
axis tight;
hold off;


%% brainmri noise metrics 
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
hold off;


%% Contrast Quality Metric 

%% Step 1 Find local measure names edge-based criterion (ECC)

ECC_CT1 = ECC(brainMRI1);
ECC_CT2 = ECC(brainMRI2);

ECC_XY = 2/(1+(1-ECC_CT2)/(1-ECC_CT1));

%% Step 2 Find Entropy 
% MaskedCT =  apply_mask(training_post_1, mask, 'ImageType');
% mask = MaskedCT(:,:,143);
brain = cast(brainMRI1, 'uint8');
entropy1 = entropy(brain);
entropy2 = entropy(cast(brainMRI2, 'uint8'));

H_XY = 2/(1+(entropy1/entropy2));


%% Step 3 Correlation coeff

imgcorr = corrcoef(brainMRI1(:), brainMRI6(:));

R_XY = imgcorr(1,2);

%% Step 4 Ambsolute Mean Brightness Error

AMBE_XY = (2.*(mean(brainMRI1(:))).*(mean(brainMRI2(:))))/((mean(brainMRI1(:)).^2)+(mean(brainMRI6(:)).^2));

%% Step 5 Linear Combination

a1 = optimvar('a1');
a2 = optimvar('a2');
a3 = optimvar('a3');
a4 = optimvar('a4');
prob = optimproblem;
prob.Objective = a1*ECC_XY + a2*H_XY + a3*R_XY + a4*AMBE_XY;
prob.Constraints.cons1 = a1+a2+a3+a4==1;
prob.Constraints.cons2 = a1>=0;
prob.Constraints.cons3 = a1<=1;
prob.Constraints.cons4 = a2>=0;
prob.Constraints.cons5 = a2<=1;
prob.Constraints.cons6 = a3>=0;
prob.Constraints.cons7 = a3<=1;
prob.Constraints.cons8 = a4>=0;
prob.Constraints.cons9 = a4<=1;

sol = solve(prob, 'Solver', 'intlinprog')
alpha = vertcat(sol.vector);
%% 
CCIQ = (sol.a1*ECC_XY) + (sol.a2*H_XY) + (sol.a3*R_XY) + (sol.a4*AMBE_XY);
 
%% 
% CCIQ = [ECC_XY H_XY R_XY AMBE_XY];
% % finding alpha 
% alpha = linprog(-CCIQ.',[],[],[1 1 1 1], 1, [0 0 0 0],[1 1 1 1]);
%% Testing 
yref = ECC_XY + H_XY + R_XY + AMBE_XY;
M = [ECC_XY(:), H_XY(:), R_XY(:), AMBE_XY(:)];
b = yref(:);
ABCD = M\b;
A = ABCD(1); B = ABCD(2); C = ABCD(3); D = ABCD(4);

%% CCIQ Values 

J = imadjust(cast(brainMRI1, 'uint8'),[],[],0.5);
J = double(J);
brainmri_gaus = imgaussfilt(cast(brainMRI1, 'uint8'));
brainMRI_stretch= histeq(brainmri_gaus);

brainMRI_stretch = double(brainMRI_stretch);
figure;
subplot(1,2,1)
imshow(J, []);
subplot(1,2,2)
imshow(brainMRI1, []);

%%

CCIQ1 = contrast_estimation(brainMRI1, brainMRI1); % 1.0799 that means 2 is higher contrast than 1
CCIQ2 = contrast_estimation(brainMRI1, brainMRI3); % 
CCIQ3 = contrast_estimation(brainMRI1, brainMRI4);
CCIQ5 = contrast_estimation(brainMRI1, brainMRI5);
CCIQ6 = contrast_estimation(brainMRI1, brainMRI6);

montage({cast(brainMRI1, 'uint8'), cast(brainMRI2, 'uint8'), cast(brainMRI3, 'uint8'), cast(brainMRI4, 'uint8'), cast(brainMRI5, 'uint8'), cast(brainMRI6, 'uint8')});