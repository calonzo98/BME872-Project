% The function [std] = noise_estimation(img_in) estimates noise of an imag using laplacian operator and adaptive
% edge detection. The function returns the standard deviatio representing
% noise estimations.
%
%
%By using an adaptive edge detector, noise can be estimated for a large
%range of noise levels including highly textured images. 
%
%A sobel operator is used for edge detection. To decide the edge map, the
%threshold value is selected when the accumulated histogram of the gradient
%magnitude reaches 10% of the entire image. The threshold values G will be
%different for each image. Therefore, we can say that it is an "adaptive"
%edge detection. The edge pixels below the threshold are elliminated, and
%the rest of the image is used for the remainder of the computation. The
%edge detector exclude image details from ontributing to the noise
%variance. 
%
%The next step is to suppress the image structures by a laplacian oeprator
%
%Lastly, the standard devition of the noise is computed 
%

function [std] = noise_estimation(img_in)

% Initialize size of image
[H, W] = size(img_in);

% Step 1: edge detection
h = fspecial('sobel');
v = rot90(fspecial('sobel'));
gh = imfilter(img_in, h);
gv = imfilter(img_in, v);
[G, gdir] = imgradient(gh,gv);

%Step 2: edge map
% computing cdf histo
[hist, normhist] = histo_norm(cast(G, 'uint8'));
cdf = hist_cum(normhist);
figure;
plot(cdf);

%determining threshold value 
[r,c] = size(cdf);
for i = 1:r
    for j = 1:c
        if (cdf(i,j)>=0.08)&&(cdf(i,j)<=0.13)
            T = j;
        else 
            disp('error');
        end
    end
end

%T =auto_thresholding(G, 1);
%computing edge map

Gth = image_threshold(G,T);
imshow(Gth,[]);

%Step 3: Apply Laplacian 
kernel = [1 -2 1; -2 4 -1; 1 -2 1];
Gth_lap = imfilter(Gth, kernel);
imshow(Gth_lap, []);

%Step 4: Compute std of noise
std = sqrt(pi/2)*(1/(6*(W-2)*(H-2))).*sum(abs(Gth_lap),'all');
end 