function [std, T] = noise_estimation(img_in)

% Step 1: edge detection
h = fspecial('sobel');
v = rot90(fspecial('sobel'));
gh = imfilter(img_in, h);
gv = imfilter(img_in, v);
[G, gdir] = imgradient(gh,gv);

%Step 2: edge map
% computing cdf histo
Gtest = cast(G, 'uint8');
[hist, normhist] = histo_norm(Gtest);
cdf = hist_cum(normhist);
figure;
plot(cdf);

%determining threshold value 
[r,c] = size(cdf);
for i = 1:r
    for j = 1:c
        if (cdf(i,j)>=0.1)&&(cdf(i,j)<=0.101)
            T = j;
        end
    end
end

%computing edge map
Gth = image_threshold(G,T);
imshow(Gth,[]);

%Step 3: Apply Laplacian 
kernel = [1 -2 1; -2 4 -1; 1 -2 1];
Gth_lap = conv2(Gth, kernel);
imshow(Gth_lap, []);

%Step 4: Compute std of noise
std = sqrt(pi/2)*(1/(6*(512-2)*(512-2))).*sum(abs(Gth_lap), [1 2]);
end 