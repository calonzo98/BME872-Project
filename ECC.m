% this function takes in an image and outputs the ECC value. 
% 
%the mean edge gray level E at pixel (i,j) is computed by the summation of
%the multiplication of the intensity values and edge pixels divided by the
%sum of the edge values in the neighboorhood. where the edge pixels are
%found using soble edge detection. 
%
%cij calculates the contrast for each pixel 
%
%ecc is then calculated by taking the average of cij
%

function EB=ECC(img)
[m,n]=size(img);
h_sob=fspecial('sobel');
img_norm=double(img./255);
img_filt_norm=imfilter(img_norm,h_sob,'replicate');
gij_xij=(img_norm).*(img_filt_norm);
h_avg=fspecial('average'); %%%default size is 3 by 3
img_filt_avg=imfilter(img_filt_norm,h_avg,'replicate');
gij_xij_avg=imfilter(gij_xij,h_avg,'replicate');
eij=(gij_xij_avg./(img_filt_avg+0.0001));
cij=abs(img_norm-eij)./abs(img_norm+eij+0.0001);
cij_u=sum(sum(uint8(round(cij*255))));
EB=cij_u/(m*n);
end