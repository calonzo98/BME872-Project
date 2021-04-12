function [img_out,t] = contrastStretch(img_in)

img_in = double(img_in);
% Find rmin and rmax
rmin = min(img_in(:));
rmax = max(img_in(:));

img_out = zeros(size(img_in)); %initializing img_out matrix 

% Transfer Function 
t=[];
for r=rmin:rmax
    t= [t (255.*((1./(rmax-rmin).*r)-(rmin./(rmax-rmin))))];
end 

for ii = 1:size(img_in, 1) % number of rows of the image
    for jj = 1:size(img_in, 2) % for number of columns in the image 
        %r=img_in(ii,jj); %get pixel value 
        %r_new = 255.*((1./(rmax-rmin).*r)-(rmin./(rmax-rmin)));
        img_out(ii,jj)=t(img_in(ii,jj)+1); % Save new pixel value in image out 
    end 
end

% Plotting Transfer Function 
% L=0:rmax;
% plot(L,t);
% title ('Contrast Stretch Transfer Function')
% xlabel ('Input Intensity Level, r')
% ylabel ('Output Intensity level, s')
end 
