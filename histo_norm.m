function [hist, hist_P]=histo_norm(im)
hist=zeros(1,256);
[r, c]=size(im);
s=r*c;
for i=1:r
    for j=1:c
        int_val=im(i,j);
        hist(int_val+1)=hist(int_val+1)+1;
    end
end
hist_P=hist./s; %Normalized Histogram(PDF)
end
