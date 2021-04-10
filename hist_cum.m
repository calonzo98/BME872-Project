function [c_hist]=hist_cum(hp)
 [r,c]=size(hp);
 c_hist=zeros(1,256);
  for j=1:c
       if(j ==1)
       c_hist(j)=hp(j);
       else
       c_hist(j)=hp(j)+c_hist(j-1);         
       end
  end
end