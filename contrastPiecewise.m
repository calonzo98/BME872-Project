function img_out = contrastPiecewise (img_in, a, b)

img_in = cast(img_in, 'double');
img_out = zeros(size(img_in)); %initializing img_out matrix 

% Defining points and slope of each line in piecewise graph
% a is the first point vector and b is the second point vector 
start=[0,0];
stop=[255,255];
m1=(start(1,2)-a(1,2))/(start(1,1)-a(1,1)); %slope of line 1 
m2=(a(1,2)-b(1,2))/(a(1,1)-b(1,1)); % Slope of line 2
m3=(b(1,2)-stop(1,2))/(b(1,1)-stop(1,1)); %Slope of line 3
c1=start(1,2)-(m1*start(1,1)); %constant for slope 1
c2=a(1,2)-(m2*a(1,1)); % constant for slope 2
c3=b(1,2)-(m3*b(1,1)); % constant for slope 3

% Transformation Function 
t=[];
for x=0:255
    if(x<a(1,1))
        t=[t (m1*x+c1)];
    end
    if(x>=a(1,1) && x<=b(1,1))
        t=[t (m2*x+c2)];
    end
    if(x>b(1,1) && x<=stop(1,1))
        t=[t (m3*x+c3)];
    end
end
    
% Filling img_out matrix 

for ii = 1:size(img_in, 1) % number of rows in image 
    for jj = 1:size(img_in, 2) % number of columng in image 
        img_out(ii,jj)=t(img_in(ii,jj)+1);
    end
end

% Plotting Transfer Function
figure;
% L=0:255;
plot(t);
title ('Piecewise Linear Transfer Function')
xlabel ('Input Intensity Level, r')
ylabel ('Output Intensity level, s')

end 