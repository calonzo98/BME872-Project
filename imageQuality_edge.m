function [MSE] = imageQuality_edge(img1,img2)

%%%%%%%Quality Assessment Metrics for Edge Detection and
%Edge-aware Filtering: A Tutorial Review
%Diana Sadykova, and Alex Pappachen James
%School of Engineering
%Nazarbayev University, Astana

%Sobel Edge Detection%

G1 = edge(img1, 'Sobel');
G2 = edge(img2, 'Sobel');

sobel_edge_img1 = cast(G1, 'uint8');
sobel_edge_img2 = cast(G2, 'uint8');

%MSE
if isempty(sobel_edge_img1) % If x is empty, y must also be empty
    MSE = [];
    return;
end

if isinteger(sobel_edge_img1)     
    sobel_edge_img1 = double(sobel_edge_img1);
    sobel_edge_img2 = double(sobel_edge_img2);
end

%Mean Square Error Equation%
MSE = (norm(sobel_edge_img1(:)-sobel_edge_img2(:),2).^2)/numel(sobel_edge_img1);

fprintf('\n The mean-squared error is %0.4f\n', MSE); 

end


