function CCIQ = contrast_estimation(img1, img2)

%Step 1 Find local measure names edge-based criterion (ECC)

ECC_CT1 = ECC(img1);
ECC_CT2 = ECC(img2);

ECC_XY = 2/(1+(1-ECC_CT2)/(1-ECC_CT1));

img1_u = cast(img1, 'uint8');
img2_u = cast(img2, 'uint8');
% Step 2 Find Entropy 
entropy1 = entropy(img1_u);
entropy2 = entropy(img2_u);
H_XY = 2/(1+(entropy1/entropy2));

% Step 3 Correlation coeff

imgcorr = corrcoef(img1(:), img2(:));

R_XY = imgcorr(1,2);

% Step 4 Ambsolute Mean Brightness Error

AMBE_XY = (2.*(mean(img1(:))).*(mean(img2(:))))/((mean(img1(:)).^2)+(mean(img2(:)).^2));

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

sol = solve(prob, 'Solver', 'intlinprog');

CCIQ = (sol.a1*ECC_XY) + (sol.a2*H_XY) + (sol.a3*R_XY) + (sol.a4*AMBE_XY);

end 


