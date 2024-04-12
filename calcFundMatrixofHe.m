function [H_inf, fund_e] = calcFundMatrixofHe(pts1, pts2, depth_img1)
% give feature correspondences and depths of target image, evaluate
% the fundamental matrix F which F=[e']_x * H_inf
% 1. x'_i ~ (H_inf*x_i + e'/z_i)
% Inputs:
% pts1: 2*N, ptx, pty,
% pts2: 2*N, ptx', pty'
% depth_img1: depth map of target image
% -------------------------------------
%% perform coordinate normalization
depth_pts1 = calcDepthofPoint(pts1, depth_img1);
[normalized_pts1, T1] = normalise2dpts([pts1; ones(1,size(pts1, 2))]);
[normalized_pts2, T2] = normalise2dpts([pts2; ones(1,size(pts2, 2))]);
normalized_pts1 = normalized_pts1(1:2,:);
normalized_pts2 = normalized_pts2(1:2,:);

%% calculation the initial H0 and e0
Equation_matrix = zeros(2*size(normalized_pts1, 2), 12);
for i=1:size(normalized_pts1, 2)
    zi = depth_pts1(i); xi=normalized_pts1(1,i); yi=normalized_pts1(2,i);
    xi_=normalized_pts2(1,i); yi_=normalized_pts2(2,i);
    tmp_coeff1 = [xi, yi, 1, 0, 0, 0, -xi*xi_, -yi*xi_, -xi_, 1/zi, 0, -xi_/zi];
    tmp_coeff2 = [0, 0, 0, xi, yi, 1, -xi*yi_, -yi*yi_, -yi_, 0, 1/zi, -yi_/zi];
    Equation_matrix(2*i-1:2*i, :) = [tmp_coeff1; tmp_coeff2];
end

[~,~,v] = svd(Equation_matrix, 0);
norm_H0 = reshape(v(1:9, 12), 3, 3)';
norm_e0 = v(10:12, 12);

H_inf0 = T2\(norm_H0*T1);   fund_e0 = T2\norm_e0;
x0 = [H_inf0(:); fund_e0 ]; 

%% calculate optimal H_inf and e'

options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt', 'Display', 'off', 'FunctionTolerance', 1e-8, 'MaxIterations', 1e3 );
[x, ~] = lsqnonlin(@(x)calcHe(x, [pts1; depth_pts1], pts2), x0,[],[],options);

H_inf = reshape(x(1:9),3,3);  fund_e = x(10:12);

end


function F = calcHe(x, pts1, pts2) % least-square cost function

F = zeros(1, 2*size(pts1,2));
H = [x(1), x(4), x(7); x(2), x(5), x(8); x(3), x(6), x(9)];
e = [x(10); x(11); x(12)];
for i=1:size(pts1,2)
    map_pts1 = (H(1:2,:)*[pts1(1:2,i);1]+e(1:2)/pts1(3,i))./(H(3,:)*[pts1(1:2,i);1]+e(3)/pts1(3,i));
    F(2*i-1) = (map_pts1(1)-pts2(1,i));
    F(2*i) = (map_pts1(2)-pts2(2,i));
end

end

