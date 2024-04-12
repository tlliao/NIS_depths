function [rmse_error] = calcDistofFund(matches1, matches2, fund_H, fund_e, depth_img)
% calculate the mapping error via fundamental matrix H_inf and epipolar e'
% x' ~ H_inf * x+e'/z
num_pts = size(matches1,2);
depth_matches1 = calcDepthofPoint(matches1, depth_img);
map_x = (fund_H(1,:)*[matches1; ones(1,num_pts)]+fund_e(1)./depth_matches1);
map_y = (fund_H(2,:)*[matches1; ones(1,num_pts)]+fund_e(2)./depth_matches1);
map_z = (fund_H(3,:)*[matches1; ones(1,num_pts)]+fund_e(3)./depth_matches1);
map_matches2 = [map_x./map_z; map_y./map_z];

rmse_error = sqrt(sum((map_matches2 - matches2).^2, 'all')/num_pts);

end

