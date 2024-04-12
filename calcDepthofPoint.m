function [depth_pts] = calcDepthofPoint(pts, depth_img)
% calculate the depths of feature points based on bilinear interpolation on
% source depth image
depth_pts = zeros(1, size(pts,2));
for i=1:size(pts,2)
    tmp_depth = depthBilinear(pts(:,i), depth_img);
    depth_pts(i) = tmp_depth;
end

end

