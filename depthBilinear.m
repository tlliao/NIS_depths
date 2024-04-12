function [ depth_pts ] = depthBilinear(pts, depth_img)
% calculate the depth of an non-integer point via bilinear interpolation
ptx = pts(1);
pty = pts(2);
tmp_fx = floor(ptx);
tmp_fy = floor(pty);
tmp_cx = floor(ptx)+1;
tmp_cy = floor(pty)+1;
left_top = [tmp_fx, tmp_fy];
left_bot = [tmp_fx, tmp_cy];
right_top = [tmp_cx, tmp_fy];
right_bot = [tmp_cx, tmp_cy];
mesh_pts = [left_top; right_top; right_bot; left_bot];
coeff_mesh = meshGridAlign(mesh_pts, [ptx, pty]);
if tmp_cx>size(depth_img,2)
    tmp_cx = tmp_fx;
end
if tmp_cy>size(depth_img,1)
    tmp_cy = tmp_fy;
end
depth_pts = coeff_mesh'*[depth_img(tmp_fy, tmp_fx); depth_img(tmp_fy, tmp_cx);...
depth_img(tmp_cy, tmp_cx); depth_img(tmp_cy, tmp_fx)]; 

end

