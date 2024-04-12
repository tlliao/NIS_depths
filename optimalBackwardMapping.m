function [warped_img1, warped_img2, average_panorama, quality_errors]=optimalBackwardMapping(img1,img2,depth_img1,depth_img2,H_inf1,fund_e1, H_inf2, fund_e2)
% rendering the warped img1 and img2 using backward mapping
% generate full warped pixels from img1 to img2 using H_inf1 and fund_e1,
% then rendered warped img1 backward using bilinear interpolation on each
% grid mesh, refining the warped result to remove occlusion pixels and
% overlapped pixels;
%% pre-process
[sz1,sz2]=size(depth_img1);
para_dist_disp=0.01*sqrt(sz1^2+sz2^2); % threshold of pixel dispalcement to refine the overlapping region
C2=sz2; C1=sz1;

%% calculate the canvas size of panorama
[X, Y]=meshgrid(linspace(1,sz2,C2),linspace(1,sz1,C1));
full_pts=[X(:),Y(:)];
depth_pts=calcDepthofPoint(full_pts', depth_img1);
tmp1=H_inf1(1,:)*[full_pts'; ones(1,size(full_pts,1))]+fund_e1(1)./depth_pts;
tmp2=H_inf1(2,:)*[full_pts'; ones(1,size(full_pts,1))]+fund_e1(2)./depth_pts;
tmp3=H_inf1(3,:)*[full_pts'; ones(1,size(full_pts,1))]+fund_e1(3)./depth_pts;
mapped_x=tmp1./tmp3;
mapped_y=tmp2./tmp3; 
mapX=reshape(mapped_x, C1, C2);
mapY=reshape(mapped_y, C1, C2);
max_x=max(mapped_x);
min_x=min(mapped_x);
max_y=max(mapped_y);
min_y=min(mapped_y);
offset=round([ 2 - min(1, min_x) ; 2 - min(1, min_y)]);
cw=max(ceil(max_x), size(img2,2)) + offset(1)-1;
ch=max(ceil(max_y), size(img2,1)) + offset(2)-1;

warped_img2=zeros(ch,cw,3);
warped_mask2=zeros(ch,cw);
warped_img2(offset(2):(offset(2)+sz1-1),offset(1):(offset(1)+sz2-1),:)=img2;
warped_mask2(offset(2):(offset(2)+sz1-1),offset(1):(offset(1)+sz2-1))=1;

%% find overlapped pixels
dist_X=sqrt((mapX(1:C1,2:C2)-mapX(1:C1,1:C2-1)).^2+(mapY(1:C1,2:C2)-mapY(1:C1,1:C2-1)).^2);
dist_Y=sqrt((mapX(1:C1-1,1:C2)-mapX(2:C1,1:C2)).^2+(mapY(1:C1-1,1:C2)-mapY(2:C1,1:C2)).^2);
length_X=(dist_X(1:end-1,:)+dist_X(2:end,:))./2;
length_Y=(dist_Y(:,1:end-1)+dist_Y(:,2:end))./2;
occlude_grid=(length_X>2*mean(length_X(:))) | (length_Y>2*mean(length_Y(:)));

epipolar_He=[reshape(H_inf2',1,9), fund_e2'];

%% backward image warping
warped_img1=optimal_backward_mapping(img1,img2,depth_img1,depth_img2,ch,cw,X,Y,mapX,mapY,offset,epipolar_He,double(occlude_grid),para_dist_disp);
warped_img1=reshape(warped_img1,[ch,cw,3]);
warped_mask1=imbinarize(rgb2gray(warped_img1),0);

%% panorama generation via average blending
average_panorama=imageBlending(warped_img1,warped_mask1,warped_img2,warped_mask2,'average');
[ssim_error, psnr_score] = full_reference_IQA(warped_img1, warped_img2);
quality_errors=[psnr_score, ssim_error];

end