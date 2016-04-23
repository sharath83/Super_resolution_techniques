% Create a ROI for every point observed
%stack = 'TheImages.tif';
stack = 'FinalBinary.tif';
info = imfinfo(stack);
load TheTemplates.mat;

final_points = csvread('points.csv');
n = size(final_points);
ROI_width = 10;
size1 = [(2*ROI_width+1)^2 1];
for p = 1:n(1)
    [xIdx, yIdx] = meshgrid(final_points(p,1)-ROI_width:final_points(p,1)+ROI_width, ...
                    final_points(p,2)- ROI_width:final_points(p,2)+ ROI_width);

    frame = imread(stack,final_points(p,5));
    temp = template(final_points(p,4),:,:);
    temp = reshape(temp,[20 20]);
    corr = normxcorr2(temp, frame);
    frame_sub = corr((xIdx(1,:))', yIdx(:,1));
    frame_sub(frame_sub < 0.1) = 0;
    %frame_sub = frame((yIdx(:,1)), xIdx(1,:));
    x = reshape(xIdx, size1);
    y = reshape(yIdx, size1);
    z = reshape(frame_sub,size1);
    x0 = final_points(p,1);
    y0 = final_points(p,2);
    [fit gof] = trippleGuassianFit(x,y,z,x0,y0);
    lobe_dist = sqrt((fit.b1-fit.b3)^2 + (fit.c1-fit.c3)^2);
    Xmid_point = fit.b2;
    Ymid_point = fit.c2;
    final_points(p,7:9) = [Xmid_point Ymid_point lobe_dist];
end

csvwrite('fit_points.csv', final_points);
%(a1/sqrt(2*pi).*exp(-((x-b1).^2/2)-((y-c1).^2/2))) + (a2/sqrt(2*pi).*exp(-((x-b2).^2/2)-((y-c2).^2/2)))