% This program reads a binary stack of images with repeated exposure of a
% cell and finds the real valued pixel positions of molecules in it using
% Super resolution techniques

stack = 'FinalBinary.tif';
info = imfinfo(stack);
load TheTemplates.mat;

% Initializing for speed
points = zeros(60,6);
fr_points = zeros(0,6);
final_points = zeros(0,6);
u = 1; %Unique point identifier
for f = 1:numel(info) % For every frame
    Frame = imread(stack,f);
    
    for t = 1:60
        temp = template(t,:,:);
        temp1 = reshape(temp,[20 20]);
        corr = normxcorr2(temp1, Frame);
        maxcorr = max(max(corr));
        [x y] = find(corr == maxcorr,1);
        
        if t > 1 %For every template create an unique point identifier
            
            if abs(prevx-x) <= 5 & abs(prevy-y) <= 5
                u = u;
            else
                u = u + 1;
            end
%         else
%             u = 1;
        end

        points(t,:) = [x y maxcorr t f u];
        prevx = x;
        prevy = y;
    end
    %[t1] = find(points(:,3) == max(points(:,3)));
    %[t1] = find(points(:,3) > .6);
    points = points(points(:,3) > 0.65,:);
    fr_points = vertcat(fr_points, points);
end

% Filter out the final Double helix candidates
i = 1;
for k = (unique(fr_points(:,6)'))
    indx = fr_points(:,6) == k;
    max_indx = fr_points(indx,3) == max(fr_points(indx,3));
    subset = fr_points(indx,:);
    final_points = vertcat(final_points,subset(max_indx,:));
    i = i+1;
%     disp(k);
end
text = sprintf('No. of points captured is %d', i);
disp(text);

%csvwrite('points.csv', final_points);

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


            
            
%Group similar coordinates
% Cluster the final points to get Molecules' positions
% We assume that the same point repeated in multiple frames belong to the
% same molecule

coord=final_points(:,[1:2,4]);
dist=pdist(coord);
tree=linkage(dist, 'single'); %linkn based on shortest distance

%Clustering using sqrt(n/2) as the numebr clusters, where 'n' is the number
%of final points
final_points(:,10)= cluster(tree,'maxclust',round(sqrt(n(1)/2)));

% Name columns of the dataset
final_points = dataset({final_points 'Xpixel','Ypixel','Corr','Template',...
    'Frame','u','Xreal','Yreal','lobeD','Cluster'});

% Group points by Cluster Index to get final position of Molecules
% Get the Orientation of the molecule -> Template with maximum mode
group = grpstats(final_points, {'Cluster'}, {'mode'},'Datavars','Template');

% Get final position coordinates
% Assuming that mean position across all the frames as the final position
% coordinates
group2 = grpstats(final_points, {'Cluster'}, {'mean'},'Datavars',{'Xreal','Yreal'});

% Prepare final table for Image reconstruction
group(:,4:5) = group2(:,3:4);
% Assuming that point appearing in only one frame could be a Noise
group = group(group.GroupCount > 1,:);

% Image reconstruction
Im = zeros(size(frame));
n = size(group);
for p = 1:n(1)
    xPoint = round(group.Var4(p));
    yPoint = round(group.Var5(p));
    tIndex = group.mode_Template(p);
    
    Im(xPoint-9:xPoint+10,yPoint-9:yPoint+10) = template(tIndex,:,:);
    
end

imagesc(Im)
title('Reconstructued Image');





