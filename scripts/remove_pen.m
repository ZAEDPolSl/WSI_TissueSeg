function [mask]=remove_pen(img,pen_color,thr_low, thr_high,thr_back,SE)

% set structuring element for morphology
SE = strel('disk',SE);

% choose thresholds based on color
if (strcmp(pen_color,'black'))
    mask = abs(img(:,:,1)-img(:,:,2))<=thr_low ...
        & abs(img(:,:,2)-img(:,:,3))<=thr_low;
elseif (strcmp(pen_color,'red'))
    mask = img(:,:,1) > thr_high ...
        & (img(:,:,2)-img(:,:,3)) < thr_low ...
        & (img(:,:,1)-(0.5*img(:,:,2)+0.5*img(:,:,3))) > 0;
elseif (strcmp(pen_color,'green'))
    mask = img(:,:,2) > thr_high ...
        & (img(:,:,1)-img(:,:,3)) < thr_low ...
        & (img(:,:,2)-(0.5*img(:,:,1)+0.5*img(:,:,3))) > 0;
elseif (strcmp(pen_color,'blue'))
    mask = img(:,:,3) > thr_high ...
        & (img(:,:,1)-img(:,:,2)) < thr_low ...
        & (img(:,:,3)-(0.5*img(:,:,1)+0.5*img(:,:,2))) > 0;
else
    warning('Wrong color selected');
    mask = false(size(img,1), size(img,2));
end

% remove background from mask
% mask = mask & img(:,:,1) < thr_back(1) & ...
%     img(:,:,2) < thr_back(2) & img(:,:,3) < thr_back(3);
mask = mask & ~((img(:,:,1) > thr_back(1) & img(:,:,2) > thr_back(2)) | ...
        (img(:,:,1) > thr_back(1) & img(:,:,3) > thr_back(3)) | ...
        (img(:,:,2) > thr_back(2) & img(:,:,3) > thr_back(3)));

% smooth mask
if (sum(sum(mask)))
    mask = bwareaopen(mask, size(SE.Neighborhood,1)*10);
    mask = imclose(mask, SE);    
    mask = imopen(mask, SE);
%     mask = imfill(mask,'holes');
end

% find individual green pixels and add to mask
if (strcmp(pen_color,'green'))
     mask = mask | img(:,:,2) > 250 ...
        & img(:,:,1) < 150 & img(:,:,3) < 100 ...
        & (img(:,:,2)-(0.5*img(:,:,1)+0.5*img(:,:,3))) > 0;
end

%     figure; imshow(mask)
%figure; imshow(img); hold on; visboundaries(mask2, 'Color', 'y','LineWidth',10)
