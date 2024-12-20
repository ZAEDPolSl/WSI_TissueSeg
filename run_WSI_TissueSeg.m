clc; clearvars; close all

% add folder with sub-functions
addpath(genpath('scripts'));

% choose initial segmentation method {'GaMRed','Otsu' or 'Peaks'}
method_name = 'GaMRed';

% load
data_name = 'Example1';
img = imread(['data/',data_name,'.png']);

% Calculate pixel counts
R = histcounts(img(:,:,1),-0.5:255.5);
G = histcounts(img(:,:,2),-0.5:255.5);
B = histcounts(img(:,:,3),-0.5:255.5);

% remove counts from artificial white pixels and lines
R(251:end) = 0;
G(251:end) = 0;
B(251:end) = 0;

% run initial segmentation
thr = zeros(1,3);
switch (method_name)
    case 'GaMRed'
        thr(1) = GaMRed_hist(0:255,R,2,0,5);
        thr(2) = GaMRed_hist(0:255,G,2,0,5);
        thr(3) = GaMRed_hist(0:255,B,2,0,5);
    case 'Otsu'
        thr(1) = 255*otsuthresh(R);
        thr(1) = round(thr(1) + (255-thr(1))*otsuthresh(R(thr(1):end)));
        thr(2) = 255*otsuthresh(G);
        thr(2) = round(thr(2) + (255-thr(2))*otsuthresh(G(thr(2):end)));
        thr(3) = 255*otsuthresh(B);
        thr(3) = round(thr(3) + (255-thr(3))*otsuthresh(B(thr(3):end)));
    case 'Peaks'
        thr = back_thr_peaks(0:255,R,G,B);
end

% plot thresholds
figure; subplot(3,1,1); hold on; box on;
bar(0:255, R, 'FaceColor', [0.9 0.9 0.9]);
plot([thr(1),thr(1)],ylim,':r','LineWidth',2);
xlabel('Red channel')
subplot(3,1,2); hold on; box on;
bar(0:255, G, 'FaceColor', [0.9 0.9 0.9]);
plot([thr(2),thr(2)],ylim,':r','LineWidth',2)
xlabel('Green channel')
subplot(3,1,3); hold on; box on;
bar(0:255, B, 'FaceColor', [0.9 0.9 0.9]);
plot([thr(3),thr(3)],ylim,':r','LineWidth',2)
xlabel('Blue channel')
saveas(gcf,['res/Thresholds/',data_name,'_thr.png'] )
close all;

% get regions above background and plot initial segmentation mask
mask = ~((img(:,:,1) > thr(1) & img(:,:,2) > thr(2)) | ...
    (img(:,:,1) > thr(1) & img(:,:,3) > thr(3)) | ...
    (img(:,:,2) > thr(2) & img(:,:,3) > thr(3)));
imwrite(mask, ['res/Masks/Mask_',method_name,'_',data_name,'.png'])

% artifacts removal
img_tmp = img;

%black pen marker removal
mask_tmp = remove_pen(img_tmp,'black',10,0,thr,3);
if(sum(sum(mask_tmp))); img_tmp = apply_mask(img_tmp,mask_tmp,1); end

%green pen marker removal
mask_tmp = remove_pen(img_tmp,'green',10,150,thr,3);
if(sum(sum(mask_tmp))); img_tmp = apply_mask(img_tmp,mask_tmp,1); end

% get regions above background
mask = ~((img_tmp(:,:,1) > thr(1) & img_tmp(:,:,2) > thr(2)) | ...
    (img_tmp(:,:,1) > thr(1) & img_tmp(:,:,3) > thr(3)) | ...
    (img_tmp(:,:,2) > thr(2) & img_tmp(:,:,3) > thr(3)));

% apply other post-processing steps
[mask_proc,mask_names] = post_proc(img, mask, method_name);
mask_final = mask_proc{4};
imwrite(mask_final, ['res/Masks/Mask_',method_name,'_proc_',data_name,'.png'])

% visualize mask on original tissue
figure; imshow(img)
hold on
visboundaries(mask_final, 'Color', 'y', 'LineWidth',4)
saveas(gcf,['res/Overlay_',method_name,'_proc_',data_name,'.png'])
close all;

