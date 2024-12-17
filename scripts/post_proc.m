function [mask_out, mask_names] = post_proc(img, mask_in, thr_name, params)
%Post-processing of initial segmentation of histopathological images for
% tissue region finding.

if nargin < 4
    % load default parameters
    params = struct;
    params.chroma_thr = 2;
    params.disk_size = 3;
    params.area_thr = 0.01;
end
mask_out = cell(4,1);
mask_names = mask_out;

% remove grey stains with low Chroma component
img = rgb2lab(img);
tmp = sqrt(double(img(:,:,2)).^2 + double(img(:,:,3)).^2);
mask_out{1} = mask_in & tmp > params.chroma_thr;
mask_names{1} = [thr_name,'_stains'];

% fill holes in mask
mask_out{2} = imfill(mask_out{1},'holes');
mask_names{2} = [thr_name,'_holes'];

% Morphological operations to clean mask
mask_out{3} = imopen(mask_out{2}, strel('disk',params.disk_size));
mask_names{3} = [thr_name,'_open'];

% remove small area objects by size (smaller than X% of all tissue mask region)
tiss_stats = struct2table(regionprops(logical(mask_out{3}),'Area'));
thr_area = round(params.area_thr*sum(tiss_stats.Area));
mask_out{4} = bwareaopen(mask_out{3}, thr_area);
mask_names{4} = [thr_name,'_proc'];

end