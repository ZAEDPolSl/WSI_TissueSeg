function img = apply_mask(img,mask,inv)

% invert mask if needed
if inv
    mask = ~mask;
end
    
% apply mask
mask = uint8(mask);
for c=1:size(img,3)
    tmp = img(:,:,c).*mask;
    tmp(tmp == 0) = 255;
    img(:,:,c) = tmp;
end