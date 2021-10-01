%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name of file: MaskSegment
% Version 1.0

% Purpose: to segment a cell in the image using k-mean square method

% Input variables: 
%   im: the image
% 
%   output variable: 
%   imout: the mask

% created: Nov 29, 2016 13:40 PM
% Author: Rosa Kaviani
% Revisions: ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [imout]= MaskSegment (im)

%im = 256* (max(max(im))-min(min(im)))/max(max(im))*im;
[~,mask] = util.kmeans(im,2);
mask = imfill(mask, 'holes');
se = strel('disk',2);
mask = imdilate(mask,se);
mask = imerode(mask,se);
mask = imfill(mask, 'holes');
temp = mask;
temp(mask ==1)=0;
temp(mask==2)=1;
imout = logical (temp);
end