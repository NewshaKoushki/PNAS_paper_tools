% Compatible with Crap Wrap
imrid = @(x) double(imread(util.gen_addr(x)));

prop.pixelsize = 0.3245;
prop.young = 20000;

rect = [650 850 639 639];
crop = @(x) imcrop(x, rect);

imFB_fin = crop(imrid('..\rawdata\Cell0518\Mark_and_Find_002\05182015_TF_9kPa_Position001_t9999_ch00.tif'));

imFB_org = imrid('..\rawdata\Cell0518\Mark_and_Find_001\05182015_TF_9kPa_Position001_t0000_ch00.tif');

imFB = crop(imFB_org);
[drift_x, drift_y] = disp_on_blocks(imFB_fin, imFB, round(size(imFB_fin), 1), 0); % TODO size

imFB = imcrop(imFB_org, rect+[drift_x drift_y 0 0]);

if isempty(imFB) || ~all(size(imFB) == size(imFB_fin))
    error('displacement trim results out of bounds of image');
end

imBigFB = imcrop(imrid('..\rawdata\Cell0518\Mark_and_Find_001\05182015_TF_9kPa_Position001_t0000_ch01.tif'), rect+[drift_y drift_x 0 0]);

[x, y, dx, dy] = beads_imcorr_bart(imFB_fin, imFB);

[xrub, yrub] = Get_Cell_Boundaries(imBigFB);

puppy_code(prop, x, y, dx , dy, xrub, yrub)

