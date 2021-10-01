function write_gif(outfile, delay, im)
% If im is passed to this function (third argument) it will be written to
% gif file, otherwise the current figure will be written

if nargin < 3
    drawnow;
    set(gcf, 'color', 'w');
    im = frame2im(getframe);
end

if ndims(im) == 2
    dat = {im};
else
    [a, map] = rgb2ind(im, 256);
    dat = {a, map};
end


if ~exist(outfile, 'file')
    imwrite(dat{:}, outfile, 'gif', 'DelayTime', delay, 'loopcount', inf);
else
    imwrite(dat{:}, outfile, 'gif', 'DelayTime', delay, 'writemode', 'append');
end

