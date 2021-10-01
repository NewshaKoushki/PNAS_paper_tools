function change_gif_delay(first_delay, delay, gif_file)

if nargin < 3
    [f, pathname] = uigetfile({'*.gif','(*.gif)'}, 'Select GIF file');
    if ~f, warning('No file selected.'); return; end
    gif_file = [pathname f];
end

[a,  map] = imread(gif_file, 'gif', 'Frames', 'all');
imwrite(a(:, :, 1, 1), map, gif_file, 'DelayTime', first_delay, 'LoopCount', inf)

for f = 2:size(a, 4)
    imwrite(a(:, :, 1, f), map, gif_file, 'DelayTime', delay, 'WriteMode', 'append')
end

fprintf('Done\n');

