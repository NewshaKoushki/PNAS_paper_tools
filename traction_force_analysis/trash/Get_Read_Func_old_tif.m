% Outputs a Read function which reads an image from the original tif file

function [read_ts, read_fin, num_channels, num_frames_ts] = Get_Read_Func(Info, varargin)

f = @(f) (cellfun(@(x) ischar(x) && strcmp(x,f), varargin));
    function out = read_pair(name, default)
        ind = find(f(name));
        if ind, out = varargin{ind + 1};
        else out = default;
        end
    end

rewrite = any(f('rewrite'));

%% Read all basic info

% make file address generic to be readable by both pc and mac
base_addr = util.gen_addr(Info.addr);
fileName_ts = [base_addr 'TimeSeries.tif'];
fileName_fin = [base_addr 'Final.tif'];

if ~exist(fileName_ts, 'file')
    error('Couldn''t find the TS file');
end
if ~exist(fileName_fin, 'file')
    error('Couldn''t find the FIN file');
end

% Time Series
image_info_ts = imfinfo(fileName_ts);
% image_info_ts(1).StripOffsets << might be time
desc = image_info_ts(1).ImageDescription;
tk = regexp(desc, 'channels=(\d+)', 'tokens');

if isempty(tk), num_channels = 3; % SHIT
else num_channels = str2double(tk{1});
end

num_frames_ts = numel(image_info_ts)/num_channels;

% Final image
image_info_fin = imfinfo(fileName_fin);
desc = image_info_fin(1).ImageDescription;
%------------------ TODO

% If user has not allocated date field in Info, read it from file creation date
if isfield(Info, 'date'), date = Info.date;
else % Get date of the file
    fileInfo = dir(fileName_ts);
    date = fileInfo.date;
end

%% Reading function

    function im = read_im_ts(idx, ch)
        im = double(imread(fileName_ts, sub2ind([num_channels, num_frames_ts], ch, idx),  'Info', image_info_ts));
    end

    function im = read_im_fin(ch)
        im = double(imread(fileName_fin, ch,  'Info', image_info_fin));
    end

read_ts = @read_im_ts;
read_fin = @read_im_fin;

%% Store data

load('im_info.mat', 'inf');
% find im_name if exists in inf
old_im_names = fieldnames(inf);
[~,im_name, ~] = fileparts(base_addr);
im_name = ['S_' regexprep(im_name, '\W', '_')]; % make im_name valid structure field name
tf = ismember(im_name, old_im_names);
if rewrite || ~tf
    fprintf('Saving info\n');
    % Load/update info
    inf.(im_name).num_channels = num_channels;
    inf.(im_name).stiffness = Info.stiffness;
    inf.(im_name).date = date;
    % Save info
    save('im_info.mat', 'inf');
end

disp(im_name)
disp(inf.(im_name))

%% Clear im_info.mat
% load('im_info.mat', 'inf');
% inf = struct;
% save('im_info.mat', 'inf');

end