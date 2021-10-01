% Outputs a Read function which reads an image from the original tif file

function [read_ts, read_fin, im_props] = Get_Read_Func(addr, varargin)

f = @(f) (cellfun(@(x) ischar(x) && strcmp(x,f), varargin));
    function out = read_pair(name, default)
        ind = find(f(name));
        if ind, out = varargin{ind + 1};
        else out = default;
        end
    end

%%

% make file address generic to be readable by both pc and mac
ts_folder = util.gen_addr([addr '/Mark_and_Find_001/']);
fin_folder = util.gen_addr([addr '/Mark_and_Find_002/']);

dr = dir(ts_folder);
init_file_name = 0;

for k = 1:length(dr)
    if dr(k).isdir, continue; end    
    tk = regexpi(dr(k).name, 'Position(\d+)_t(\d+)_ch(\d+).tif', 'tokens');
    if isempty(tk), continue; end
    stuff_idx(k,:) = str2double(tk{:}); % position, timestamp and channel
    if ~init_file_name
        tk = regexp(dr(k).name, '(\w+)_Position\w+.tif', 'tokens');
        init_file_name = tk{1}{1};
    end
end

im_props.max_pos = max(stuff_idx(:,1)); % number of positions
im_props.max_ts = max(stuff_idx(:,2)); % number of timestamps

%% Reading function

    function im = read_im_ts(pos, ts, ch)
        if pos > im_props.max_pos || pos < 1
            warning('pos %d does not exist', pos); return;
        end
        if ts > im_props.max_ts || ts < 0
            warning('ts %d does not exist', ts); return;
        end
        if ch > 3 || ch < 0
            warning('channel %d does not exist', ch); return;
        end
        im = double(imread(sprintf('%s%s_Position%03d_t%04d_ch%02d.tif', ts_folder, init_file_name, pos, ts, ch)));
    end

    function im = read_im_fin(pos, ch)
        if pos > im_props.max_pos || pos < 1
            warning('pos %d does not exist', pos); return;
        end
        if ch > 3 || ch < 0
            warning('channel %d does not exist', ch); return;
        end
        im = double(imread(sprintf('%s%s_Position%03d_t9999_ch%02d.tif', fin_folder, init_file_name, pos, ch)));
    end

read_ts = @read_im_ts;
read_fin = @read_im_fin;

end