function ax = Image_Array(s, ims, varargin)

f = @(f) (cellfun(@(x) ischar(x) && strcmp(x,f), varargin));
    function out = read_pair(name,def)
        ind = find(f(name));
        if ind, out = varargin{ind + 1};
        else out = def;
        end
    end

r = s(1); c = s(2); % array size
% ims: image size
sp = read_pair('sp',[5 5]); % Spacing between images

set(0,'units','pixels');
ss = get(0,'screensize'); ss = ss([4 3]);
scale = min(1, min((ss - (s+1).*sp - [200 200])./(s.*ims)));

ims = floor(scale*ims);

fig_s = (s.*ims + (s+1).*sp)+[20 35];

figure('units','pixels','position', [fliplr(floor((ss - fig_s)/2))-[0 20] fliplr(fig_s)]);

ax = zeros(s);
imh = zeros(s);
for k1 = 1:c
    for k2 = 1:r
        ax(r+1-k2,k1) = axes('units','pixels','position',round([25+k1*sp(2)+(k1-1)*ims(2), k2*sp(1)+(k2-1)*ims(1), ims(2), ims(1)]));
        axis off
        imh(r+1-k2,k1) = imshow(uint8(zeros(ims)), []);
    end
end

end
%% normalized
% r = s(1); c = s(2);
%
% % figure('units','normalized','position',[.2 .4 .6 .3]);
% ax = zeros(s);
% for k1 = 1:c
%     for k2 = 1:r
%         ax(r+1-k2,k1) = axes('units','normalized','position',[.1+(k1-1)/c, .1+(k2-1)/r, 1/c, 1/r]*.8, varargin{:});
%         axis off
%     end
% end

