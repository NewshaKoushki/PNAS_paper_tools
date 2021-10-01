function multi_im_show(varargin)

num_im = nargin;

ax = util.Image_Array([1 num_im], size(varargin{1}));
for k = 1:num_im
    axes(ax(1,k)); imshow(varargin{k}, []);
end

end
