function [xrub, yrub] = Get_Cell_Boundaries(im, imBright)

ax = util.Image_Array([1 2], size(im));
axes(ax(1,1)); imshow(imBright, []); hold on
axes(ax(1,2)); imshow(im, []);

[~, xrub, yrub] = roipoly(uint8(im));

close
drawnow

%%
% set(gcf, 'keypressfcn', @key_cb, 'CloseRequestFcn', @close_cb)
% run = 1;
% function key_cb(~,e)
%     if strcmp(e.Key, 'escape')
%         run = 0;
%     end
% end
% 
% function close_cb(obj,~)
%     run = 0;
%     delete(obj);
% end
% 
% %%
% 
% while run
%     axes(ax(1,2));
%     [~, xrub, yrub] = roipoly(uint8(im));
%     if ~isempty(xrub) x = xrub; y = yrub; end % MESS
%     
%     try
%         axes(ax(1,1));
%         plot(xrub, yrub,'r.-', 'linewidth', 3);
%     end
% end
% 

end