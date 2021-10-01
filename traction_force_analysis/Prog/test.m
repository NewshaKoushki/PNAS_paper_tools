% function test(xq,yq,xv,yv,in)
% figure
% 
% plot(xv,yv) % polygon
% axis equal
% 
% hold on
% plot(xq(in),yq(in),'r+') % points inside
% plot(xq(~in),yq(~in),'bo') % points outside
% hold off

%%

% for k = 1:10
%     util.write_gif('test_g.gif', 0.01, uint8(rand(100, 100, 3)*255));
% end

%%

[im_ts, npos] = lif_load(filename, pos);

