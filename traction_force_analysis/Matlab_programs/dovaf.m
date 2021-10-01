function v = dovaf(x,y)
% compute variance accounted for % between two signals
% v = vaf(x,y) where x is observed and y is predicted
%
% 22/2/94 rek
% 27/4/96 GNM changed from cov to var
% v=100*(1-cov(x-y)/cov(x));   % REK

% v= 100*(1 - sum((x - y).^2)./sum(x.^2)); % this works only with means of x (and y?) removed !
% GNM changed to std !
% 
% v0 =100*(1-cov(x-y)/cov(x))					% REK
% v1 = 100*(1 - sum((x - y).^2)./sum(x.^2))		% GNM error

% v = 100*( 1 - (std(x - y)/std(x)).^2);			% GNM fix same as REK
ss = size(x);
if min(ss) > 1,
   disp('error calculating vaf on matrices, muxt be a vector');
   return
end;
v =100*(1-cov(x-y)/cov(x));						% REK method
