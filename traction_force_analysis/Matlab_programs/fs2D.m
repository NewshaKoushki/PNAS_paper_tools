%image shift using 2D Fourier transform
function rms=fs2D(T,ayi,axi,shift,fta,ima,imb);

%T = length(ima);
%Nyquist = (T/2+1);

%for x = 1:T,
%    for y = 1:T,
%        if (x<=Nyquist) & (y<=Nyquist),
%            shift(y,x) = exp(-i*2*pi/T*(ax*(x-1)+ay*(y-1)));
%        elseif (x>Nyquist) & (y<=Nyquist),
%            shift(y,x) = exp(-i*2*pi/T*(-ax*(T-(x-1))+ay*(y-1)));
%        elseif (x<=Nyquist) & (y>Nyquist),
%            shift(y,x) = exp(-i*2*pi/T*(-ay*(T-(y-1))+ax*(x-1))); 
%        else
%            shift(y,x) = exp(-i*2*pi/T*(-ax*(T-(x-1))-ay*(T-(y-1))));     
%        end;
%    end;
%end;
%fta = fft2(ima);

shifta=zeros(T,T);
shifta(:,:)=shift(ayi,axi,:,:);
ftas = shifta.*fta;
imas = ifft2(ftas);
%imrms(ayi,axi) = sqrt(sum(sum((real(imas(3:T-2,3:T-2))-imb).^2)));
rms = sqrt(mean(mean((real(imas(3:T-2,3:T-2))-imb).^2)));