function dsp=smoothdispl(displ,ni)

%displ is the original displacement field,
%n is to smooth how many peak point
if nargin<2, nd=1; else nd=ni; end

N=sqrt(size(displ,1));
xv=displ(:,1);
yv=displ(:,2);
   dxv=displ(:,3);
   dyv=displ(:,4);
   zv=sqrt(dxv.^2+dyv.^2);
x=reshape(xv,N,N);
y=reshape(yv,N,N);
dx=reshape(dxv,N,N);
dy=reshape(dyv,N,N);
z=reshape(zv,N,N);


   for i=1:nd
      mi=find(zv==max(zv));
      mi=max(mi);
      yi=floor(mi/N)+1;
      xi=mod(mi,N);
%sqrt(dxv(mi)^2+dyv(mi)^2)      
%sqrt(dx(xi,yi)^2+dy(xi,yi)^2)
      if xi>=N-1, xd=N-1; 
      elseif xi<=2, xd=2;
      else xd=xi;
      end
      
      if yi>=N-1, yd=N-1; 
      elseif yi<=2, yd=2;
      else yd=yi;
      end
      
      dxv(mi)=mean([dx(xd+1,yi),dx(xd-1,yi),dx(xi,yd+1),dx(xi,yd-1)]);
      dyv(mi)=mean([dy(xd+1,yi),dy(xd-1,yi),dy(xi,yd+1),dy(xi,yd-1)]);
      zv(mi)=sqrt(dyv(mi)^2+dxv(mi)^2);
      %invmi=[1:mi-1,mi+1:1];
      %dxv_interp=interp2(xv(invmi),yv(invmi),dxv(invmi),xv(mi),yv(mi),'cubic');
      %dyv_interp=interp2(xv(invmi),yv(invmi),dyv(invmi),xv(mi),yv(mi),'cubic');
      %dxv(mi)=dxv_interp;
      %yv(mi)=dyv_interp;
   end
   %xv=reshape(xv,N*N,1);
	%yv=reshape(yv,N*N,1);
	%dxv=reshape(dxv,N*N,1);
	%dyv=reshape(dyv,N*N,1);

   
    displ(:,3)=dxv;
    displ(:,4)=dyv;
    
    dsp=[xv,yv,dxv,dyv];

      
               