    
    
folds={'Cell137'};

for fi1 = 1: length(folds)      
        cellfolder	= folds{fi1}
   
     eval([strcat('cd C:\Singlecell\displacement\',cellfolder)]); 
  
      d2 = dir; si = 1;
   d2b(1)=d2(1);			%d2b is part of d2 and will include [.][..] and original cell files
   d2b(2)=d2(2);
   for d2i = 3 : length(d2),
      if ~isempty(findstr(d2(d2i).name,'mat')) 	% name has to contain 'mat'
      		d2b(si)	=d2(d2i);
		      %name3 		= d2(d2i).name;
          	si = si + 1;
        end;
   end;
   
   for file = 1:length(d2b),
       load(d2b(file).name);
        
   
   try
    size_im = 896;
    m = size_im;
    n = size_im;
    r = zeros(m,n);
    for i = 1:m
        r(i,:) = i;
    end
    x = reshape(r, m*n,1);
    c = zeros(m,n);
    for j = 1:n
        c(:,j) = j;
    end
    
  %  area_cell = polyarea(xrub,yrub);
  %  radius = sqrt(polyarea/3.14159);
    
    
    y = reshape(c, m*n,1);
    in = inpolygon(x,y,xrub,yrub);
    x_in  = x(find(in > 0));
    y_in  = y(find(in > 0));
    x_bar = sum(x_in)/length(x_in);
    y_bar = sum(y_in)/length(y_in);
    Ixx = sum(y_in.*y_in);
    Iyy = sum(x_in.*x_in);
    Ixy = sum(x_in.*y_in);
    Ixx_bar = Ixx - (y_bar^2)*length(x_in);
    Iyy_bar = Iyy - (x_bar^2)*length(x_in);
    Ixy_bar = Ixy - x_bar*y_bar*length(x_in);
    Imatr = [Ixx_bar Ixy_bar; Ixy_bar Iyy_bar];
    [Vec I_prim] = eig(Imatr);
    th_eig = atan(Vec(1,2)/Vec(2,2));
    Imin = I_prim(1,1)
    Imax = I_prim(2,2)
    thet = th_eig*180/pi;
   POSITION(file)=str2num(d2b(file).name(4:6));
   CELLNAME(file)=str2num(d2b(file).name(18));
   SIDE (file) = str2num(d2b(file).name(11));
   TIME(file) = str2num(d2b(file).name(14:15));
   THETA(file) = thet;
   POLARITY(file)=Imax/Imin;
   catch
       continue
   end;
   end;
   
   % Write output to an excel file
   data(:,1)=POSITION;
   data(:,2)=SIDE;
   data(:,3)=TIME;
   data(:,4)=CELLNAME;
 
data(:,5)=THETA;
data(:,6)=POLARITY;

   header = '';
   colnames = {'Position','Side','Time(s)','Cell','Angle(areaMI)','Polarity(shape)'};
   filename = 'AreaMI.xlsx';
   eval(['cd C:\Singlecell\results']); 
   xlswrite(data,header,colnames,filename,cellfolder);
   close all;
   clear all;

end;