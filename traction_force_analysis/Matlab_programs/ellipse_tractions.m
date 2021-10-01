 % Plot ellipses from Tx and Ty (values of moment matrix). Tx represents the
% major axis and Ty the minor axis. The principal moment angle (theta),
% also obtained from moment matrix, is used to determine angle of
% orientation of the ellipse. The ellipses are spaced a set distance apart.
% Each ellipse is plotted from data obtained 5 minutes apart.
%05/11/10; Rama and Andreea

 sheet = {'Cell81'}
 
 pixelsize=0.462;

 % Plot traction ellipses
 eval(['cd C:\Singlecell\results' ]); 
% 
figure (1); hold on;
            displ = xlsread('Tractionellipse.xlsx',sheet{1}); 
%                    for i = 1:length(displ)
                    for i = 1:7
%                        fig=figure(i);                     
           scrsz = get(0,'ScreenSize');
           axis equal; 
%           figure('Position',[500 600 scrsz(3)/8 scrsz(4)/8])
           %figure('Position',[600 600 scrsz(4)/8 scrsz(4)/8])
       
                   h = ellipse(displ(i,10)/2,displ(i,11)/2,0,0); axis([-15 15 -15 15]);
                     set(h,'LineWidth',1,'Color','blue');
                   rotate(h,[0,0,1],(displ(i,12)));
       set(gca,'Color',[1 1 1],'ColorOrder',[0 0 0],'XColor',[0 0 0],'YColor',[0 0 0],'FontSize',8);
%                   set(text(0,10,strcat(num2str(displ(i,1)),'s')),'Color',[0 0 0],'FontSize',8);
                   set(text(0,10,'t00s'),'Color',[0 0 0],'FontSize',8);
                   box on
 
                     eval(['cd C:\Singlecell\traction\Cell81']); 
                   saveas(gcf, sprintf('tracfigure%d',i), 'emf')
               end;
            
              
% % % Plot cell boundaries           
% 
%   eval(['cd C:\Singlecell\displacement\Cell71']);
%     
%           displ1 = xlsread('C:\Singlecell\displacement\Cell71\Pos001_S001_t00_ch00.xlsx');
% 
%   eval(['cd C:\Singlecell\displacement\Cell71\areaMI']); 
%   
%       d2 = dir; si = 1;
%    d2b(1)=d2(1);			%d2b is part of d2 and will include [.][..] and original cell files
%    d2b(2)=d2(2);
%    for d2i = 3 : length(d2),
%       if ~isempty(findstr(d2(d2i).name,'mat')) 	% name has to contain 'mat'
%       		d2b(si)	=d2(d2i);
% 		      %name3 		= d2(d2i).name;
%           	si = si + 1;
%         end;
%    end;
% %    
% %      for j = 1:length(d2b),
% 
%           for j = 1:18
%  load(d2b(j).name)            
% 
% 
%             xrub = (xrub - mean(displ1(:,1))) * pixelsize;
%             yrub = (yrub - mean(displ1(:,2))) * pixelsize; %!!!! The -1 was added to fix discrepancy
%             % between y axis from traction program and y axis from a
%             % regular matlab plot
%             
% 
%              ln = plot(xrub, yrub,'r-'); axis([-200 200 -200 200]);
%           set(ln,'LineWidth',1.5);
%           
%           set(gca,'Color',[1 1 1],'ColorOrder',[0 0 0],'XColor',[0 0 0],'YColor',[0 0 0],'FontSize',8)
%           
% %          time = strcat(d2b(j).name(13:15),'s');
%             time = 't00s';
%           set(text(0,150,time),'Color',[0 0 0],'FontSize',8);
%            saveas(gcf, sprintf('morphfigure%d',j), 'emf')
%       end;
%   
            
           
