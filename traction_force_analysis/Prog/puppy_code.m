
function [rmst_iterative,max_traction, theta0, Trace_moment, Uecm, sumforce, prestress, area_cell] = ...
    puppy_code(prop, x, y, dx , dy, xrub, yrub, ResultName, ts)

%  INPUT: PIXELSIZE, YOUNG'S MODULUS, POISSON'S RATIO

% pixelsize = 0.189; %0.92 for Rama's Inverted 10x ; 1.297 for Xinyong 10x; 0.462 for Bart - 20x
young = prop.young; % Pa
pixelsize = prop.pixelsize;

pois = 0.48;

young_e12         = young * 1e-12; % why? (Rose wonders)
RMSTmin       = 10; % don't change

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

base_gif_name = @(name) [ResultName name '.gif'];

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

% Settings for an animated movie

% mvfile1 = strcat(savedirectory,'Constraint_Traction',cellfolder,'_',num2str(1),'.avi');
% mvfile2 = strcat(savedirectory,'UnConstr_Traction',cellfolder,'_',num2str(1),'.avi');
% mvfile3 = strcat(savedirectory,'Displacement',cellfolder,'_',num2str(1),'.avi');
% aviobj1 = VideoWriter(mvfile1,'Uncompressed AVI');
% aviobj2 = VideoWriter(mvfile2,'Uncompressed AVI');
% aviobj3 = VideoWriter(mvfile3,'Uncompressed AVI');
%
% aviobj1.FrameRate = 2;
% aviobj2.FrameRate = 2;
% aviobj3.FrameRate = 2;

%==========================================================================
%  CONSTANTS

a1   = (1.0 + pois) * (1.0 - pois) / (pi * young_e12);
b1   = (1.0 + pois) * pois / (pi * young_e12);
c1   = (1.0 + pois) * pois / (pi * young_e12);

% LOOP OVER FILES


%====================================================================================
ncells = 1;
for   cellno = 1:ncells
    
    %====================================================================================
    %  CENTER DISPLACEMENTS, PIXELS TO MICRONS, SPACING
    x_org = x;
    y_org = y;
    
    x  = (x - mean(x(:))) * pixelsize;
    y  = (y - mean(y(:))) * pixelsize;
    ux = dx * pixelsize;
    uy = dy * pixelsize;
    
    xv = x(:);
    yv = y(:);
    xv_orig = xv;
    yv_orig = yv;
    
    spacing = yv(2) - yv(1);
    
    N   = sqrt(length(xv));
    
    %=====================================================================================
    % Filter for the displacement matrix
    
    ux  = myfilter_exp(ux, 15); % f0 = 15, the cut-off frequency.
    uy  = myfilter_exp(uy, 15);
    
    %====================================================================================
    %  kx, ky, AND  k_abs
    
    clear mala*;
    for i = 1 : (N/2)
        malax(i,:) = 0 : ((N/2)-1);
        malay(i,1:(N/2)) = (N/2)-i;
    end;
    
    kx = [ malax  malax-(N/2);  malax  malax-(N/2) ];
    ky = [ malay-(N/2)  malay-(N/2);  malay  malay ];
    ky = flipud(ky);
    kx(:,(N/2+1)) =  kx(:,(N/2+1));
    ky((N/2+1),:) =  ky((N/2+1),:);
    
    k_abs = sqrt(kx.^2 + ky.^2);
    
    %====================================================================================
    %  ALPHA
    
    alpha = atan2(ky, kx);
    if kx(1,1) == 0 && ky(1,1) == 0,
        alpha(1,1) = 1.57080;
    end;
    
    %====================================================================================
    %  C AND D
    
    Cx = ((k_abs * young_e12) / (2 * (1 - pois^2))) .* (1 - pois + pois .* (cos(alpha)).^2);
    Cy = ((k_abs * young_e12) / (2 * (1 - pois^2))) .* (1 - pois + pois .* (sin(alpha)).^2);
    D  = ((k_abs * young_e12) / (2 * (1 - pois^2))) .* (pois .* sin(alpha) .* cos(alpha));
    
    D(:,(N/2+1)) = zeros(N,1);
    D((N/2+1),:) = zeros(1,N);
    
    %====================================================================================
    %  CALCULATE THE TRACTIONS
    
    Dx = fft2(ux * 2 * pi / (N * spacing));
    Dy = fft2(uy * 2 * pi / (N * spacing));
    
    Tx = Cx .* Dx  +  D  .* Dy;
    Ty = D  .* Dx  +  Cy .* Dy;
    
    tx = real(ifft2(Tx));
    ty = real(ifft2(Ty));
    
    %====================================================================================
    % RMS TRACTION, THETA, TRACE, U
    
    rmst_onestep = sqrt(mean2(tx.^2 + ty.^2));
    
    % average over the large region = sum_over_region(fraction_i* force_i)/sum_over_region(force_i)
    
    clear i;
    Dxx = imag(Tx(1,2)) * (N * spacing / (2*pi)) * spacing^2 * 1e-6;
    Dyy = imag(Ty(2,1)) * (N * spacing / (2*pi)) * spacing^2 * 1e-6;
    Dxy = imag(Ty(1,2)) * (N * spacing / (2*pi)) * spacing^2 * 1e-6;
    Dyx = imag(Tx(2,1)) * (N * spacing / (2*pi)) * spacing^2 * 1e-6;
    
    % [a_sym, lambda, theta] = rotate1([Dxx Dxy; Dyx Dyy]);
    % Trace_moment = Dxx + Dyy;
    
    Uonestep = 0.5 * sum(sum([tx ty] .* [ux uy])) * spacing^2 * 1e-6;
    
    interior = find(x ~= min(x(:)) & x ~= max(x(:)) & y ~= min(y(:)) & y ~= max(y(:)));
    Uonestepint = 0.5 * sum(sum([tx(interior) ty(interior)].*[ux(interior) uy(interior)])) * spacing^2 * 1e-6;
    
    %====================================================================================
    %  Save figures parameters
    
    disptitle = strcat('Displacement(Microns)');
    tractitle = strcat('Traction map (Pa)');
    untractitle = strcat('Unconstrained FTTC: Tractions (Pa)');
    
    %====================================================================================
    %%  FIGURES
    
    set(gca,'FontName','Times','fontsize',18)
    fig1 = figure(1);
    
    rect = [100, 100, 800, 560];
    set(gcf,'Position',rect,'Renderer','zbuffer'); clf; hold on;
    
    % caxis([0 0.3]) % This the scale bar for the displacement (microns)
%     im_displacement.x = x(2:end-2,2:end-2);
%     im_displacement.y = y(2:end-2,2:end-2);
%     im_displacement.z = sqrt(ux(2:end-2,2:end-2).^2+uy(2:end-2,2:end-2).^2);
%     
    surf(x,y,sqrt(ux.^2+uy.^2));
    view(2); colormap jet; shading interp;
    
    cbh = colorbar; set(cbh,'FontSize',18,'FontWeight','bold','Color','k');
    caxis([0 3])
    m1 = max(max(sqrt(ux.^2+uy.^2)));
    h = quiver3(x,y,m1*ones(size(x)),ux,uy,zeros(size(x)),1);
    set(h,'Color','w','LineWidth',1,'Color','k');
    axis equal; axis tight; axis ij;
    set(gca, 'XColor', 'k');
    set(gca, 'YColor', 'k');
    axis([-100 100 -100 100]);
    %axis([min(x(2:end-2)) max(x(2:end-2)) min(y(2:end-2)) max(y(2:end-2))])
    %   xlabel('x (\mum)','FontSize',18,'FontWeight','bold'); ylabel('y (\mum)','FontSize',18,'FontWeight','bold');
    set(gca,'FontSize',18,'FontWeight','bold','Color','k'); box on;
    %     disptime = str2num(savefilename(14:15))*5;
    %     headerdisp = strcat(disptitle,';',savefilename(1:6),';',...
    %         num2str(disptime),'min');
    headerdisp = '';
    title( 'Displacement','FontSize',18,'FontWeight','bold','Color','k');
    export_fig(fig1,[ResultName 't' num2str(ts) '_Displacement.tif'],'-transparent');
    util.write_gif(base_gif_name('Displacement'), 0.5 + (ts == 1)*1.5);
    
    close(fig1);
    
    
    %%
    set(gca,'FontName','Times','fontsize',18)
    fig2 = figure(2); set(gcf,'Position',rect,'Renderer','zbuffer'); clf; hold on;
    %  caxis([0 400]) % This is the scale bar for unconstrained traction (Pa)
    surf(x,y,sqrt(tx.^2+ty.^2)*1e12);
    %   view(2); colormap jet; shading faceted;
    view(2); colormap jet; shading interp;
    caxis([0 500])
    cbh = colorbar; set(cbh,'FontSize',18,'FontWeight','bold','Color','k');
    m2 = max(max(sqrt(tx.^2+ty.^2)*1e12));
    h = quiver3(x,y,m2*ones(size(x)),tx,ty,zeros(size(x)),1);
    set(h,'Color','w','LineWidth',1,'Color','k');
    axis equal; axis tight; axis ij;
    set(gca, 'XColor', 'k');
    set(gca, 'YColor', 'k');
     axis([-100 100 -100 100]);
    %axis([min(x(2:end-2)) max(x(2:end-2)) min(y(2:end-2)) max(y(2:end-2))])
    %   xlabel('x (\mum)','FontSize',18,'FontWeight','bold'); ylabel('y (\mum)','FontSize',18,'FontWeight','bold');
    set(gca,'FontSize',18,'FontWeight','bold'); box on;
    title(untractitle,'FontSize',18,'FontWeight','bold', 'Color', 'k');
    export_fig(fig2,[ResultName 't' num2str(ts) '_untractitle.tif'],'-transparent')
    util.write_gif(base_gif_name('Untractitle'), 0.5 + (ts == 1)*1.5);
    close(fig2);
    %%
    %        figure(4); set(gcf,'Position',rect);
    %        clf; hold on; set(gca,'Visible','off');
    %        left = 0.2; left2 = 0.5; mov = 0.05; top = 1; difv = 0.05;
    %        text(left,top,'Pixel to \mum:','FontSize',18,'FontWeight','bold');
    %        text(left2,top,num2str(pixelsize,'%10.3f'),'FontSize',18,'FontWeight','bold');
    %        text(left,top-difv,'Young''s modulus (Pa):','FontSize',18,'FontWeight','bold');
    %        text(left2,top-difv,num2str(young*1e12,'%10.0f'),'FontSize',18,'FontWeight','bold');
    %        text(left,top-2*difv,'Poisson''s ratio:','FontSize',18,'FontWeight','bold');
    %        text(left2,top-2*difv,num2str(pois,'%10.3f'),'FontSize',18,'FontWeight','bold');
    %        left = -0.1; top = top-5*difv;
    %        h = text(left+mov,top+difv,'Unconstrained FFTC'); set(h,'FontSize',18,'FontWeight','bold');
    %        text(left,top-difv,'RMS traction (Pa):','FontSize',18,'FontWeight','bold');
    %        text(left+mov,top-2*difv,num2str(rmst_onestep*1e12,'%10.5f'),'FontSize',18,'FontWeight','bold');
    %        text(left,top-3*difv,'Orientation of principle tractions (^o):','FontSize',18,'FontWeight','bold');
    %        text(left+mov,top-4*difv,num2str(theta*180/pi,'%10.5f'),'FontSize',18,'FontWeight','bold');
    %        text(left,top-5*difv,'Net contractile moment (pJ):','FontSize',18,'FontWeight','bold');
    %        text(left+mov,top-6*difv,num2str(-Trace_moment*1e12,'%10.5f'),'FontSize',18,'FontWeight','bold');
    %        text(left,top-7*difv,'Total strain energy (pJ):','FontSize',18,'FontWeight','bold');
    %        text(left+mov,top-8*difv,num2str(Uonestepint*1e12,'%10.5f'),'FontSize',18,'FontWeight','bold');
    
    
    %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    %     CONSTRAINED FTTC
    %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
    %====================================================================================
    %  N, RESHAPE DISPLACEMENTS
    
    N = sqrt(numel(xv)) * 2;
    
    iteration = 100;
    the_same  = 1e-6;
    
    ux    = ux * 2 * pi / (N * spacing);
    uy    = uy * 2 * pi / (N * spacing);
    
    [x,y]  = meshgrid([(min(xv) - spacing*N/4 : spacing : min(xv)), x(1,2:size(x,2)-1),...
        (max(xv) : spacing : max(xv) + spacing*N/4)],...
        [(min(yv) - spacing*N/4 : spacing : min(yv)), y(2:size(y,1)-1,1)',...
        (max(yv) : spacing : max(yv) + spacing*N/4)]);
    ux = [zeros(N/4,N); zeros(N/2,N/4), ux, zeros(N/2,N/4); zeros(N/4,N)];
    uy = [zeros(N/4,N); zeros(N/2,N/4), uy, zeros(N/2,N/4); zeros(N/4,N)];
    
    xv     =  reshape(x,N^2,1);
    yv     =  reshape(y,N^2,1);
    uxv    =  reshape(ux,N^2,1);
    uyv    =  reshape(uy,N^2,1);
    
    %====================================================================================
    %  CELL BOUNDARY
    
    xrub = (xrub - mean(x_org(:))) * pixelsize;
    yrub = (yrub - mean(y_org(:))) * pixelsize;
    
    area_cell = polyarea(xrub,yrub);
    
    cell1 = find(xv >= min(xv_orig) & xv <= max(xv_orig) &...
        yv >= min(yv_orig) & yv <= max(yv_orig));
    cell  = inpolygon(xv(cell1),yv(cell1),xrub,yrub);
    cell  = cell1(cell > 0);                                % POINTS IN THE INTERIOR OF THE CELL
    
    %====================================================================================
    %  kx, ky, AND k_abs
    
    clear mala*;
    for i = 1:(N/2)
        malax(i,:) = 0:((N/2)-1);
        malay(i,1:(N/2)) = (N/2)-i;
    end;
    
    kx = [ malax  malax-(N/2);  malax  malax-(N/2) ];
    ky = [ malay-(N/2)  malay-(N/2);  malay  malay ];
    ky = flipud(ky);
    kx(:,(N/2+1)) =  kx(:,(N/2+1));
    ky((N/2+1),:) =  ky((N/2+1),:);
    
    k_abs = sqrt(kx.^2 + ky.^2);
    
    %====================================================================================
    %  ALPHA
    
    alpha = atan2(ky,kx);
    if kx(1,1) == 0 & ky(1,1) == 0,
        alpha(1,1) = 1.57080;
    end;
    
    %====================================================================================
    %  C AND D
    
    Cx = ((k_abs * young_e12) / (2 * (1 - pois^2))) .* (1 - pois + pois .* (cos(alpha)).^2);
    Cy = ((k_abs * young_e12) / (2 * (1 - pois^2))) .* (1 - pois + pois .* (sin(alpha)).^2);
    D  = ((k_abs * young_e12) / (2 * (1 - pois^2))) .* (pois .* sin(alpha) .* cos(alpha));
    
    D(:,(N/2+1)) = zeros(N,1);
    D((N/2+1),:) = zeros(1,N);
    
    %====================================================================================
    %  A AND B
    
    if k_abs(1,1) == 0,
        k_abs(1,1) = 1;
    end;
    
    Ax = a1 * ( 2 * pi ./ k_abs) + b1 * (2 * pi ./ k_abs) .* (sin(alpha)).^2;
    Ay = a1 * ( 2 * pi ./ k_abs) + b1 * (2 * pi ./ k_abs) .* (cos(alpha)).^2;
    B  = c1 * (-2 * pi ./ k_abs) .* sin(alpha) .* cos(alpha);
    
    B(:,(N/2+1)) = zeros(N,1);
    B((N/2+1),:) = zeros(1,N);
    
    %====================================================================================
    %  CALCULATE THE FIRST TRACTIONS
    
    Dx = fft2(ux);
    Dy = fft2(uy);
    
    Tx = Cx .* Dx  +  D  .* Dy;
    Ty = D  .* Dx  +  Cy .* Dy;
    
    tx = real(ifft2(Tx));
    ty = real(ifft2(Ty));
    
    %====================================================================================
    %  ITERATIONS
    
    step = 1;
    are_the_same = 0;
    
    while step <= iteration && are_the_same == 0,
        
        %====================================================================================
        %  FORM NEW (MIXED) TRACTION MATRIX
        
        tx2 = zeros(N,N);
        ty2 = zeros(N,N);
        
        tx2(cell) = tx(cell);
        ty2(cell) = ty(cell);
        
        %====================================================================================
        %  CALCULATE THE INDUCED DISPLACEMENTS
        
        Tx = fft2(tx2);
        Ty = fft2(ty2);
        
        Dx = Ax .* Tx  +  B  .* Ty;
        Dy = B  .* Tx  +  Ay .* Ty;
        
        dx2 = real(ifft2(Dx));
        dy2 = real(ifft2(Dy));
        
        %====================================================================================
        %  FORM NEW (MIXED) DISPLACEMENT MATRIX
        
        dx3 = dx2;
        dy3 = dy2;
        
        dx3(cell) = ux(cell);
        dy3(cell) = uy(cell);
        
        %====================================================================================
        %  CALCULATE THE TRACTIONS
        
        Dx = fft2(dx3);
        Dy = fft2(dy3);
        
        Tx = Cx .* Dx  +  D  .* Dy;
        Ty = D  .* Dx  +  Cy .* Dy;
        
        tx3 = real(ifft2(Tx));
        ty3 = real(ifft2(Ty));
        
        %====================================================================================
        %  FORM NEW (MIXED) TRACTION MATRIX
        
        tx = 0.5 * tx3 + 0.5 * tx2;
        ty = 0.5 * ty3 + 0.5 * ty2;
        
        if spacing <= 4,
            tx = 0.2 * tx3 + 0.8 * tx2;
            ty = 0.2 * ty3 + 0.8 * ty2;
        end;
        
        maxt(step) = max(max(tx3.^2 + ty3.^2));
        
        %====================================================================================
        % CHECK IF THE TRACTIONS HAVE CONVERGED
        
        if step > 3,
            if abs(maxt(step) - maxt(step-1)) <= the_same% * maxt(step) % Rose's MESS
                are_the_same = 1;
            end
        end
        step = step + 1;
    end
    if are_the_same ~= 1, error('Convergence was not reached!'); end;
    
    %====================================================================================
    %  SUBTRACT mean(t(cell)), BACK TO ORIGINAL ux
    
    tx = real(tx3);
    ty = real(ty3);
    tx(cell) = tx(cell) - mean(tx(cell));
    ty(cell) = ty(cell) - mean(ty(cell));
    
    
    tx4      = tx(cell);
    ty4      = ty(cell);
    
    
    
    ux = ux * N * spacing / (2 * pi);
    uy = uy * N * spacing / (2 * pi);
    
    % % Save tx(cell) and ty(cell) into a file
    %         tracfilename = strcat('C:\C2C12\traction\',cellfolder,'\tracmap_',num2str(time));
    %         fid = fopen(tracfilename,'wt');     % 'wt' means "write text"
    %         if (fid < 0)
    %             error('could not open file "tracmap_time.txt"');
    %         end;
    %         tracval = [sqrt(tx(cell).^2+ty(cell).^2)*1e12];
    %         fprintf(fid,'%10.5f\n',tracval);
    %         fclose(fid);
    
    
    %====================================================================================
    % %  RMS TRACTION, TOTAL FORCE, THETA, A_CUTS, TRACE, U
    
    
    count = 1;
    for loop = 1: size(cell,1),
        if or(abs(tx4(loop)*1e12) >=RMSTmin,abs(ty4(loop)*1e12) >=RMSTmin),
            tx5(count)        = tx4(loop);
            ty5(count)        = ty4(loop);
            count             = count + 1;
            data6(loop,:) =  {tx4(loop), ty4(loop)};
            
        end;
    end;
    
    
    median_trac = median(sqrt(tx(cell).^2 + ty(cell).^2));
    rmst_iterative = sqrt(mean2(tx(cell).^2 + ty(cell).^2));
    max_traction = max(sqrt(tx(cell).^2 + ty(cell).^2));
    % clear tx4,ty4,tx5,ty5
    
    
    tot_force   = sum(sum(sqrt(tx(cell).^2 + ty(cell).^2)))*spacing^2;
    
    Dxx = sum(sum([tx] .* [x])) * spacing^2 * 1e-6;
    Dxy = sum(sum([ty] .* [x])) * spacing^2 * 1e-6;
    Dyx = sum(sum([tx] .* [y])) * spacing^2 * 1e-6;
    Dyy = sum(sum([ty] .* [y])) * spacing^2 * 1e-6;
    
    [a_sym, lambda, theta0] = rotate1([Dxx Dxy; Dyx Dyy]);
    
    if abs(lambda(1,1)) > abs(lambda(2,2)), theta = -theta0;
    else theta = pi/4 - theta0;
    end;
    
    a_tractions = tan(theta);
    a_cuts = -1 / tan(theta);
    
    Trace_moment = Dxx + Dyy;
    
    Uecm = 0.5 * sum(sum([tx(cell) ty(cell)] .* [ux(cell) uy(cell)])) * spacing^2 * 1e-6;
    
    
    %======================================================================================
    %  CALCULATION OF PRESTRESS
    
    max_height = 6;      % height of the cell at the nucleus
    min_height = 0.5;    % height of the cell at the edges of the cell
    
    %======================================================================================
    %  POINTS OF CUTS ON THE LINE THROUGH (0,0)
    
    number = 0;
    clear point_x point_y
    
    for i = floor(-sqrt(2)*N/2) : ceil(sqrt(2)*N/2),
        
        p = [i * spacing * cos(theta), i * spacing * sin(theta)];
        axbo = max([max(xrub),max(yrub),-min(xrub),-min(yrub)]);
        if inpolygon(p(1),p(2),[-axbo axbo axbo -axbo],...
                [-axbo -axbo axbo axbo]) > 0,
            number = number + 1;
            point_x(number) = p(1);
            point_y(number) = p(2);
        end; %(if)
        
    end; %(for)
    
    %======================================================================================
    %  HEIGHTS IN POINTS
    
    height = zeros(number,1);
    middle = find(point_x == 0 & point_y == 0);
    height(middle) = max_height;
    for i = 1 : middle - 1,
        height(i) = min_height + (i-1) * (max_height - min_height) / (middle - 1);
    end;
    for i = number : (-1) : middle + 1,
        height(i) = min_height + (number-i) * (max_height - min_height) / (number - middle);
    end;
    
    %======================================================================================
    %  PRESTRESS
    
    clear sumforce diameter area_cross prestress
    
    for i = 1 : number,
        
        %===================================================================================
        %  B_CUTS, SUMFORCE
        
        b_cuts      = point_y(i) - a_cuts * point_x(i);
        if a_cuts > 0, znak = '>'; else znak = '<'; end;
        eval(['halfplane = find(y ' znak ' a_cuts * x + b_cuts);']);
        cutleft     = intersect(halfplane,cell);
        if ~isempty(cutleft) && ~isempty(setdiff(cell,cutleft))
            vec         = [sqrt(1/(1+a_tractions^2)) sqrt(1/(1+a_tractions^2))*a_tractions];
            clear vec_line;
            for j = 1 : length(cutleft),
                vec_line(j,:) = vec;
            end;
            
            sumforce(i) = sum(sum([tx(cutleft)  ty(cutleft)]  .* vec_line)) * spacing^2 * 1e12;
            %===================================================================================
            %  DIAMETER
            
            if abs(a_cuts) < 1,
                xdots    = min(xrub) : 0.1 : max(xrub);
                ydots    = a_cuts .* xdots + b_cuts;
            else
                ydots    = min(yrub) : 0.1 : max(yrub);
                xdots    = (ydots - b_cuts) ./ a_cuts;
            end
            indots      = inpolygon(xdots,ydots,xrub,yrub);
            xindots     = xdots(indots > 0);
            yindots     = ydots(indots > 0);
            jump        = 0;
            for j = 1 : length(xindots)-1,
                if abs((abs(a_cuts) < 1) * (xindots(j) - xindots(j+1)) + ...
                        (abs(a_cuts) >= 1) * (yindots(j) - yindots(j+1))) > 0.11,
                    jump = [jump, j];
                end;
            end;
            if length(jump) > 1,
                if jump(2) == 1, jump(2) = []; end;
            end;
            jump = [jump, length(xindots)+1];
            diameter(i)   = 0;
            for j = 1 : length(jump)-1,
                diameter(i) = diameter(i) + sqrt((xindots(jump(j)+1) - xindots(jump(j+1)-1))^2 +...
                    (yindots(jump(j)+1) - yindots(jump(j+1)-1))^2);
            end;
            
            %===================================================================================
            %  AREA_CROSS
            
            radius = (diameter(i)^2 / 4 + height(i)^2) / (2 * height(i));
            alpha  = 2 * atan2(diameter(i)/2, radius - height(i));
            area_cross(i) = 0.5 * radius^2 * (alpha - sin(alpha));
            
            %===================================================================================
            %  PRESTRESS
            
            prestress(i) = sumforce(i) / area_cross(i);
        else
            prestress(i) = NaN;
        end; %(if ~isempty(cutleft))
    end; %(for i)
    
    
    %%  FIGURES 3
    
    fig3 = figure(3);
    rect = [100, 100, 800, 560];
    load('MyColormaps','mycmap')
    set(gcf,'Position',rect,'Renderer','zbuffer'); clf; hold on;
    %      set(gcf,'color',[1 1 1],'Position',rect,'Renderer','zbuffer'); clf; hold on;
    set(figure(3),'Colormap',mycmap);
    set(gca,'color',[0 0 0]);
    
    caxis([0 500]) % This is the scale bar for constrained traction (Pa)
    surf(x(2:end-2,2:end-2),y(2:end-2,2:end-2),sqrt(tx(2:end-2,2:end-2).^2+ty(2:end-2,2:end-2).^2)*1e12);
    whitebg('black');
    
    view(2); colormap jet; shading interp;
    %     view(2); shading flat;
    cbh = colorbar;
    set(gca,'FontSize',18,'FontWeight','bold', 'Color', 'K'); box on;
    set(cbh,'FontSize',18,'FontWeight','bold', 'Color', 'K');
    m3 = max(max(sqrt(tx.^2+ty.^2)*1e12));
    h = plot3(xrub,yrub,m3*ones(size(xrub)),'w-'); % Plots cell boundary
    set(h,'LineWidth',2.0);
    h = quiver3(x,y,m3*ones(size(x)),tx,ty,zeros(size(x)),1);
    set(h,'Color','w','LineWidth',1,'Color', 'K');
    axis image;
    axis tight;
    axis ij;
    %     axis([-115  115  -115  115]);
     axis([-115 115 -115 115]);
    %axis([min(xv_orig) max(xv_orig) min(yv_orig) max(yv_orig)]);
    %   set (gca, 'XTickLabel',{'','','','',''},'YTickLabel',{'','','','',''});
    xlabel('x (\mum)','FontSize',18,'FontWeight','bold','Color', 'K'); ylabel('y (\mum)','FontSize',18,'FontWeight','bold','Color', 'K');
    %     tractime = str2num(savefilename(14:15))*5;
    %     headertrac = strcat(tractitle,';',savefilename(1:6),';',...
    %         num2str(tractime),'min');
    title('traction (Pa)','FontSize',18,'FontWeight','bold','Color', 'K'); % headertrac
    export_fig(fig3,[ResultName 't' num2str(ts) '_Traction.tif'],'-transparent')    %  Save the constrained traction figure to a frame
    util.write_gif(base_gif_name('Traction'), 0.5 + (ts == 1)*1.5);
    
    close(fig3);
    %         open(aviobj1);
    %         frame = getframe(figure(3));
    %         writeVideo(aviobj1,frame);
    %%
    %figure(1);
    % h = plot3(xrub,yrub,m1*ones(size(xrub)),'w-'); set(h,'LineWidth',2.0);
    
    %         open(aviobj3);
    %         frame1 = getframe(figure(1));
    %         writeVideo(aviobj3,frame1);
    %%
    %figure(2);
    %h = plot3(xrub,yrub,m2*ones(size(xrub)),'w-'); set(h,'LineWidth',2.0);% Plots cell boundary
    %         open(aviobj2);
    %         frame2 = getframe(figure(2));
    %         writeVideo(aviobj2,frame2);
    %%
    %        figure(4);
    %        left = 0.5;
    %        h = text(left+mov,top+difv,'Constrained FFTC'); set(h,'FontSize',18,'FontWeight','bold');
    %        text(left,top-difv,'RMS traction (Pa):','FontSize',18,'FontWeight','bold');
    %        text(left+mov,top-2*difv,num2str(rmst_iterative*1e12,'%10.5f'),'FontSize',18,'FontWeight','bold');
    %        text(left,top-3*difv,'Orientation of principle tractions (^o):','FontSize',18,'FontWeight','bold');
    %        text(left+mov,top-4*difv,num2str(theta0*180/pi,'%10.5f'),'FontSize',18,'FontWeight','bold');
    %        text(left,top-5*difv,'Net contractile moment (pNm):','FontSize',18,'FontWeight','bold');
    %        text(left+mov,top-6*difv,num2str(-Trace_moment*1e12,'%10.5f'),'FontSize',18,'FontWeight','bold');
    %        text(left,top-7*difv,'Total strain energy (pJ):','FontSize',18,'FontWeight','bold');
    %        text(left+mov,top-8*difv,num2str(Uecm*1e12,'%10.5f'),'FontSize',18,'FontWeight','bold');
    %        text(left,top-9*difv,'Max cumulative force (nN):','FontSize',18,'FontWeight','bold');
    %        text(left+mov,top-10*difv,num2str(max(sumforce)*1e-3,'%10.5f'),'FontSize',18,'FontWeight','bold');
    %        text(left,top-11*difv,'Prestress (Pa):','FontSize',18,'FontWeight','bold');
    %        text(left+mov,top-12*difv,num2str(stnanmean(prestress),'%10.5f'),'FontSize',18,'FontWeight','bold');
    %        text(left,top-13*difv,'Projected area of the cell (\mum^2):','FontSize',18,'FontWeight','bold');
    %        text(left+mov,top-14*difv,num2str(area_cell,'%10.5f'),'FontSize',18,'FontWeight','bold');
    %
    %
    %====================================================================================
    %  Save program outputs
    
    %     Rose's MESS: search for savefilename and fix it
    %     Pos      = str2num(savefilename(26:28)); % Rose's MESS
    %     time     = str2num(savefilename(31:32)); % Rose's MESS
    Pos = 0;
    time = 0;
    
    %cellname = num2str(cellno);
    %key = str2num(strcat(savefilename(4:6),savefilename(9:11),savefilename(14:15),num2str(cellname)));
    
    % Entries in traction file are: 1. RMS traction, 2. Median traction, 3. CM, 4. Total force, and 5. Cell area
    %tdata          = [Pos, time, cellno, rmst_iterative*1e12,-Trace_moment*1e12,area_cell,-lambda(1,1)*1e12,-lambda(2,2)*1e12,theta*180/pi, Uecm*1e12];
    
    %====================================================================================
    
    %cd('/Singlecell/Traction');
    %filenametrac=strcat('/HTSBart/Traction/',cellfolder ,'/Traction_',savefilename(1:6),'_',num2str(cellno),'.dat');
    %save(filenametrac,'tdata','-ASCII');
    
end; %(cellno)


% close(aviobj1);
% close(aviobj2);
% close(aviobj3);

end

function [a_sym,lambda,th] = rotate1(a)

a_sym  = 0.5 * (a + a');

th = atan2((a(1,2) + a(2,1)) , (a(2,2) - a(1,1))) / 2;

U  = [cos(th) sin(th);
    -sin(th) cos(th)];

lambda = U' * a_sym * U;

end