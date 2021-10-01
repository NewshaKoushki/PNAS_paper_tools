%omtcana_d
%For analysing beads displacement files from optical magnetic twisting cytometry
%GNM Setp. 1999 Revised 8 Dec. 1999 for new microscope stage (no more metal)
%this version (omtcana_d3) should be used with the Leica microscope in the Butler lab


clear all;
figures = [1 1 1 1 1 1 1];

%% New parameters for the new machine 
%size image usually 600*800
size_image_x = 800;
size_image_y = 600;

% Pixel nm correspondance
nmPerPixel = 327;               % Pete's machine: 10x: 652, 20x: 327, 40x: 162

% Gauss per A correspondance
G_per_amp = 33.8983;

% Bead constant:
C_bead = 3.332;

% Reconstruct the time course and valulate vaf between data and linear model fit from G' and G'' with T
do.RefBead = 0;
do.Ben = 0;		                % doBen: nomalize results for f = 1 Hz or maximum torque
do.Benswitch = 0;	            % parameter that is switched to 1 if doing a drug file to normalize the beads to themselves.
do.Stupid = 0;	                % correct for self inductance of the optical stage


% Parameter to select ''Good beads"
% rejection will be: if e.type = e.freq, if VAF < 90 for freq just < 100 reject
% rejection will be: if e.type = e.torq, if VAF < 90 for torq just > 20  reject
Lim.X_moved_low_nm = -10;        % 10 nm eliminate all beads with amplitudes below this number
					    		% below Lim.f_high_Hz (<20 Hz)
Lim.vaf_low = -90000;				% 90 eliminate all beads with vaf below this number
Lim.Harmonic_frac = -0.18;       % -0.18 square wave has a harmonic fraction of 0.3512, so 0.18 is halfway between sinusoid and square 
% eliminate some beads (see next line) which have a harmonic fraction greater than this number
Lim.X_high_harm_frac_nm = 200000;	        % 200nm Only apply the harmonic fraction limit for beads that move with amplitude grater than this number
Lim.X_3rd_harmonic_low_nm = 100000;        % 100nm eliminate all beads with a 2nd harmonic greater than this amount of nanometers
Lim.X_irreproducible_factor = 2000;	    % 2 If amplitude changes on the decreaseing freqeuncies by more than a factor of 2 at any frequency <=10 Hz, eliminate it.
xy_criterion = 0;	                    % 1: stdev of x-displacement must be larger than sdev of y-displacement


% Other parameters
pick.coding = 0;
e.simu = 0;
torque_sign = -1;		                % -1 for normal, -1 for inverted (if wrong E is negative).
points_per_cycle = 16;                  
tdrug = 55; 			                % limtis exclude from tdrug - 10 to tdrug + 20
tshutter = 73e-6;                       % sec
stdI_max = 1e-2;                        %for automatic detection of experiments with changing torque; default value is 1e-2



if (~pick.coding) & (~e.simu), 
   if ~exist('do_allfiles'),
      start_dir = pwd;
      [fileName.raw,directory] = uigetfiles('*.txt','Select the data file for analysis');
      if isempty(fileName.raw), return; else end;
      cd(directory);
      [fileName.image,directory] = uigetfiles('*.bmp','Select the image file');
      if isempty(fileName.image) return; else end;
      fileName.raw = char(fileName.raw);
		fileName.image = char(fileName.image);
      temp = fileName.image;
      if temp(length(temp)-2:length(temp)) ~= 'bmp', return; else end;
      if figures(1),
         f1 = figure(1); set(f1,'Position',[1,400,450,300]);clf;
         set(f1,'numbertitle','off');
         set(f1,'name',fileName.image);
         colormap(gray);
         [image, map] = imread(fileName.image,'bmp');
         imagesc(image); title(fileName.image);
         drawnow;
      end;
   else
      if figures(1),
         f1 = figure(1); set(f1,'Position',[1,400,450,300]); clf;
         set(f1,'numbertitle','off');
         set(f1,'name',fileName.raw);
      end;
   end;
	temp = dir;      
   for k = 1:length(temp),
      testfile = strcmp(fileName.raw,temp(k).name);
      if testfile,
         thisfile = k;
      end;
   end;
   [fileName.CreateDate fileName.CreateTime] = strtok(temp(thisfile).date,' ');
   fileName.bytes = temp(thisfile).bytes;
   


   set(f1,'numbertitle','off');
   set(f1,'name',fileName.raw);
   %m.m = matlab11_dlmread(fileName.raw,';');
   m.m = matlab11_dlmread(fileName.raw,';');
   disp('________________________________');
   disp(['read file: ' fileName.raw]);  
   disp(['Bead constant = ' num2str(C_bead)]);
  
  	
elseif e.simu,
   f1 = figure(1); clf;
   C_bead = 3.49;
   nmPerPixel = 1;
   simu.I = 0 + j*1;
   simu.E = 5;     % 0.01 amp/nm
   simu.Gpp = 2.5;
   nb = 103;
   fosc = 10;
   dur = nb/fosc;
   t = [0:nb*16-1]'/(nb*16)*dur;
   m.m(:,1) = t;
   m.m(:,2) = log10(fosc)*ones(size(t));
   m.m(:,3) =  real(simu.I)*cos(2*pi*fosc*t) + imag(simu.I)*sin(2*pi*fosc*t); 
   % G = T/X,  X = T/G;
   simu.T = torque_sign*simu.I*G_per_amp*C_bead
   simu.X = simu.T/[simu.E + j*simu.Gpp];
   m.m(:,4) = 100 + cos(2*pi*fosc*t)*real(simu.X) + sin(2*pi*fosc*t)*imag(simu.X);  
   m.m(:,5) = 100 + 0.01*cos(2*pi*fosc*t);
   fileName.raw = 'simu';
   warning off
end;
titleStr = fileName.raw;
tic;

e.type = 0;		% type of experiment;
e.freq = 1;
e.torq = 2;
e.drug = 3;

[m.r m.c] = size(m.m);

t = m.m(:,1); 	% time;
dur = max(t);
logf = m.m(:,2); % logf;
i = -m.m(:,3); 	% current;
x = [];
y = [];
for k = 4:3:m.c-1,
   x = [x k];
   y = [y k+1];
end;
disp('file gets converted in nm');
x = m.m(:,x).*nmPerPixel;
y = m.m(:,y).*nmPerPixel;
clear m; %Ben: clear m to conserve memory
disp('file is converted in nm');

beads.totalN = size(x)*[0;1];
L = length(t);
kept1 = find(x(L,:)>0);		% last point if zero, lost the bead during imaging
kept2 = find(y(L,:)>0);
beads.kept = intersect(kept1,kept2);
if isempty(beads.kept),
   disp('no beads, no analysis done');
   disp('');
   return;
end

xo = x(17,:);			% start positions
yo = y(17,:);
x = x - ones(L,1)*xo;
y = y - ones(L,1)*yo;
disp('start position  substracted from all beads');


if do.RefBead,
   disp('high pass filter x-coordinates');  
	filt.N = 48;  
	[filt.fil,filt.a] = fir1(filt.N,1/18,'high');		% the filter
	xf = filtfilt(filt.fil,filt.a,x);
	% remove the first and last 16 points
	xf = xf(17:L-16,:);
	 
	disp('now find reference bead');
	x_sdb = std(xf);
	bead_ref = find(x_sdb == min(x_sdb));
	% reference bead
	for k=1:beads.totalN,
   	x(:,k) = x(:,k) -x(:,bead_ref);
	end;
   beads.kept = setdiff(beads.kept,bead_ref);
end;


pick.do_filt = 1;
if pick.do_filt, 
   filt.N = 48; 
   if filt.N > length(x)/3, 
      disp(['file too short to filter, not analyzing it']); disp('');
      return
   end
	[filt.fil,filt.a] = fir1(filt.N,1/18,'high');		% the filter
	xf = filtfilt(filt.fil,filt.a,x);
	yf  = filtfilt(filt.fil,filt.a,y);
    % i  = filtfilt(filt.fil,filt.a,i);   %caution: this is only good for Ben's superposition experiment!!!!!!!!!!!!!!!!!!!!
end;

x = x(:,beads.kept);
y = y(:,beads.kept);
xf = xf(:,beads.kept);
yf = yf(:,beads.kept);
xo = xo(:,beads.kept);
yo = yo(:,beads.kept);
beads.keptN = length(beads.kept);
beads.lostN = beads.totalN - beads.keptN;

% remove the first and last 16 points
logf = logf(17:L-16);
x = x(17:L-16,:);
y = y(17:L-16,:);
xf = xf(17:L-16,:);
yf = yf(17:L-16,:);
i = i(17:L-16);
t = t(17:L-16);
disp('file filtered');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% decide if independent variable is the frequency %%%%%%%%
L = length(t);
dlogf = logf(2:L)-logf(1:L-1);	% differentiate logf
b.fi = find(abs(dlogf)>0.05);  	% can only discriminate changes in logf > 0.05
b.fi = b.fi + 1;						% have index for each change in frequency, excluding 1 and logf(L);

% also check if Torque amplitude is constant or changing
temp = floor(length(i)/16);
temp = reshape(i(1:16*temp),16,temp);
temp = fft(temp)/length(temp)*2;
stdI = std(abs(temp(2,:)));


if length(b.fi) ~= 0,
   e.type = e.freq;				
   b.fi = [1; b.fi; length(t)+1]; % include the first point and last point
   % The independent variable is the frequency
   nb = length(b.fi)-1;			% number of blocks
   disp(['Varying frequency experiment, number of frequencies = ' num2str(nb)]);
elseif stdI > stdI_max, 		% torques are changing
   e.type = e.torq;
    % it's a changing Torque file, (fixed freqeuncy);
    fosc = 10^mean(logf);
    ncycles = floor(L/16) - 1;
    sincycle = sin(2*pi*[0:(ncycles+1)*16-1]/16)';  % only use cos for i
    Isum = cumsum(sincycle.*i(1:length(sincycle)))/16;
    Iindices = [1:ncycles]*16;
    Is = Isum(Iindices);
    Isnew = [Is(1);  Is(2:length(Is))-Is(1:length(Is)-1)];
	 dIs = Isnew(2:length(Isnew))-Isnew(1:length(Isnew)-1);    
    b.fi = find((abs(dIs) > 0.005));
    b.fi = b.fi*16 + 1;
    b.fi = [1; b.fi; length(t)+1]; % include the first point and last point
    nb = length(b.fi) - 1;
%   nb = floor(L/16) - 1;
   % nb = 98;     % currently this many Torques [1:98];
%   b.fi = [0:points_per_cycle:(nb)*points_per_cycle] + 1;
   % The torque changes every cycle, but may not do that in future
   % since every cycle is 16 points long, look every 4 + 16*ci, where ci is the cycle number
   disp(['Varying Torque experiment, number of Torques = ' num2str(nb)]);
else
   e.type = e.drug
   nb = floor(length(i)/points_per_cycle);
   b.fi = [0:points_per_cycle:(nb)*points_per_cycle] + 1;
   fosc = 10^mean(logf);
   % it's a constant torque, constant frequency experiment, probably drug delivery
end;
clear kept1 kept2 dlogf  

L = b.fi(nb+1)-1;
if L > length(t), disp([' truncated file: exiting ']); return;  end;

t = t(1:L); i = i(1:L);
x = x(1:L,:); y = y(1:L,:);
xf = xf(1:L,:); 
yf = yf(1:L,:);
logf = logf(1:L);

% analyze each block one at a time
p.I = zeros(1,nb);
p.T = zeros(1,nb);
p.X = zeros(nb,beads.keptN);
p.Y = zeros(nb,beads.keptN);


disp('fourier-analysis begins');
for k = 1:nb,				 % skips first cycle and last logf value
   a1 = b.fi(k);
   a2 = b.fi(k+1) - 1;
   b.t = t(a1:a2);				% b for block
   b.dur = t(a2)-t(a1);
   b.fs = 1/(b.t(2)-b.t(1));
   b.L = length(b.t);
   b.f = [0:b.L-1]*b.fs/(b.L-1);
   b.fosc = 10^mean(logf(a1:a2));
   b.i = i(a1:a2);
   b.x = xf(a1:a2,:);
   % b.y = y(a1:a2,:);
   b.I = fft(b.i)/b.L*2;
   b.t = t(round(a1 + (a2-a1)/2));
   % fill parameter arrays
	p.fss(k) = b.fs;   
   p.Cycles(k) = round(b.dur*b.fosc);
   p.durs(k) = b.dur;
   p.points(k) = b.L;    % points_per_condition
   if e.type == e.freq,  % Ben: take fundamental by taking min (first spike) of fourier
      p.Cycles_rec(k) = min(find(abs(b.I(2:length(b.I)/2)) > 0.01));
      p.foscs(k) = b.fosc;
   else
      p.Cycles_rec(k) = 1;	 % one cycle for each torque condition or drug
      p.foscs(k) = fosc;
   end;      
   % I know how many points, but now, how many cycles ?
   % Currently the setting is 2 cycles for fosc < 0.2 Hz, 5 cycles for large foscs
   % but I could get the number of cycles by finding the peak in the fft of b.I
   
	%if e.type == e.freq,
      angle_i = 2*pi*[0:p.points(k)-1]/points_per_cycle;
      anglex = angle_i + 2*pi*tshutter*p.foscs(k);
      sincyclei = sin(angle_i)';
      coscyclei = cos(angle_i)';
      sincyclex = sin(anglex)'*ones(1,beads.keptN);
   	coscyclex = cos(anglex)'*ones(1,beads.keptN);
      sincyclex3 = sin(3*anglex)'*ones(1,beads.keptN);
	   coscyclex3 = cos(3*anglex)'*ones(1,beads.keptN);
	%end;

   % I don't actually need the cumsum, it might be faster just to use sum
   % then b.isum is a single value, not an array.
   b.isum = 2*((cumsum(coscyclei.*b.i) + j*cumsum(sincyclei.*b.i)));
   b.xsum = 2*((cumsum(coscyclex.*b.x) + j*cumsum(sincyclex.*b.x)));
   %b.ysum = 2*((cumsum(coscyclex.*b.y) + j*cumsum(sincyclex.*b.y)));
   b.xsum3 = 2*((cumsum(coscyclex3.*b.x) + j*cumsum(sincyclex3.*b.x)));   
   
   % This is faster: extract the fundamental at the end of the 2 or 5 cycle condition and divide by the full length of the multiple cycles
   % kk = points_per_cycle;  error
   kk = b.L;
   p.t(k) = b.t;
   p.I(k) = b.isum(kk)./kk;   
   if do.Stupid
      i_phase_bef = angle(p.I(k));
      i_amp_bef = abs(p.I(k));
      % omega = 2*pi*p.foscs(k)
      i_phase_corr(k) = 1.054/20-0.05128*log10(p.foscs(k));
      if i_phase_corr(k) > 0.0 
         i_phase_corr(k) = 0.0;
      end;
      i_amp_corr(k) = 1.054-0.05128*log10(p.foscs(k));
      if i_amp_corr(k) > 1.0 
         i_amp_corr(k) = 1.0;
      end;
      
      i_phase_aft = i_phase_bef + i_phase_corr(k);
      i_amp_aft = i_amp_bef*i_amp_corr(k);
      
      p.I(k) = i_amp_aft*(cos(i_phase_aft) + sqrt(-1)*sin(i_phase_aft));
   end;   
   p.T(k) = torque_sign*p.I(k)*C_bead*G_per_amp*0.6;
   %%%% added 18 September
   if sqrt(p.T(k).*conj(p.T(k))) < 5e-4,
       p.T(k) = 5e-4;
   end;
   p.X(k,:) = b.xsum(kk,:)./kk;
   %p.Y(k,:) = b.ysum(kk,:)./kk;
	% Extract the harmonic at the end of the 2 or 5 cycle condition
   p.X3(k,:) = b.xsum3(kk,:)./kk;

   p.SRate(k,:) = abs(p.X(k,:))*p.foscs(k);
   % OK now we need to calculate the E, R, G", and eta, the motion is in the x_direction
   p.Gp(k,:) = real((p.T(k)*ones(1,beads.keptN))./p.X(k,:));
   p.Gpp(k,:) = -imag((p.T(k)*ones(1,beads.keptN))./p.X(k,:));
   p.eta(k,:) = p.Gpp(k,:)./p.Gp(k,:);
   p.R(k,:) = p.Gpp(k,:)./p.foscs(k);

end;
clear den kk;
disp('fourier-analysis done');


% OK lets sort the beads
if e.type == e.freq, 		
   % [temp fi] = min(p.foscs(1:floor(nb/2)));  	% Geoff: sort based on amplitude at the lowest frequency
   fi = find(p.foscs>0.5);	 % Ben: sort based on amplitude at ~ 1 Hz 
   fi = fi(1);
   [temp is] = sort(abs(p.X(fi,:)));
elseif e.type == e.torq,		% sort based on amplitude at the maximum torque
   [temp ti] = max(abs(p.T));
   [temp is] = sort(abs(p.X(ti,:)));
elseif e.type == e.drug,		% sort based on amplitude before drug delivery (at 60 seconds);
   [temp tdrugi] = max(find(t < tdrug));
   ti = floor(tdrugi/points_per_cycle);
   [temp is] = sort(abs(p.X(ti,:)));
else
   disp('Invalid experiement type, exiting'); return;
end;
disp('sorting based on amplitude (into a temp variable) done');

% now sort all intermediate results based on "is"
p.X = p.X(:,is);
% p.X2 = p.X2(:,is);
p.X3 = p.X3(:,is);
% p.FT = p.FT(:,:,is);
p.Y = p.Y(:,is);
p.SRate = p.SRate(:,is);
p.Gp = p.Gp(:,is);
p.Gpp = p.Gpp(:,is);
p.eta = p.eta(:,is);
p.R = p.R(:,is);

x = x(:,is);
y = y(:,is);
xf = xf(:,is);
yf = yf(:,is);
yo = yo(is);
xo = xo(is);
disp('all beads sorted'); 

% Now let's calculate the model bead motion based on the determined G, for each bead


% generate time course for each bead (without transients). G = T/X, X = G*T
% x = p.I*c_bead*G_per_amp*(E*sin(2*pi*fosc*t) + Gpp*cos(2*pi*fosc*t))
% where t goes from 
% x_mod = x*c_bead*G_per_amp
prev = 1;
j = (-1)^0.5;
tcount1 = clock;
for k = 1:nb,
   p.X_mod(k,:) = (p.T(k)*ones(1,beads.keptN))./[p.Gp(k,:) - j*p.Gpp(k,:)];
   sincyclex = sin(2*pi*[0:p.points(k)-1]/points_per_cycle)';
   coscyclex = cos(2*pi*[0:p.points(k)-1]/points_per_cycle)';
   cycle = [prev:prev + p.points(k)-1];
   prev = prev + p.points(k);
   % check signs
   x_mod(cycle,:) = coscyclex*real(p.X_mod(k,:)) + sincyclex*imag(p.X_mod(k,:));
   
   % now a more efficient method to compute dovaf
   x_diff =xf(cycle,:)-x_mod(cycle,:);
   cov_xf = diag(cov(xf(cycle,:)));
   cov_diff = diag(cov(x_diff));
   p.vaf(k,:)=100*(1-cov_diff./cov_xf);
   %old method, very slow
   %for kk = 1:beads.keptN,
   %   p.vaf(k,kk) = dovaf(xf(cycle,kk),x_mod(cycle,kk));
      %zz(k,kk) = dovaf(xf(cycle,kk),x_mod(cycle,kk));
   %end;
   tcount2 = clock;
   if etime(tcount2,tcount1)>1;
      disp(k/nb*100);
      tcount1 = tcount2;
   end;
   % ssr(k,:) = dossr(xf(cycle,:),x_mod(cycle,:));
end;
disp('model for bead movement computed'); 

% now calculate mean vaf for freqs up to 100 Hz
% and T > 50;
clear sincycle coscycle sincyclex coscyclex cycle sincyclex3 coscyclex3 prev

% examine which beads are below the resolution, (did not move);	
beads.moved = find(mean(abs(p.X)) > Lim.X_moved_low_nm);  % across all conditions
beads.OK_reproducible = -1;
if e.type == e.freq | e.type == e.drug,
  Lim.f_high_Hz = 20; %ignore frequencies > 10 Hz because beads usually misbehave
  Lim.f_low_Hz = 0.001; %ignore frequencies <0.001 Hz because beads usually misbehave      
  ind = find(abs(p.foscs) < Lim.f_high_Hz);
  ind1 = find(abs(p.foscs) > Lim.f_low_Hz);
  ind = intersect(ind,ind1);
  if e.type == e.drug,
     ind2 = find(p.t < tdrug - 5); %do not apply vaf during drug delivery
     ind3 = find(p.t > tdrug + 20);
     %ind = ind2;   %do not apply vaf after drug delivery
     ind = [ind2 ind3];  %do not apply vaf during drug delivery
  end;
  Q = mean(abs(p.X3(ind,:))./abs(p.X(ind,:)));  % average 2nd harmonic fraction across Torques
  i1bead = find(Q < Lim.Harmonic_frac);         % accept
  Q = mean(abs(p.X(ind,:)));
  i2bead = find(Q < Lim.X_3rd_harmonic_low_nm);    % I only want to apply this limit on beads with large amplitude X
  [beads.OK_harmonic_frac] = union(i1bead,i2bead);
  [beads.OK_harmonic_frac] = intersect(beads.moved,beads.OK_harmonic_frac);
  
  %2nd harmonic amplitude test      
  i1bead = find(max(abs(p.X3(ind,:))) < Lim.X_high_harm_frac_nm);	% across all conditions
  beads.OK_3rd_harmonic = intersect(beads.OK_harmonic_frac,i1bead);
  
  %variance accounted for test
  % ind = find(abs(p.T) > lim.T_high));  % calculated above
  beads.OK_vaf = find(mean(p.vaf(ind,:)) > Lim.vaf_low);
  beads.OK_vaf = intersect(beads.OK_3rd_harmonic,beads.OK_vaf);
  
  beads.OK_reproducible = beads.OK_vaf;
  % reproducibility test
  if e.type == e.freq
     [high_f max_fi] = max(p.foscs);
     up_fi = 1:max_fi;
     down_fi = max_fi+1:length(p.foscs);
     up_f = fix(log10(p.foscs(up_fi))*100)/100;
     down_f = fix(log10(p.foscs(down_fi))*100)/100;
     [common_f common_fi_up common_fi_down] = intersect(up_f,down_f);
     ind = find(10.^common_f < Lim.f_high_Hz);		% restrict test to below this frequency
     common_fi_up = common_fi_up(ind);
     common_fi_down = common_fi_down(ind);
     
     % thus common_f = up_f(common_fi_up);
     % and common_f = down_f(common_fi_down);
     f_diff = abs(p.foscs(common_fi_up)) - abs(p.foscs(down_fi(common_fi_down)));
     % f_oscs_for_test = abs(p.foscs(common_fi_up))
     amp_change = abs(p.X(common_fi_up,:))./abs(p.X(down_fi(common_fi_down),:));
     % semilogx(p.foscs(common_fi_up),amp_change);
     [s1 s2] = size(amp_change);
     if s1 > 1
        beads.OK_reproducible = find(max(amp_change) < Lim.X_irreproducible_factor);
     else
        beads.OK_reproducible = find(amp_change < Lim.X_irreproducible_factor);
		end
     % semilogx(p.foscs(common_fi_up),amp_change(:,beads.OK_reproducible));
  end;
  beads.OK_reproducible = intersect(beads.OK_reproducible,beads.OK_vaf);
  beads.OK = beads.OK_reproducible;

else % e.type == e.torq,
  % harmonic test:
  Lim.T_high_dynpercm2 = max(abs(p.T))/5;%   analyze only torques higher than 1/5 max torque
  ind = find(abs(p.T) > Lim.T_high_dynpercm2);
  Q = mean(abs(p.X3(ind,:))./abs(p.X(ind,:)));  % average across Torques
  i1bead = find(Q < Lim.Harmonic_frac);  % accept
  Q = mean(abs(p.X(ind,:)));
  i2bead = find(Q < Lim.X_3rd_harmonic_low_nm);    % I only want to apply this limit on beads with large amplitude X      
  [beads.OK_harmonic_frac] =  union(i1bead,i2bead);
  [beads.OK_harmonic_frac] =  intersect(beads.moved,beads.OK_harmonic_frac);
      
  % 2nd harmonic amplitude test
  i1bead = find(max(abs(p.X3)) < Lim.X_high_harm_frac_nm);	% across all conditions
  beads.OK_3rd_harmonic = intersect(beads.OK_harmonic_frac,i1bead);      
  
  %variance accounted for test
  % ind = find(abs(p.T) > lim.T_high));  % calculated above
  beads.OK_vaf = find(mean(p.vaf(ind,:)) > Lim.vaf_low);
  beads.OK_vaf = intersect(beads.OK_3rd_harmonic,beads.OK_vaf);
  beads.OK = beads.OK_vaf;
end;

if xy_criterion == 1,
	disp('exclude beads that move too much in the y-direction or in the altogether wrong direction');
	x_sdb = std(xf(ind,:));
	y_sdb = std(yf(ind,:));
	beads.OK_XY = find(x_sdb>y_sdb);		% y-movement should not be larger than x-movement
    beads.OK_XY = intersect(beads.OK_XY,beads.OK);
    % test for wrong direction
    if e.type == e.freq,
        fi = find(p.foscs>0.5);
        fi = fi(1);
    elseif  e.type == e.torq,  
        fi = find(abs(p.T) >= max(abs(p.T)));
    elseif  e.type == e.drug,  
        ind2 = find(p.t < tdrug - 5); %do not apply vaf during drug delivery
        ind3 = find(p.t > tdrug + 20);
        fi = [ind2 ind3];  %do not apply vaf during drug delivery        
    end;
    if length(fi)>1,   
        Gp_sign=median(median(p.Gp(fi,beads.OK)));
        Gpp_sign=median(median(p.Gpp(fi,beads.OK)));
        %beads.OK_XY=find(sign(median(p.Gp(fi,:))) == sign(Gp_sign)) ;   
        beads.OK_XY=find(mean(sign(p.Gp(fi,:))/sign(Gp_sign)) > 0.9); %Gp must be 90% of the time positive (or negative)
        beads.OK_XY = intersect(beads.OK_XY,beads.OK);
        beads.OK = beads.OK_XY;
        beads.OK_XY=find(mean(sign(p.eta(fi,:))) > 0.9);   %eta must be 90% of the time positive
        beads.OK_XY = intersect(beads.OK_XY,beads.OK);      
    end;                    
    beads.OK = beads.OK_XY;
else beads.OK_XY = beads.OK;
end;

disp('bad beads excluded');
beads.stats.OK = length(beads.OK);
% beads.percentIncluded = beads.stats.OK/beads.totalN*100;

beads.stats.total = beads.keptN;
beads.stats.moved = length(beads.moved);
beads.stats.OK_harmonic_frac = length(beads.OK_harmonic_frac);
beads.stats.OK_3rd_harmonic = length(beads.OK_3rd_harmonic);
beads.stats.OK_vaf = length(beads.OK_vaf);
beads.stats.OK = length(beads.OK);
beads.stats.OK_XY = length(beads.OK_XY);
beads.stats.OK_reproducible = length(beads.OK_reproducible);
beads.stats.percentIncluded = beads.stats.OK/beads.stats.total*100;
beads.statsArray = [beads.stats.total beads.stats.moved beads.stats.OK_harmonic_frac beads.stats.OK_3rd_harmonic beads.stats.OK_vaf beads.stats.OK_XY beads.stats.OK_reproducible beads.stats.OK beads.stats.percentIncluded];
a = beads.statsArray;
beads.stats_str = ['beads:' num2str(a(1)) ', moved:' num2str(a(2)) ', Hrm frac:' num2str(a(3)) ', 3rd Hrm:' num2str(a(4)) ', vaf:' num2str(a(5)) ', X>Y:' num2str(a(6)) ', reprod:' num2str(a(7)) ', OK:' num2str(a(8)) ', OK:' num2str(fix(a(9)*10)/10) '%'];

disp([num2str(beads.stats.OK) ' beads within limits of ' num2str(beads.stats.total) ' total number of beads ' num2str(beads.stats.percentIncluded) ' %']);

clear i1bead i2bead i3bead fi ti Q

if beads.stats.OK == 0, 
   figure(1);
   
   title(beads.stats_str);
   
   figure(2);  clf;
   if e.type == e.freq,
      xx = abs(p.foscs); xtext = 'fosc Hz';
   elseif e.type == e.torq,
      xx = abs(p.T); xtext = 'T dyn/cm^2';
   elseif e.type == e.drug,
		xx = abs(p.t); xtext = 't sec';
   end
   yy = beads.OK_harmonic_frac;
   if yy > 0,
      subplot(221); waterfall(xx,yy,abs(p.X(:,yy))'); ylabel('beads'); xlabel(xtext); zlabel('|X| nm');
   end;   
   yy = beads.OK_3rd_harmonic;
   if yy > 0,
      subplot(222); waterfall(xx,yy,abs(p.X(:,yy))'); ylabel('beads'); xlabel(xtext); zlabel('|X| nm');
   end;
   yy = beads.OK_vaf;
   if yy > 0,
      subplot(223); waterfall(xx,yy,abs(p.X(:,yy))'); ylabel('beads'); xlabel(xtext); zlabel('|X| nm');
   end;
   yy = length(beads.kept);
   subplot(224); 
   if e.type == e.freq,
 		semilogx(xx,p.vaf(:,yy)'); axis([0.01 1000 50 100]);
   else
      plot(xx,p.vaf(:,yy)');  v = axis; axis([v(1:2) 50 100]);
   end;
   xlabel(xtext); ylabel('vaf% maxX bead');    % waterfall(abs(p.X(:,beads.OK))'); ylabel('beads'); xlabel('condition'); zlabel('|X| nm');
   title(['NO ' beads.stats_str]);
   pause(1)
   return;
end;
drawnow;
% params by mean of params;
% not as good, since the noise is amplified

GpMedian1Hz_beadsOK = 1;
GppMedian1Hz_beadsOK = 1;
% Ben: nomalize results for f = 1 Hz or maximum torque
if e.type == e.drug, 
   if do.Ben == 0, 
      disp('switching to normalization');
      do.Benswitch = 1
    end;
end;
if do.Ben,
   if e.type == e.freq,
      fi = find(p.foscs>0.5);
      fi = fi(1);
   elseif  e.type == e.torq,  
      fi = find(abs(p.T) >= max(abs(p.T)));
   elseif  e.type == e.drug,  
      fi = find(p.t<tdrug);  %all cycles before drug addition
   end;
   if length(fi)>1,
   	Gp_norm =  abs(median(p.Gp(fi,:)));	% median over all cycles before drug addition
      Gpp_norm =  abs(median(p.Gpp(fi,:)));
   	GpMedian1Hz_beadsOK = median(mean(p.Gp(fi,beads.OK)));
   	GppMedian1Hz_beadsOK = median(mean(p.Gpp(fi,beads.OK)));      
   else
   	Gp_norm =  abs(p.Gp(fi,:));	% 
      Gpp_norm =  abs(p.Gpp(fi,:));
   	GpMedian1Hz_beadsOK = median(p.Gp(fi,beads.OK));
   	GppMedian1Hz_beadsOK = median(p.Gpp(fi,beads.OK));      
   end;

   for k=1:nb,
      p.Gp(k,:) = p.Gp(k,:) ./ Gp_norm;
      p.Gpp(k,:) = p.Gpp(k,:) ./ Gpp_norm;
   end;
   disp('Gp and Gpp normalized');
else
   Gp_norm = ones(size(p.Gp(1,:)));
   Gpp_norm =  Gp_norm;
end;
if do.Benswitch == 1,
   do.Benswitch = 0;
end;
   
k = beads.OK;							% do means of only the OK beads
pick.mean_of_params = 1;
if pick.mean_of_params, 
   if beads.stats.OK > 1, 
      p.X_mb = median(abs(p.X(:,k)'));  
      p.SRate_mb = median(p.SRate(:,k)'); 
      p.Gp_mb = median(p.Gp(:,k)')*GpMedian1Hz_beadsOK;
      p.Gpp_mb = median(p.Gpp(:,k)')*GppMedian1Hz_beadsOK;
      p.eta_mb = median(p.eta(:,k)');
      p.R_mb = median(p.R(:,k)');
      p.X_sdb = std(p.X(:,k)');
      p.Gp_sdb = std(p.Gp(:,k)')*GpMedian1Hz_beadsOK;
      p.Gpp_sdb = std(p.Gpp(:,k)')*GppMedian1Hz_beadsOK;
      p.eta_sdb = std(p.eta(:,k)');
      p.R_sdb = std(p.R(:,k)');
      p.SRate_sdb = std(p.SRate(:,k)');
   end;
end; 
if beads.stats.OK == 1,
   p.X_mb = p.X(:,k)';   
   p.SRate_mb = p.SRate(:,k)';
   p.Gp_mb = p.Gp(:,k)';
   p.Gpp_mb = p.Gpp(:,k)';
   p.eta_mb = p.eta(:,k)';
   p.R_mb = p.R(:,k)';
   p.X_sdb = zeros(size(p.X_mb));
   p.Gp_sdb = zeros(size(p.X_mb));
   p.Gpp_sdb = zeros(size(p.X_mb));
   p.eta_sdb = zeros(size(p.X_mb));
   p.R_sdb = zeros(size(p.X_mb));
   p.SRate_sdb = zeros(size(p.X_mb));
end;   
disp('median and stdev computed');

%this is to exclude "temporarily misbehaving" beads from being included in the *.ave file
pick.mean_by_X = 0;
if pick.mean_by_X
   for k = 1:nb,	% number of blocks (conditions)
      beads.OK_Gp = find(sign(p.Gp(k,:)) == sign(Gp_sign)) ;   
      beads.OK_Gpp = find(sign(p.Gpp(k,:)) == sign(Gpp_sign)) ; 
      beads.OK_temp = find(p.vaf(k,:) > Lim.vaf_low);
      beads.OK_temp = intersect(beads.OK,beads.OK_temp);
      beads.OK_temp = intersect(beads.OK_Gp,beads.OK_temp);
      beads.OK_temp = intersect(beads.OK_Gpp,beads.OK_temp);
      if length(beads.OK_temp) > 1
         p.X_mb(k) = median(p.X(k,beads.OK_temp));   
      	p.SRate_mb(k) = median(p.SRate(k,beads.OK_temp)); 
      	p.Gp_mb(k) = median(p.Gp(k,beads.OK_temp))*GpMedian1Hz_beadsOK;
      	p.Gpp_mb(k) = median(p.Gpp(k,beads.OK_temp))*GppMedian1Hz_beadsOK;
      	p.eta_mb(k) = median(p.eta(k,beads.OK_temp));
      	p.R_mb(k) = median(p.R(k,beads.OK_temp));
      	p.X_sdb(k) = std(p.X(k,beads.OK_temp));
      	p.Gp_sdb(k) = std(p.Gp(k,beads.OK_temp))*GpMedian1Hz_beadsOK;
      	p.Gpp_sdb(k) = std(p.Gpp(k,beads.OK_temp))*GppMedian1Hz_beadsOK;
      	p.eta_sdb(k) = std(p.eta(k,beads.OK_temp));
      	p.R_sdb(k) = std(p.R(k,beads.OK_temp));
      	p.SRate_sdb(k) = std(p.SRate(k,beads.OK_temp));
      elseif  length(beads.OK_temp) == 1 
         p.X_mb(k) = median(p.X(k,beads.OK_temp));   
      	p.SRate_mb(k) = median(p.SRate(k,beads.OK_temp)); 
      	p.Gp_mb(k) = median(p.Gp(k,beads.OK_temp))*GpMedian1Hz_beadsOK;
      	p.Gpp_mb(k) = median(p.Gpp(k,beads.OK_temp))*GppMedian1Hz_beadsOK;
      	p.eta_mb(k) = median(p.eta(k,beads.OK_temp));
      	p.R_mb(k) = median(p.R(k,beads.OK_temp));
      	p.X_sdb(k) = 0;
      	p.Gp_sdb(k) = 0;
      	p.Gpp_sdb(k) = 0;
      	p.eta_sdb(k) = 0;
      	p.R_sdb(k) = 0;
      	p.SRate_sdb(k) = 0;         
      else
         p.X_mb(k) = 0;   
      	p.SRate_mb(k) = 0; 
      	p.Gp_mb(k) = 0;
      	p.Gpp_mb(k) = 0;
      	p.eta_mb(k) = 0;
      	p.R_mb(k) = 0;
      	p.X_sdb(k) = 0;
      	p.Gp_sdb(k) = 0;
      	p.Gpp_sdb(k) = 0;
      	p.eta_sdb(k) = 0;
      	p.R_sdb(k) = 0;
      	p.SRate_sdb(k) = 0;           
      end;      
   end;
   disp('median and stdev computed based on individual blocks (exclude all blocks with low vaf)');
end;

if figures(1),
   % show which beads are analyzed on the figure
   % first get the loop from the largest Torque or the lowest frequency;
   disp('draw figure 1'); 
   figure(1); hold on;
   if e.type == e.freq,
      [temp k] = min(p.foscs(1:floor(nb/2)) < 1);
   else
      [temp k] = max(abs(p.T));
   end;
   a1 = b.fi(k);
   a2 = b.fi(k+1);   % for just one loop
   a1 = 1;
   a2 = b.fi(length(b.fi))-1;		% for all loops
   gain = 0.1;
   for k = 1:beads.keptN,   
      px = (ones(a2-a1+1,1)*xo(k) + gain*xf(a1:a2,k))/nmPerPixel;
      py = (ones(a2-a1+1,1)*yo(k) + gain*yf(a1:a2,k))/nmPerPixel;
      if ismember(k,beads.OK),
         %% define zero at lower corner
         plot(px, py,'r');
         h = text(xo(k)/nmPerPixel,yo(k)/nmPerPixel,[' ' num2str(k) '+']);
         set(h,'color','red');
      else
         plot(px,py,'b');
         h = text(xo(k)/nmPerPixel,yo(k)/nmPerPixel,[' ' num2str(k) '-']);
         set(h,'color','blue');
      end;
   end;
   a = beads.stats;
   title(beads.stats_str);
   clear a a1 a2 px py h;
   drawnow;
end;


% output the average results
disp('file output average results'); 
fileName.root = strtok(fileName.raw,'.');
fileName.ave = [fileName.root '.ave'];
fileName.AnalyzeDate = date;
temp = clock; 
fileName.AnalyzeTime = [num2str(temp(4)) ':' num2str(temp(5)) ':' num2str(fix(temp(6)))];

fid = fopen(fileName.ave,'w');
fprintf(fid,'%s\t',fileName.ave);  fprintf(fid,'%s\t','Created'); fprintf(fid,'%s\t',fileName.CreateDate); fprintf(fid,'%s\t',fileName.CreateTime); 
fprintf(fid,'%s\t','Analyzed'); fprintf(fid,'%s\t',fileName.AnalyzeDate); fprintf(fid,'%s\t',fileName.AnalyzeTime); fprintf(fid,'%s\n',' ');  
fields.beadsstats = fieldnames(beads.stats);
fields.Lim = fieldnames(Lim);
for k = 1:length(fields.beadsstats), fprintf(fid,'%s\t',char(fields.beadsstats(k))); end;  %names
fprintf(fid,'%s\t','C_bead dyn/cm^2/G'); fprintf(fid,'%s\n',' ');  
for k = 1:length(fields.beadsstats), fprintf(fid,'%g\t',getfield(beads.stats,char(fields.beadsstats(k)))); end;  %values
fprintf(fid,'%g\t',C_bead); fprintf(fid,'%s\n',' '); 
for k = 1:length(fields.Lim), fprintf(fid,'%s\t',char(fields.Lim(k))); end;
fprintf(fid,'%s\n',' '); 
for k = 1:length(fields.Lim), fprintf(fid,'%g\t',getfield(Lim,char(fields.Lim(k)))); end;
fprintf(fid,'%s\n',' '); 

fprintf(fid,'%s\n',' '); % n is for a new line
% the minus sign is for left justified, t is for a tab
fprintf(fid,'%-s\t','t sec','fosc Hz','I mA','T dyn/cm^2','amp nm','amp_se','Gp','Gp_se','Gpp','Gpp_se','eta','eta_se','R','R_se','Srate','Srate_se','time','mean square displacement','frequ','fft_x','fft_y'); 
fprintf(fid,'%s\n',' ');
for k = 1:length(p.foscs),
   Nsqr = beads.stats.OK^0.5;
   fprintf(fid,'%g\t',p.t(k));
   fprintf(fid,'%g\t',p.foscs(k));					
   fprintf(fid,'%g\t',abs(p.I(k)));				
   fprintf(fid,'%g\t',abs(p.T(k)));
   if beads.stats.OK > 0,
      fprintf(fid,'%g\t',abs(p.X_mb(k)));				
      fprintf(fid,'%g\t',p.X_sdb(k)/Nsqr);		
      fprintf(fid,'%g\t',p.Gp_mb(k));				
      fprintf(fid,'%g\t',p.Gp_sdb(k)/Nsqr);		
      fprintf(fid,'%g\t',p.Gpp_mb(k));		
      fprintf(fid,'%g\t',p.Gpp_sdb(k)/Nsqr);		
      fprintf(fid,'%g\t',p.eta_mb(k));		
      fprintf(fid,'%g\t',p.eta_sdb(k)/Nsqr);		
      fprintf(fid,'%g\t',p.R_mb(k));		
      fprintf(fid,'%g\t',p.R_sdb(k)/Nsqr);
      fprintf(fid,'%g\t',p.SRate_mb(k));
      fprintf(fid,'%g\t',p.SRate_sdb(k)/Nsqr);
   end;
   fprintf(fid,'%s\n',' ');
end;
fclose(fid);

% output the results beads by bead %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('file output bead by bead'); 
fileName.beads = [fileName.root '.beads'];
fid = fopen(fileName.beads,'w');

fprintf(fid,'%s\t',fileName.beads); fprintf(fid,'%s\t','Created'); fprintf(fid,'%s\t',fileName.CreateDate); fprintf(fid,'%s\t',fileName.CreateTime);
fprintf(fid,'%s\t','Analyzed'); fprintf(fid,'%s\t',fileName.AnalyzeDate); fprintf(fid,'%s\t',fileName.AnalyzeTime); fprintf(fid,'%s\n',' ');  
fields.beadsstats = fieldnames(beads.stats);
fields.Lim = fieldnames(Lim);
for k = 1:length(fields.beadsstats), fprintf(fid,'%s\t',char(fields.beadsstats(k))); end;  %names
fprintf(fid,'%s\t','C_bead dyn/cm^2/G'); fprintf(fid,'%s\n',' ');  
for k = 1:length(fields.beadsstats), fprintf(fid,'%g\t',getfield(beads.stats,char(fields.beadsstats(k)))); end;  %values
fprintf(fid,'%g\t',C_bead); fprintf(fid,'%s\n',' '); 
for k = 1:length(fields.Lim), fprintf(fid,'%s\t',char(fields.Lim(k))); end;
fprintf(fid,'%s\n',' '); 
for k = 1:length(fields.Lim), fprintf(fid,'%g\t',getfield(Lim,char(fields.Lim(k)))); end;
fprintf(fid,'%s\n',' '); 

% values
do_new_bead_file_format = 1;		% bead xo yo amp Gp Gpp vaf
if do_new_bead_file_format & (e.type == e.freq) | (e.type == e.torq),
   fprintf(fid,'%-s\t','','','','','');   
   for k = 1:5,
      fprintf(fid,'%s\t','fosc Hz'); fprintf(fid,'%g\t',p.foscs);
	end;
	fprintf(fid,'%s\n',' ');
	fprintf(fid,'%-s\t','','','','','');
	for k = 1:5,
   	fprintf(fid,'%s\t','T dyn/cm2'); fprintf(fid,'%g\t',abs(p.T));
	end;
	fprintf(fid,'%s\n',' ');
	fprintf(fid,'%-s\t','bead#','xo pix','yo pix','Gp_norm','Gpp_norm',' ');
	fprintf(fid,'%-s\t','amp nm');    for k = 1:length(p.foscs), fprintf(fid,'%s\t',' '); end;
	fprintf(fid,'%-s\t','Gp');        for k = 1:length(p.foscs), fprintf(fid,'%s\t',' '); end;
   fprintf(fid,'%-s\t','Gpp');       for k = 1:length(p.foscs), fprintf(fid,'%s\t',' '); end;
   fprintf(fid,'%-s\t','eta');       for k = 1:length(p.foscs), fprintf(fid,'%s\t',' '); end;
	fprintf(fid,'%-s\t','vaf%');      for k = 1:length(p.foscs), fprintf(fid,'%s\t',' '); end;
	fprintf(fid,'%s\n',' ');
	
	which_beads = beads.kept; 			% for all beads not lost
	which_beads = beads.stats.OK; 	% for all OK beads
   for k = 1:which_beads,
      beadix = beads.OK(k);
      fprintf(fid,'%g\t',beadix);
      fprintf(fid,'%g\t',xo(beadix)/nmPerPixel);	% pixels
      fprintf(fid,'%g\t',yo(beadix)/nmPerPixel);	% pixels
      fprintf(fid,'%g\t',Gp_norm(beadix));	% 
      fprintf(fid,'%g\t',Gpp_norm(beadix));	% 
      fprintf(fid,'%s\t',' ');
      fprintf(fid,'%g\t',abs(p.X(:,beadix))); fprintf(fid,'%s\t',' ');
      fprintf(fid,'%g\t',p.Gp(:,beadix)); fprintf(fid,'%s\t',' ');
      fprintf(fid,'%g\t',p.Gpp(:,beadix)); fprintf(fid,'%s\t',' ');
      fprintf(fid,'%g\t',p.eta(:,beadix)); fprintf(fid,'%s\t',' ');
      fprintf(fid,'%g\t',p.vaf(:,beadix)); fprintf(fid,'%s\t',' ');
      fprintf(fid,'%s\n',' ');
   end;
   fclose(fid);
else

   fprintf(fid,'%-s\t','bead#'); fprintf(fid,'%s\t',' '); fprintf(fid,'%s\t',' '); fprintf(fid,'%s\t',' ');
   for k = 1:beads.stats.OK,
      fprintf(fid,'%g\t',beads.OK(k)); fprintf(fid,'%s\t',' '); fprintf(fid,'%s\t',' '); fprintf(fid,'%s\t',' ');
   end;
   fprintf(fid,'%s\n',' ');
   
   fprintf(fid,'%-s\t','t sec','fosc Hz','I mA','T dyn/cm^2');
   for k = 1:beads.stats.OK,
      fprintf(fid,'%s\t','amp nm','Gp','Gpp','vaf%'); 
   end;
   fprintf(fid,'%s\n',' ');
   for k = 1:length(p.t),
      fprintf(fid,'%g\t',p.t(k));
      fprintf(fid,'%g\t',p.foscs(k));
      fprintf(fid,'%g\t',abs(p.I(k)));				
      fprintf(fid,'%g\t',abs(p.T(k)));
      for kk = 1:length(beads.OK),
         a1 = beads.OK(kk);
         fprintf(fid,'%g\t',abs(p.X(k,a1)));
         fprintf(fid,'%g\t',(p.Gp(k,a1)));
         fprintf(fid,'%g\t',(p.Gpp(k,a1)));
         fprintf(fid,'%g\t',(p.vaf(k,a1)));
      end;
      fprintf(fid,'%s\n',' ');
   end;   
   fclose(fid);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% output the results beads by bead only Gp for Excel  %%%%%%%%%%%%%%%%%%%%%%% 
disp('file output bead by bead Gp'); 
fileName.beads = [fileName.root '.bdsgp'];
fid = fopen(fileName.beads,'w');

fprintf(fid,'%s\t',fileName.beads); fprintf(fid,'%s\t','Created'); fprintf(fid,'%s\t',fileName.CreateDate); fprintf(fid,'%s\t',fileName.CreateTime);
fprintf(fid,'%s\t','Analyzed'); fprintf(fid,'%s\t',fileName.AnalyzeDate); fprintf(fid,'%s\t',fileName.AnalyzeTime); fprintf(fid,'%s\n',' ');  
fields.beadsstats = fieldnames(beads.stats);
fields.Lim = fieldnames(Lim);
for k = 1:length(fields.beadsstats), fprintf(fid,'%s\t',char(fields.beadsstats(k))); end;  %names
fprintf(fid,'%s\t','C_bead dyn/cm^2/G'); fprintf(fid,'%s\n',' ');  
for k = 1:length(fields.beadsstats), fprintf(fid,'%g\t',getfield(beads.stats,char(fields.beadsstats(k)))); end;  %values
fprintf(fid,'%g\t',C_bead); fprintf(fid,'%s\n',' '); 
for k = 1:length(fields.Lim), fprintf(fid,'%s\t',char(fields.Lim(k))); end;
fprintf(fid,'%s\n',' '); 
for k = 1:length(fields.Lim), fprintf(fid,'%g\t',getfield(Lim,char(fields.Lim(k)))); end;
fprintf(fid,'%s\n',' '); 

% values
do_new_bead_file_format = 1;		% bead xo yo amp Gp Gpp vaf
if do_new_bead_file_format & (e.type == e.freq) | (e.type == e.torq),
   fprintf(fid,'%-s\t','','','','','');   
   for k = 1:5,
      fprintf(fid,'%s\t','fosc Hz'); fprintf(fid,'%g\t',p.foscs);
	end;
	fprintf(fid,'%s\n',' ');
	fprintf(fid,'%-s\t','','','','','');
	for k = 1:5,
   	fprintf(fid,'%s\t','T Pa'); fprintf(fid,'%g\t',abs(p.T));
	end;
	fprintf(fid,'%s\n',' ');
	fprintf(fid,'%-s\t','bead#','xo pix','yo pix','Gp_norm','Gpp_norm',' ');
	fprintf(fid,'%-s\t','amp nm');    for k = 1:length(p.foscs), fprintf(fid,'%s\t',' '); end;
	fprintf(fid,'%-s\t','Gp');        for k = 1:length(p.foscs), fprintf(fid,'%s\t',' '); end;
   fprintf(fid,'%-s\t','Gpp');       for k = 1:length(p.foscs), fprintf(fid,'%s\t',' '); end;
   fprintf(fid,'%-s\t','eta');       for k = 1:length(p.foscs), fprintf(fid,'%s\t',' '); end;
	fprintf(fid,'%-s\t','vaf%');      for k = 1:length(p.foscs), fprintf(fid,'%s\t',' '); end;
	fprintf(fid,'%s\n',' ');
	
	which_beads = beads.kept; 			% for all beads not lost
	which_beads = beads.stats.OK; 	% for all OK beads
   for k = 1:which_beads,
      beadix = beads.OK(k);
      fprintf(fid,'%g\t',beadix);
      fprintf(fid,'%g\t',xo(beadix)/nmPerPixel);	% pixels
      fprintf(fid,'%g\t',yo(beadix)/nmPerPixel);	% pixels
      fprintf(fid,'%g\t',Gp_norm(beadix));	% 
      fprintf(fid,'%g\t',Gpp_norm(beadix));	% 
      fprintf(fid,'%s\t',' ');
      fprintf(fid,'%g\t',abs(p.X(:,beadix))); fprintf(fid,'%s\t',' ');
      fprintf(fid,'%g\t',p.Gp(:,beadix)); fprintf(fid,'%s\t',' ');
      fprintf(fid,'%g\t',p.Gpp(:,beadix)); fprintf(fid,'%s\t',' ');
      fprintf(fid,'%g\t',p.eta(:,beadix)); fprintf(fid,'%s\t',' ');
      fprintf(fid,'%g\t',p.vaf(:,beadix)); fprintf(fid,'%s\t',' ');
      fprintf(fid,'%s\n',' ');
   end;
   fclose(fid);
else

   fprintf(fid,'%-s\t','bead#'); fprintf(fid,'%s\t',' '); fprintf(fid,'%s\t',' '); fprintf(fid,'%s\t',' ');
   for k = 1:beads.stats.OK,
      fprintf(fid,'%g\t',beads.OK(k)); 
   end;
   fprintf(fid,'%s\n',' ');
   
   fprintf(fid,'%-s\t','t sec','fosc Hz','I mA','T Pa');
   for k = 1:beads.stats.OK,
      fprintf(fid,'%s\t','Gp'); 
   end;
   fprintf(fid,'%s\n',' ');
   for k = 1:length(p.t),
      fprintf(fid,'%g\t',p.t(k));
      fprintf(fid,'%g\t',p.foscs(k));
      fprintf(fid,'%g\t',abs(p.I(k)));				
      fprintf(fid,'%g\t',abs(p.T(k)));
      for kk = 1:length(beads.OK),
         a1 = beads.OK(kk);
         %fprintf(fid,'%g\t',abs(p.X(k,a1)));
         fprintf(fid,'%g\t',(p.Gp(k,a1)));
         %fprintf(fid,'%g\t',(p.Gpp(k,a1)));
         %fprintf(fid,'%g\t',(p.vaf(k,a1)));
      end;
      fprintf(fid,'%s\n',' ');
   end;   
   fclose(fid);
end;




% output the results beads by bead only Gpp for Excel  %%%%%%%%%%%%%%%%%%%%%%% 
disp('file output bead by bead Gpp'); 
fileName.beads = [fileName.root '.bdsgpp'];
fid = fopen(fileName.beads,'w');

fprintf(fid,'%s\t',fileName.beads); fprintf(fid,'%s\t','Created'); fprintf(fid,'%s\t',fileName.CreateDate); fprintf(fid,'%s\t',fileName.CreateTime);
fprintf(fid,'%s\t','Analyzed'); fprintf(fid,'%s\t',fileName.AnalyzeDate); fprintf(fid,'%s\t',fileName.AnalyzeTime); fprintf(fid,'%s\n',' ');  
fields.beadsstats = fieldnames(beads.stats);
fields.Lim = fieldnames(Lim);
for k = 1:length(fields.beadsstats), fprintf(fid,'%s\t',char(fields.beadsstats(k))); end;  %names
fprintf(fid,'%s\t','C_bead dyn/cm^2/G'); fprintf(fid,'%s\n',' ');  
for k = 1:length(fields.beadsstats), fprintf(fid,'%g\t',getfield(beads.stats,char(fields.beadsstats(k)))); end;  %values
fprintf(fid,'%g\t',C_bead); fprintf(fid,'%s\n',' '); 
for k = 1:length(fields.Lim), fprintf(fid,'%s\t',char(fields.Lim(k))); end;
fprintf(fid,'%s\n',' '); 
for k = 1:length(fields.Lim), fprintf(fid,'%g\t',getfield(Lim,char(fields.Lim(k)))); end;
fprintf(fid,'%s\n',' '); 

% values
do_new_bead_file_format = 1;		% bead xo yo amp Gp Gpp vaf
if do_new_bead_file_format & (e.type == e.freq) | (e.type == e.torq),
   fprintf(fid,'%-s\t','','','','','');   
   for k = 1:5,
      fprintf(fid,'%s\t','fosc Hz'); fprintf(fid,'%g\t',p.foscs);
	end;
	fprintf(fid,'%s\n',' ');
	fprintf(fid,'%-s\t','','','','','');
	for k = 1:5,
   	fprintf(fid,'%s\t','T Pa'); fprintf(fid,'%g\t',abs(p.T));
	end;
	fprintf(fid,'%s\n',' ');
	fprintf(fid,'%-s\t','bead#','xo pix','yo pix','Gp_norm','Gpp_norm',' ');
	fprintf(fid,'%-s\t','amp nm');    for k = 1:length(p.foscs), fprintf(fid,'%s\t',' '); end;
	fprintf(fid,'%-s\t','Gp');        for k = 1:length(p.foscs), fprintf(fid,'%s\t',' '); end;
   fprintf(fid,'%-s\t','Gpp');       for k = 1:length(p.foscs), fprintf(fid,'%s\t',' '); end;
   fprintf(fid,'%-s\t','eta');       for k = 1:length(p.foscs), fprintf(fid,'%s\t',' '); end;
	fprintf(fid,'%-s\t','vaf%');      for k = 1:length(p.foscs), fprintf(fid,'%s\t',' '); end;
	fprintf(fid,'%s\n',' ');
	
	which_beads = beads.kept; 			% for all beads not lost
	which_beads = beads.stats.OK; 	% for all OK beads
   for k = 1:which_beads,
      beadix = beads.OK(k);
      fprintf(fid,'%g\t',beadix);
      fprintf(fid,'%g\t',xo(beadix)/nmPerPixel);	% pixels
      fprintf(fid,'%g\t',yo(beadix)/nmPerPixel);	% pixels
      fprintf(fid,'%g\t',Gp_norm(beadix));	% 
      fprintf(fid,'%g\t',Gpp_norm(beadix));	% 
      fprintf(fid,'%s\t',' ');
      fprintf(fid,'%g\t',abs(p.X(:,beadix))); fprintf(fid,'%s\t',' ');
      fprintf(fid,'%g\t',p.Gp(:,beadix)); fprintf(fid,'%s\t',' ');
      fprintf(fid,'%g\t',p.Gpp(:,beadix)); fprintf(fid,'%s\t',' ');
      fprintf(fid,'%g\t',p.eta(:,beadix)); fprintf(fid,'%s\t',' ');
      fprintf(fid,'%g\t',p.vaf(:,beadix)); fprintf(fid,'%s\t',' ');
      fprintf(fid,'%s\n',' ');
   end;
   fclose(fid);
else

   fprintf(fid,'%-s\t','bead#'); fprintf(fid,'%s\t',' '); fprintf(fid,'%s\t',' '); fprintf(fid,'%s\t',' ');
   for k = 1:beads.stats.OK,
      fprintf(fid,'%g\t',beads.OK(k)); 
   end;
   fprintf(fid,'%s\n',' ');
   
   fprintf(fid,'%-s\t','t sec','fosc Hz','I mA','T Pa');
   for k = 1:beads.stats.OK,
      fprintf(fid,'%s\t','Gpp'); 
   end;
   fprintf(fid,'%s\n',' ');
   for k = 1:length(p.t),
      fprintf(fid,'%g\t',p.t(k));
      fprintf(fid,'%g\t',p.foscs(k));
      fprintf(fid,'%g\t',abs(p.I(k)));				
      fprintf(fid,'%g\t',abs(p.T(k)));
      for kk = 1:length(beads.OK),
         a1 = beads.OK(kk);
         %fprintf(fid,'%g\t',abs(p.X(k,a1)));
         fprintf(fid,'%g\t',(p.Gpp(k,a1)));
         %fprintf(fid,'%g\t',(p.Gpp(k,a1)));
         %fprintf(fid,'%g\t',(p.vaf(k,a1)));
      end;
      fprintf(fid,'%s\n',' ');
   end;   
   fclose(fid);
end;




% output the results beads by bead only eta for Excel  %%%%%%%%%%%%%%%%%%%%%%% 
disp('file output bead by bead eta'); 
fileName.beads = [fileName.root '.bdseta'];
fid = fopen(fileName.beads,'w');

fprintf(fid,'%s\t',fileName.beads); fprintf(fid,'%s\t','Created'); fprintf(fid,'%s\t',fileName.CreateDate); fprintf(fid,'%s\t',fileName.CreateTime);
fprintf(fid,'%s\t','Analyzed'); fprintf(fid,'%s\t',fileName.AnalyzeDate); fprintf(fid,'%s\t',fileName.AnalyzeTime); fprintf(fid,'%s\n',' ');  
fields.beadsstats = fieldnames(beads.stats);
fields.Lim = fieldnames(Lim);
for k = 1:length(fields.beadsstats), fprintf(fid,'%s\t',char(fields.beadsstats(k))); end;  %names
fprintf(fid,'%s\t','C_bead dyn/cm^2/G'); fprintf(fid,'%s\n',' ');  
for k = 1:length(fields.beadsstats), fprintf(fid,'%g\t',getfield(beads.stats,char(fields.beadsstats(k)))); end;  %values
fprintf(fid,'%g\t',C_bead); fprintf(fid,'%s\n',' '); 
for k = 1:length(fields.Lim), fprintf(fid,'%s\t',char(fields.Lim(k))); end;
fprintf(fid,'%s\n',' '); 
for k = 1:length(fields.Lim), fprintf(fid,'%g\t',getfield(Lim,char(fields.Lim(k)))); end;
fprintf(fid,'%s\n',' '); 

% values
do_new_bead_file_format = 1;		% bead xo yo amp Gp Gpp vaf
if do_new_bead_file_format & (e.type == e.freq) | (e.type == e.torq),
   fprintf(fid,'%-s\t','','','','','');   
   for k = 1:5,
      fprintf(fid,'%s\t','fosc Hz'); fprintf(fid,'%g\t',p.foscs);
	end;
	fprintf(fid,'%s\n',' ');
	fprintf(fid,'%-s\t','','','','','');
	for k = 1:5,
   	fprintf(fid,'%s\t','T Pa'); fprintf(fid,'%g\t',abs(p.T));
	end;
	fprintf(fid,'%s\n',' ');
	fprintf(fid,'%-s\t','bead#','xo pix','yo pix','Gp_norm','Gpp_norm',' ');
	fprintf(fid,'%-s\t','amp nm');    for k = 1:length(p.foscs), fprintf(fid,'%s\t',' '); end;
	fprintf(fid,'%-s\t','Gp');        for k = 1:length(p.foscs), fprintf(fid,'%s\t',' '); end;
   fprintf(fid,'%-s\t','Gpp');       for k = 1:length(p.foscs), fprintf(fid,'%s\t',' '); end;
   fprintf(fid,'%-s\t','eta');       for k = 1:length(p.foscs), fprintf(fid,'%s\t',' '); end;
	fprintf(fid,'%-s\t','vaf%');      for k = 1:length(p.foscs), fprintf(fid,'%s\t',' '); end;
	fprintf(fid,'%s\n',' ');
	
	which_beads = beads.kept; 			% for all beads not lost
	which_beads = beads.stats.OK; 	% for all OK beads
   for k = 1:which_beads,
      beadix = beads.OK(k);
      fprintf(fid,'%g\t',beadix);
      fprintf(fid,'%g\t',xo(beadix)/nmPerPixel);	% pixels
      fprintf(fid,'%g\t',yo(beadix)/nmPerPixel);	% pixels
      fprintf(fid,'%g\t',Gp_norm(beadix));	% 
      fprintf(fid,'%g\t',Gpp_norm(beadix));	% 
      fprintf(fid,'%s\t',' ');
      fprintf(fid,'%g\t',abs(p.X(:,beadix))); fprintf(fid,'%s\t',' ');
      fprintf(fid,'%g\t',p.Gp(:,beadix)); fprintf(fid,'%s\t',' ');
      fprintf(fid,'%g\t',p.Gpp(:,beadix)); fprintf(fid,'%s\t',' ');
      fprintf(fid,'%g\t',p.eta(:,beadix)); fprintf(fid,'%s\t',' ');
      fprintf(fid,'%g\t',p.vaf(:,beadix)); fprintf(fid,'%s\t',' ');
      fprintf(fid,'%s\n',' ');
   end;
   fclose(fid);
else

   fprintf(fid,'%-s\t','bead#'); fprintf(fid,'%s\t',' '); fprintf(fid,'%s\t',' '); fprintf(fid,'%s\t',' ');
   for k = 1:beads.stats.OK,
      fprintf(fid,'%g\t',beads.OK(k)); 
   end;
   fprintf(fid,'%s\n',' ');
   
   fprintf(fid,'%-s\t','t sec','fosc Hz','I mA','T Pa');
   for k = 1:beads.stats.OK,
      fprintf(fid,'%s\t','eta'); 
   end;
   fprintf(fid,'%s\n',' ');
   for k = 1:length(p.t),
      fprintf(fid,'%g\t',p.t(k));
      fprintf(fid,'%g\t',p.foscs(k));
      fprintf(fid,'%g\t',abs(p.I(k)));				
      fprintf(fid,'%g\t',abs(p.T(k)));
      for kk = 1:length(beads.OK),
         a1 = beads.OK(kk);
         %fprintf(fid,'%g\t',abs(p.X(k,a1)));
         fprintf(fid,'%g\t',(p.eta(k,a1)));
         %fprintf(fid,'%g\t',(p.Gpp(k,a1)));
         %fprintf(fid,'%g\t',(p.vaf(k,a1)));
      end;
      fprintf(fid,'%s\n',' ');
   end;   
   fclose(fid);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc;

if figures(2),
   disp('draw figure 2'); 
   figure(2); set(2,'Position',[1,30,450,300]); clf;
   if e.type == e.freq, 
      xtext = 'condition';    xtext2 = 'f Hz'; 
      subplot(331); semilogy(p.foscs,'-'); ylabel('f Hz'); xlabel(xtext);
      xx = p.foscs;
      subplot(332); semilogx(xx,abs(p.T),'-'); ylabel('|T| dyn/cm^2'); v=axis; axis([0.01 1000 0 v(4)]); set(gca,'XTick',[0.1  10 1000]);
      subplot(333); semilogx(xx,[abs(p.X_mb)' p.X_sdb'/Nsqr]); ylabel('|X| & se nm'); v=axis; axis([0.01 1000 0 v(4)]); set(gca,'XTick',[0.1  10 1000]);
      subplot(334); loglog(xx,[p.Gp_mb' p.Gp_sdb'/Nsqr]); ylabel('G'' & se '); v=axis; axis([0.01 1000 0 v(4)]); set(gca,'XTick',[0.1  10 1000]);
      subplot(335); loglog(xx,[p.Gpp_mb' p.Gpp_sdb'/Nsqr]); ylabel('G" & se'); v=axis; axis([0.01 1000 0 v(4)]); set(gca,'XTick',[0.1  10 1000]);
      subplot(336); semilogx(xx,[p.eta_mb' p.eta_sdb'/Nsqr]); ylabel('eta & se '); v=axis; axis([0.01 1000 0 v(4)]); set(gca,'XTick',[0.1  10 1000]);
      subplot(337); loglog(xx,[p.R_mb' p.R_sdb'/Nsqr]); ylabel('R & se'); xlabel(xtext2); v=axis; axis([0.01 1000 0 v(4)]); set(gca,'XTick',[0.1  10 1000]);
      subplot(338); loglog(xx,[p.SRate_mb' p.SRate_sdb'/Nsqr]); ylabel('shear rate nm/s'); xlabel(xtext2); v=axis; axis([0.01 1000 0 v(4)]); set(gca,'XTick',[0.1  10 1000]);
      title(titleStr);
   elseif e.type == e.torq,
      xtext = 'condition';    xtext2 = 'T dyn/cm^2'; 
      subplot(331); plot(p.foscs,'-'); ylabel('f Hz'); xlabel(xtext);
      xx = abs(p.T);
      subplot(332); plot(abs(p.T),'-'); ylabel('|T| dyn/cm^2'); xlabel(xtext); 
      subplot(333); plot([abs(p.X_mb)' p.X_sdb'/Nsqr]); ylabel('|X| & se nm'); xlabel(xtext); 
      subplot(334); plot([p.Gp_mb' p.Gp_sdb'/Nsqr]); ylabel('G'' & se '); xlabel(xtext); v=axis; axis([v(1:2) 0 v(4)]);
      subplot(335); plot([p.Gpp_mb' p.Gpp_sdb'/Nsqr]); ylabel('G" & se'); xlabel(xtext); v=axis; axis([v(1:2) 0 v(4)]);
      subplot(336); plot([p.eta_mb' p.eta_sdb'/Nsqr]); ylabel('eta & se '); xlabel(xtext); v=axis; axis([v(1:2) 0 v(4)]);
      subplot(337); plot([p.R_mb' p.R_sdb'/Nsqr]); ylabel('R & se'); xlabel(xtext); v=axis; axis([v(1:2) 0 v(4)]);
      subplot(338); plot([p.SRate_mb' p.SRate_sdb'/Nsqr]); ylabel('shear rate nm/s'); xlabel(xtext);
      title(titleStr);
   elseif e.type == e.drug,
      xtext = 'condition'; xtext2 = 'time s'; 
      subplot(331); plot(p.t,p.foscs,'-'); ylabel('f Hz');
      subplot(332); plot(p.t,abs(p.T)); ylabel('|T| dyn/cm^2'); 
      subplot(333); plot(p.t,[abs(p.X_mb)' p.X_sdb'/Nsqr]); ylabel('|X| & se nm'); v=axis; axis([v(1:2) 0 v(4)]);
      subplot(334); plot(p.t,[p.Gp_mb' p.Gp_sdb'/Nsqr]); ylabel('G'' & se '); v=axis; axis([v(1:2) 0 v(4)]);
      subplot(335); plot(p.t,[p.Gpp_mb' p.Gpp_sdb'/Nsqr]); ylabel('G" & se'); v=axis; axis([v(1:2) 0 v(4)]);
      subplot(336); plot(p.t,[p.eta_mb' p.eta_sdb'/Nsqr]); ylabel('eta & se '); xlabel(xtext); v=axis; axis([v(1:2) 0 2*median(p.eta_mb)]);
      subplot(337); plot(p.t,[p.R_mb' p.R_sdb'/Nsqr]); ylabel('R & se'); xlabel(xtext); v=axis; axis([v(1:2) 0 v(4)]);
      subplot(338); plot(p.t,[p.SRate_mb' p.SRate_sdb'/Nsqr]); ylabel('shear rate nm/s'); xlabel(xtext);
      title(titleStr);
   end;
end;
drawnow;
if figures(3)==0,
   disp('draw figure 3'); 
   figure(3); clf;
   if beads.keptN > 1,
      mesh(abs(p.X)'); ylabel('beads'); xlabel('condition'); zlabel('|X| nm');
   else
      plot(abs(p.X)'); xlabel('condition'); ylabel('|X| nm');
   end;
   title(titleStr);
end;
if figures(4)==0,
   disp('draw figure 4'); 
   figure(4); clf;
   if beads.OK > 0,
      yy = beads.OK;
      zz = abs(p.X(:,yy))';
      if e.type == e.torq, 
         xx = abs(p.T); xtext = 'T dyn/cm^2';
         waterfall(xx,yy,zz); ylabel('beads'); xlabel(xtext); zlabel('|X| nm');
      elseif e.type == e.freq,
         xx = abs(p.foscs); xtext = 'f Hz';
         waterfall(xx,yy,zz); ylabel('beads'); xlabel(xtext); zlabel('|X| nm'); set(gca,'xScale','log'); v = axis; axis([0.01 1000 v(3:6)]); set(gca,'XTick',[0.1  10 1000]);
      elseif e.type == e.drug,
         xx = abs(p.t); xtext = 't sec';
         waterfall(xx,yy,zz); ylabel('beads'); xlabel(xtext); zlabel('|X| nm');
      else
         disp(['invalid experiment type, cannot plot figure 4: exiting']); return;
      end;
      % waterfall(abs(p.X(:,beads.OK))'); ylabel('beads'); xlabel('condition'); zlabel('|X| nm');
   end;
   title(titleStr);
end;
if figures(5)==0,
   disp('draw figure 5'); 
   figure(5); clf;
   if beads.OK > 0,
      yy = beads.OK;
      if e.type == e.torq, 
         xtext = 'T dyn/cm^2'; xx = abs(p.T);
         subplot(221); waterfall(xx,yy,abs(p.Gp(:,yy))'); ylabel('beads'); xlabel(xtext); zlabel('G'''); view([105 10]);
         subplot(222); waterfall(xx,yy,abs(p.Gpp(:,yy))'); ylabel('beads'); xlabel(xtext); zlabel('G"'); view([105 10]);
         subplot(223); waterfall(xx,yy,abs(p.R(:,yy))'); ylabel('beads'); xlabel(xtext); zlabel('R'); view([105 10]);
         subplot(224); waterfall(xx,yy,abs(p.eta(:,yy))'); ylabel('beads'); xlabel(xtext); zlabel('eta'); view([105 10]);
      elseif e.type == e.freq,
         xtext = 'f Hz'; xx = p.foscs;
         subplot(221); waterfall(xx,yy,abs(p.Gp(:,yy))'); ylabel('beads'); xlabel(xtext); zlabel('G'''); set(gca,'xScale','log'); v = axis; axis([0.01 1000 v(3:6)]); view([-25 25]); set(gca,'Ydir','reverse'); set(gca,'zScale','log');  set(gca,'XTick',[0.1  10 1000]);
         subplot(222); waterfall(xx,yy,abs(p.Gpp(:,yy))'); ylabel('beads'); xlabel(xtext); zlabel('G"'); set(gca,'xScale','log'); v = axis; axis([0.01 1000 v(3:6)]); view([-25 25]); set(gca,'Ydir','reverse'); set(gca,'zScale','log'); set(gca,'XTick',[0.1  10 1000]);
         subplot(223); waterfall(xx,yy,abs(p.R(:,yy))'); ylabel('beads'); xlabel(xtext); zlabel('R'); set(gca,'xScale','log'); v = axis; axis([0.01 1000 v(3:6)]); view([-25 25]); set(gca,'Ydir','reverse'); set(gca,'zScale','log'); set(gca,'XTick',[0.1  10 1000]);
         subplot(224); waterfall(xx,yy,abs(p.eta(:,yy))'); ylabel('beads'); xlabel(xtext); zlabel('eta'); set(gca,'xScale','log'); v = axis; axis([0.01 1000 v(3:6)]); view([-25 25]); set(gca,'Ydir','reverse'); set(gca,'zScale','log'); set(gca,'XTick',[0.1  10 1000]);
      elseif e.type == e.drug,
         xtext = 't sec';   xx = abs(p.t);
         subplot(221); waterfall(xx,yy,abs(p.Gp(:,yy))'); ylabel('beads'); xlabel(xtext); zlabel('G'''); view([-25 25]); set(gca,'Ydir','reverse'); set(gca,'XTick',[0.1  10 1000]);
         subplot(222); waterfall(xx,yy,abs(p.Gpp(:,yy))'); ylabel('beads'); xlabel(xtext); zlabel('G"'); view([-25 25]); set(gca,'Ydir','reverse'); set(gca,'XTick',[0.1  10 1000]);
         subplot(223); waterfall(xx,yy,abs(p.R(:,yy))'); ylabel('beads'); xlabel(xtext); zlabel('R'); view([-25 25]); set(gca,'Ydir','reverse'); set(gca,'XTick',[0.1  10 1000]);
         subplot(224); waterfall(xx,yy,abs(p.eta(:,yy))'); ylabel('beads'); xlabel(xtext); zlabel('eta'); view([-25 25]); set(gca,'Ydir','reverse'); set(gca,'XTick',[0.1  10 1000]);
      else
         disp(['invalid experiment type, cannot plot figure 5: exiting']); return;
      end;
      title(titleStr);
   end;
end;
if figures(6)==0,
   disp('draw figure 6'); 
   figure(6); 
   subplot(211); plot(p.vaf); v = axis; axis([v(1:2) 60 100]);
   ylabel('vaf% all beads'); title(titleStr);
   subplot(212); plot(p.vaf(:,beads.OK)); v = axis; axis([v(1:2) 80 100]);
   ylabel('vaf% kept beads'); xlabel('condition');
end;

if figures(7),
   disp('draw figure 7'); 
   figure(7);
   set(7,'Position',[501,30,450,300])
   % plot the best bead with model
   
   [temp is] = sort(mean(p.vaf(:,beads.OK)));
   l=(length(is));
   while l>0,
   besti = beads.OK(is(l));
   al=temp(l);
   if e.type == e.freq, 
        ti = [1:floor(L/2)];
% average loops of same frequency
		
		for k = 1:max_fi,				 % skips first cycle and last logf value
         [common_f common_fi_up common_fi_down] = intersect(up_f,down_f);
         a1 = b.fi(k);
         a2 = b.fi(k+1) - 1;   
         if (length(common_fi_up) == length(common_fi_down)),
         	a3 = b.fi(nb-k+1);
            a4 = b.fi(nb-k+2) - 1; 
         else 
            a3 = a1;
            a4 = a2;
         end;
      	b.i = [i(a1:a2) i(a3:a4)];
      	b.xf = [xf(a1:a2,besti) xf(a3:a4,besti)];
      	number_of_cycles = floor(length(b.i)/16);
      	b_resh.i = reshape(b.i(1:16*number_of_cycles),16,number_of_cycles);
      	b_resh.xf = reshape(b.xf(1:16*number_of_cycles),16,number_of_cycles);
         if k == 1,
            b_mean.i = mean(b_resh.i');
            b_mean.xf = mean(b_resh.xf');
			else   
            b_mean.i = [b_mean.i mean(b_resh.i')];
            b_mean.xf = [b_mean.xf mean(b_resh.xf')];
         end;

   	end;
  		
   elseif e.type == e.torq,
      t1i = floor(10*L/50);
      t2i = t1i + 16*8;
      ti = [t1i:t2i];
   elseif e.type == e.drug,
      t1i = tdrugi - 5*points_per_cycle;
      t2i = max(find(t < 50+tdrug));
      ti = [t1i:t2i];
   else
      disp(['invalid experiement type, cannot plot time course']);
   end;
   if e.type == e.freq, 
      xx = p.foscs;
      subplot(221); semilogx(xx,p.eta(:,besti)); ylabel('eta'); v=axis; axis([0.01 1000 0 v(4)]); set(gca,'XTick',[0.1  10 1000]);
   else
      subplot(221); plot(t(ti),i(ti)); ylabel('i mA');
   end;
   subplot(223); plot(t(ti),xf(ti,besti),t(ti),x_mod(ti,besti)); 
   xlabel('time s'); ylabel('x & x_mod nm');
   if e.type == e.freq, 
      subplot(222); plot(b_mean.i,b_mean.xf); xlabel('i mA'); ylabel('xf nm'); hold on;
      subplot(222); plot(b_mean.i(1:16),b_mean.xf(1:16),'r-'); xlabel('i mA'); ylabel('xf nm'); drawnow;
      fileName.loop = [fileName.root '.loop'];
      fid = fopen(fileName.loop,'w');
      for k=1:length(b_mean.i),
         fprintf(fid,'%g\t',b_mean.i(k));
         fprintf(fid,'%g\t',b_mean.xf(k));
         fprintf(fid,'%s\n',' ');
         if mod(k,16) == 0,	%print first point of loop again, then insert empty line
         	fprintf(fid,'%g\t',b_mean.i(k-15));
         	fprintf(fid,'%g\t',b_mean.xf(k-15));
         	fprintf(fid,'%s\n',' ');            
            fprintf(fid,'%s\n',' ');
         end;
      end;
      fclose(fid);
   else 
       subplot(222); plot(i(ti),xf(ti,besti)); xlabel('i mA'); ylabel('xf nm');
   end;
   if e.type == e.freq,       
      subplot(224); loglog(p.foscs,p.Gp(:,besti)); ylabel('Gp'); v=axis; axis([0.01 1000 0 v(4)]); set(gca,'XTick',[0.1  10 1000]);
   else
      subplot(224); plot(i(ti),x_mod(ti,besti)); xlabel('i mA'); ylabel('x_mod nm');
   end;
   title([char(titleStr) ' mean VAF= ' num2str(al) '%, bead ' num2str(besti)]);
   kkk = 2;
   %kkk = input('press return to continue, 2 to terminate, or 1 to go a bead back ');
   l = l - 1;
   if kkk == 2,
      l = 0;
   elseif kkk == 1
      l = l + 2;
   end;
   end;
end;   

clear Q k kk v xx yy zz
drawnow;
warning on;
return;


