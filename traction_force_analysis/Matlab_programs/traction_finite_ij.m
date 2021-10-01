%function [tzx, tzy] = traction_finite_ij(uin,vin,h,d,sig,E)
function [tzx, tzy, Tzx, Tzy] = traction_finite_ij(uin,vin,h,d,sig,E)
% TRACTION_FINITE returns a 2D traction field from a 2D displacement field.
% 
% tzx and tzy are the tractions in x and y directions
% uin and vin are the displacement fields in x and y directions
% h is the z position at which the tzx and tzy are computed (usually the height of the gel)
% d is the size of one pixel of the PIV analysis (distance between 2
% displacement points)
% sig is the poisson modulus
% E is the Young's modulus

%
% UNITS:    uin, vin, h and d should be in microns
%           E in Pascal
%           
%
% The program implements the solution provided by Del Alamo et al (PNAS, 2007)
% A number of typos in Del Alamo et al have been identified and fixed.
%
% Xavier Trepat and Jim Butler 09/2008
%

Nalpha = size(uin,2); 
Nbeta = size(uin,1); 
s1=1-sig; % useful factors that appear later
s2=1-2*sig;
s34=3-4*sig;

U = fft2(uin).'; % Fourier transform the displacement field
V = fft2(vin).';

for kalpha=1:Nalpha, 
  %  fprintf('.')
    for kbeta=1:Nbeta,        
        
        if kalpha <= Nalpha/2, alpha = (kalpha-1)*2*pi/(Nalpha*d);
        else alpha =(-Nalpha-1+kalpha)*2*pi/(Nalpha*d);
        end
        
        if kbeta <=  Nbeta/2, beta = (kbeta-1)*2*pi/(Nbeta*d);
        else beta =(-Nbeta-1+kbeta)*2*pi/(Nbeta*d);
        end
        
        k=sqrt(alpha^2+beta^2);
        if k==0, % DC
            Tzx(1,1)=E/(2*(1+sig))*U(1,1)/h;
            Tzy(1,1)=E/(2*(1+sig))*V(1,1)/h;
        else
            uh0=U(kalpha,kbeta); % these are the fourier coeficients of the displacement field
            vh0=V(kalpha,kbeta);
            
 
            Tzx(kalpha,kbeta) = -E*beta *cosh(k*h)/(2*(1+sig)*k*sinh(k*h))*(vh0*alpha-uh0*beta) ...
                + E*alpha/(2*(1-sig^2)*k)*((s34*cosh(k*h)^2)+s2^2+(k*h)^2)/(s34*sinh(k*h)*cosh(k*h)+k*h)*(alpha*uh0+beta*vh0);
            
            Tzy(kalpha,kbeta) = +E*alpha*cosh(k*h)/(2*(1+sig)*k*sinh(k*h))*(vh0*alpha-uh0*beta) ...
                + E*beta /(2*(1-sig^2)*k)*((s34*cosh(k*h)^2)+s2^2+(k*h)^2)/(s34*sinh(k*h)*cosh(k*h)+k*h)*(alpha*uh0+beta*vh0);
               
         
        end % if k==0;
    end %for kalpha=1:Nalpha,
end % for kbeta=1:Nbeta,

% invert
tzx=real(ifft2(Tzx).');
tzy=real(ifft2(Tzy).');
