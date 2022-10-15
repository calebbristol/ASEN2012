function [V1,dVdh] = get_Volume(h0,L)
% CODE CHALLENGE 5 - Given Function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT A SINGLE LINE IN THIS FUNCTION!!!!!! %
%   Take a peak inside to see how something like   %
% this can be done, but I will be using this       %
% function AS IS to test with your code.           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(h0 < 0.15) || any(h0 > 21.386) % checking if input is within domain
   error('h is out of allowed domain, h = [0.15, 21.386].\nCannot give accurate volume measurement');
end
h0 = 21.3858 - h0; %reversing where h is measured from

a = 0.019; %function coefficients
b = -5e-6;
c = 3.3e-10;
d = 3.4;
f = 0.0015;

dh = 1e-10; % dh for estimating dV/dh
for trial = 1:2 % two measurements will be taking 
    if trial == 2 %on second measurement, implement dh
        h0 = h0 - dh;
    end
h_x = @(x)h_of_x(x,h0,a,b,c,d,f); % function relating x position versus depth
x0 = 1000; % first initial guess point
x1 = fzero(h_x,x0); % zero finder
if h0 < 10 % if depth of interest is less than 10 ft
    x0 = 6000; % second initial point
else
    x0 = 3480; % second initial point
end
x2 = fzero(h_x,x0); % zero finder
    if trial == 1 % for first measurement
        V1 = V_of_x1x2(x1,x2,L,a,b,c,d,f,h0); %get Volume at requested depth
    else
        V2 = V_of_x1x2(x1,x2,L,a,b,c,d,f,h0); %get Volume dh away from requested depth
        dVdh = (V2 - V1)/dh; % get dVdh
    end
end

%% Subsubfunctions
function [h] = h_of_x(x,h0,a,b,c,d,f)
    h = a.*x+b.*x.^2+c.*x.^3+d.*sin(f.*x)-h0; %relate h with x
end
function [V] = V_of_x1x2(x1,x2,L,a,b,c,d,f,h0) %relate Vol with x1 and x2
    A = a.*x2^2/2+b.*x2.^3/3+c.*x2.^4/4-d./f.*cos(f.*x2);
    B = a.*x1^2/2+b.*x1.^3/3+c.*x1.^4/4-d./f.*cos(f.*x1);
    R = h0*abs(x2-x1);
    V = L*(A-B-R);
end

end