function [dhdt] = get_dhdt(h,L,alpha,dV_in)

dV_out = alpha*h; % calculate dV_out
dVdt = dV_in - dV_out; % calculate net dV/dt
[~,dVdh] = get_Volume(h,L) ; % get current dV/dh
dhdt = dVdt / dVdh; % convert dV/dt to dh/dt
end

