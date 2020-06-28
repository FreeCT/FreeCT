function f=ramp(k,ds,c,a)
% Simple function for creating filters
%
% Inputs:
%    k - list of indices (like [-200 -199 -198 ... 198 199])
%    ds - detector spacing (r_f*sin(fan_angle_increment/2))
%    c - [0,1]
%    a - [1, 0.5 0.54]
    
    f=(c^2/(2*ds))*(a*r(c*pi*k)+((1-a)/2)*r(pi*c*k+pi)+((1-a)/2)*r(pi*c*k-pi));
end

function vec=r(t)
vec=(sin(t)./t)+(cos(t)-1)./(t.^2);
vec(t==0)=0.5;
end