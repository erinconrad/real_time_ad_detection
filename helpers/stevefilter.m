function [xout,zf] = stevefilter(xin,zi,alpha)
% STEVEFILTER  Causal high-pass-ish filter used throughout the AD pipeline.
% Removes a slow leaky-integrator estimate of the signal (one-pole),
% leaving the fast component. Maintains filter state zi across chunks for
% true real-time / streaming operation.
%
%   [xout,zf] = stevefilter(xin,zi)         uses default alpha = 0.99
%   [xout,zf] = stevefilter(xin,zi,alpha)
%
% zi/zf are the per-channel filter states (1 x nChannels).

if nargin < 3 || isempty(alpha)
    alpha = 0.99;
end

[xf,zf] = filter(1-alpha,[1, -alpha],xin,zi);
xout = xin - xf;

end
