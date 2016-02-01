function [gam, sig, fitRes] = getInitParams(F, dt)

% this function will provide a reasonable guess for the exponential decay
% constant to be used in foopsi (Vogelstein's deconvolution algorithm)
% it also estimates the Gaussian noise variance (sig).
%
% Works under the assumption that the calcium trace is a convolution of
% Poisson spike train with exponential decay kernel with additive
% Gaussian White Noise.

% selecting nLags is tricky, should be done adaptively, especially when
% taking into account that we have different GCaMPs
if nargin>1
    nLags = max(30, ceil(5/dt)-1); % at least 5 seconds, or 30 samples
else
    nLags = 30;
    dt = 1;
end
xc = xcov(F, nLags, 'unbiased');
% xc = xc-xc(end);
tauAxis = dt*(-nLags:nLags)';
zeroLag = nLags+1;
dLag = 1;
myFun = @(a, b, c, x) a*exp(b*x)+c;

% taking the first 10 points to avoid noisy data
lags = zeroLag+dLag:zeroLag+dLag+10;
% 'guessing' some IC for the exponential fit
b0 = log(mean(xc(lags+1)./xc(lags)))/dt;
a0 = xc(zeroLag+dLag)*exp(-b0*dLag*dt);
% this is the initial fit, to get a good guess for the initial conditions
% for the full fit. Should be quite fast, because is optimized (analytical solution?)
weights = linspace(1, 0, nLags-dLag+1);
f = fit(tauAxis(zeroLag+dLag:end), double(xc(zeroLag+dLag:end)-xc(end)), 'exp1', 'StartPoint', [a0, b0], 'Weights', weights);
% and this one is the full fit
[f, gof] = fit(tauAxis(zeroLag+dLag:end), double(xc(zeroLag+dLag:end)), myFun, 'StartPoint', [f.a, f.b, xc(end)], 'Weights', weights);

gam = exp(f.b*dt);
sig = sqrt(xc(zeroLag)-f.a-f.c);

fitRes.tauAxis = tauAxis;
fitRes.xc = xc;
fitRes.f = f;
fitRes.gof = gof;

return
%%
figure


plot(tauAxis, xc, 'b.');
hold on;
plot(tauAxis(zeroLag:end), f.a*exp(f.b*tauAxis(zeroLag:end))+f.c, 'r');
