function [n_best, P_best, V, C]=runFoopsi(F,V,P)

% [n_best, P_best, V, C]=runFoopsi(F,V,P)
% Input:
% runFoopsi(F) - F - fluorescence trace, one cell at a time
% runFoopsi(F, V) - V - a structure of algorithm running parameters
% runFoopsi(F, dt) - dt - sampling interval of the trace
% runFoopsi(F, V, P) - P - calcium trace parameters (exp decay, scaling etc.)
% Outputs:
% n_best - inferred firing rate
% P_best - estimated P
% V - parameters used to run the algorithm
% C - inferred calcium trace (sum of exponential kernels).
%     F = P.a*C+P.b+noise

if ~isequal(class(F), 'double')
    F = double(F);
end

if nargin == 1
    V = struct;
    V.dt = 1;
end

if nargin >= 2
    if ~isstruct(V)
        dt = V;
        clear V;
        V = struct;
        V.dt = dt;
    elseif ~isfield(V, 'dt')
        V.dt = 1;
    end
end

if nargin < 3
    P = struct;
end

if ~isfield(P, 'gam')
    % if not defined, guess the exponential decay and noise level
    [gam, sig, ~] = getInitParams(F, V.dt);
    if gam<1 && gam>0
        P.gam = gam;
    end
    if isreal(sig) && sig>0
        P.sig = sig;
    end
end

% now, call the foopsi function
[n_best, P_best, V, C]=fast_oopsi(F,V,P);

