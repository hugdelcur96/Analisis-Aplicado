function [pC] = pCauchy(B, g, delta)
% In:   B     ... (symmetric matrix) approximates the hessian of f in xk
%       g     ... (vector) gradient of f in xk
%       delta ... trust region radius
% 
% Out:  pC    ... The Cauchy point

    if dot(g, B * g) <= 0
        tau = 1;
    else
        tau = min([norm(g)^3 / (delta * dot(g, B * g)), 1]);
    end
    
    pC = -tau * delta / norm(g) * g;
end