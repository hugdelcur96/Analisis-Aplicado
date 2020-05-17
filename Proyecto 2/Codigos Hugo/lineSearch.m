function [alpha, gnew] = lineSearch( f, xk, dk, gk )
% In : f  ... objectve function (handle)
%      xk ... current point
%      dk ... chosen direction of descent
%      gk ... gradient of f in xk
%
% Out: alpha ... a parameter satisfying (W1) and (W2) 
%      gnew  ... gradient of f in xk+alpha*dk
    
    c1 = 1e-4;
    c2 = 0.99;
    alpha0 = 0;
    alpha1 = 1;
    alphaMax = 500;
    
    phi = @(alpha) f(xk + alpha * dk);
    Dphi = @(alpha) dot(apGrad(f, xk + alpha * dk), dk);
    L = @(alpha) phi(0) + c1 * alpha * Dphi(0);
    
    while alpha1 < alphaMax
        if phi(alpha1) > L(alpha1) || phi(alpha1) >= phi(alpha0)
            alpha = zoom(alpha0, alpha1, f, xk, dk);
            break
        elseif abs(Dphi(alpha1)) <= -c2 * Dphi(0)
            alpha = alpha1;
            break
        elseif Dphi(alpha1) >= 0
            alpha = zoom(alpha1, alpha0, f, xk, dk);
            break
        else
            alpha0 = alpha1;
            alpha1 = 2 * alpha1;
        end
    end
    gnew = apGrad(f, xk + alpha * dk);
end

function [alpha] = zoom(alphaLo, alphaHi, f, xk, dk)
    c1 = 1e-4;
    c2 = 0.99;
    
    phi = @(alpha) f(xk + alpha * dk);
    Dphi = @(alpha) dot(apGrad(f, xk + alpha * dk), dk);
    L = @(alpha) phi(0) + c1 * alpha * Dphi(0);
    
    while 1
        alphaj = (alphaLo + alphaHi)/2; 
        
        if phi(alphaj) > L(alphaj) || phi(alphaj) >= phi(alphaLo)
            alphaHi = alphaj;
        elseif abs(Dphi(alphaj)) <= -c2 * Dphi(0)
            alpha = alphaj;
            break
        elseif Dphi(alphaj) * (alphaHi - alphaLo) >= 0
            alphaHi = alphaLo;
            alphaLo = alphaj;
        else
            alphaLo = alphaj;
        end
    end
end