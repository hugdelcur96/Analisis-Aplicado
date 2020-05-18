function [alpha, gnew] = lineSearch(f, xk, dk, gk)
% In : f  ... objectve function (handle)
%      xk ... current point
%      dk ... chosen direction of descent
%      gk ... gradient of f in xk
%
% Out: alpha ... a parameter satisfying (W1) and (W2) 
%      gnew  ... gradient of f in xk+alpha*dk

    c1 = 1e-4;
    c2 = 0.99;
    alphaO = 0;
    alphaN = 1;
    alphaMax = 500;
    
    Dphi0 = dot(gk, dk);
    phi = @(alpha) f(xk + alpha * dk);
    L = @(alpha) phi(0) + c1 * alpha * Dphi0; 
    Dphi = @(alpha) dot(apGrad(f, xk + alpha * dk),dk);
    
    while alphaN > 0 && alphaN < alphaMax
        if phi(alphaN) > L(alphaN) || phi(alphaN) >= phi(alphaO)
            alphaAux = zoom(alphaO, alphaN, f, xk, gk, dk);
            break
        elseif abs(Dphi(alphaN)) <= -c2*Dphi0
            alphaAux = alphaN;
            break
        elseif Dphi(alphaN) >= 0
            alphaAux = zoom(alphaN, alphaO, f, xk, gk, dk);
            break
        else
            alphaO = alphaN;
            alphaN = 2 * alphaN;
        end
    end
    alpha = alphaAux;
    gnew = apGrad(f, xk + alpha * dk);
end


function [alpha] = zoom(alphaLo, alphaHi, f, xk, gk, dk)

    c1 = 1e-4;
    c2 = 0.99;

    Dphi0 = dot(gk, dk);
    phi = @(alpha) f(xk + alpha * dk);
    Dphi = @(alpha) dot(apGrad(f, xk + alpha * dk), dk);
    L = @(alpha) f(xk) + c1 * alpha * Dphi0;
    
    while 1
        alphaj = (alphaLo + alphaHi) / 2;
        if phi(alphaj) > L(alphaj) || phi(alphaj) >= phi(alphaLo)
            alphaHi = alphaj;
        elseif abs(Dphi(alphaj)) <= -c2 * Dphi0
            alphaAux = alphaj;
            break
        elseif Dphi(alphaj)*(alphaHi - alphaLo) >= 0
            alphaHi = alphaLo;
            alphaLo = alphaj;
        else
            alphaLo = alphaj;
        end
    end

    alpha = alphaAux;
end