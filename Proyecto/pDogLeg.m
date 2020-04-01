function [p] = pDogLeg(B, g, delta)
    alfaU = dot(g,g) / dot(g, B * g);
    
    if alfaU >= delta / norm(g)
        p = - delta / norm(g) * g;
        break
    else
        pB = linsolve(-B, g);
        
        if norm(pB) <= delta
            p = pB;
            break
        else
            pU = -0.99 * alfaU * g;
            a = norm(pB - pU)^2;
            b = 2 * dot(pU, pB - pU);
            c = norm(pU)^2 - delta^2;
            alfa = (-b + sqrt(b^2 - 4*a*c)) / (2 * a);
            p = pU + alfa * (pB - pU);
end