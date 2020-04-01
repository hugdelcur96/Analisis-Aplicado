function [x, msg, iter] = mRC1(f, x0, itmax)
    
    x = x0;
    delta = 0.5;
    eta = 0.1;
    tol = 1e-5;
    deltaMax = 1.5;
    
    for k = 1:itmax
        fk = f(x);
        gk = apGrad(f, x);
        Bk = apHess(f, x);
        mk = @(p) fk + dot(gk, p) + 0.5 * dot(p, Bk*p);
        pk = pCauchy(Bk, gk, delta);
        
        rhok = (f(x) - f(x + pk)) / (mk(zeros(length(pk), 1)) - mk(pk));
        
        if rhok < 0.25
            delta = 0.25 * delta;
        elseif rhok > 0.75 && norm(pk) == delta
            delta = min(2 * delta, deltaMax);
        else
            delta = delta;
        end
        
        if rhok > eta
            x = x + pk;
        else
            x = x;
        end
        
        if norm(gk) < tol
            xk = x;
            msg = 'Se encontró el mínimo :)';
            iter = k;
            break
        else
            xk = x;
            msg = 'No se encontró el mínimo :(';
            iter = k;
        end
    end
end
    