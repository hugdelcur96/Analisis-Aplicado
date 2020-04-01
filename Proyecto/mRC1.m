function [x, msg, iter] = mRC1(f, x0, itmax)
    
    xk = x0;
    delta = 0.5;
    eta = 0.1;
    tol = 1e-5;
    deltaMax = 1.5;
    
    for k = 1:itmax
        fk = f(xk);
        gk = apGrad(f, xk);
        Bk = apHess(f, xk);
        mk = @(p) fk + dot(gk, p) + 0.5 * dot(p, Bk*p);
        pk = pCauchy(Bk, gk, delta);
        
        rhok = (f(xk) - f(xk + pk)) / (mk(zeros(length(pk), 1)) - mk(pk));
        
        if rhok < 0.25
            delta = 0.25 * delta;
        elseif rhok > 0.75 && norm(pk) == delta
            delta = min(2 * delta, deltaMax);
        end
        
        if rhok > eta
            xk = xk + pk;
        end
        
        if norm(gk) < tol
            x = xk;
            msg = 'Se encontro el minimo :)';
            iter = k;
            break
        elseif k == itmax
            x = xk;
            k = itmax;
            msg = 'No se encontro el minimo :(';
        end
    end
end
    