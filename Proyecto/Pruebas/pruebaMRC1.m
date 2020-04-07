function [x, msg, iter] = pruebaMRC1(f, x0, itmax)
    xk = x0;
    delta = 0.5;
    eta = 0.1;
    tol = 1e-5;
    deltaMax = 1.5;
    
    k = 0;
    gk = apGrad(f, xk);
    
    while k < itmax && norm(gk) > tol
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
        k = k + 1;
    end
    
    x = xk;
    iter = k;
    
    if norm(gk) < tol
        msg = 'Se encontro el minimo :)';
    else
        msg = 'No se encontro el minimo :)';
end