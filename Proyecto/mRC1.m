function [xk, msg] = mRC1(f, x0, itmax)
    
    xk = x0
    deltak = 1
    eta = 0.1
    tol = 1e-5
    deltaMax = 1.5
    
    for k = 1:itmax
        fk = f(xk)
        gk = apGrad(f, xk)
        Bk = apHess(f, xk)
        mk = @(p) fk + dot(gk, p) + 0.5 * dot(p, Bk*p)
        pk = pCauchy(Bk, gk, deltak)
        
        rhok = (f(xk) - f(xk + pk)) / (mk(zeros(length(x0))') - mk(pk))
        
        if rhok < 0.25
            deltak = 0.25 * deltak
        elseif rhok > 0.75 && norm(pk) == deltak
            deltak = min(2 * deltak, deltamax)
        else
            deltak = deltak
        end
        
        if rhok > eta
            xk = xk + pk
        else
            xk = xk
        end

    end
    
    xk = xk
    msg = 'Se encontró el mínimo :)'
    
end
    