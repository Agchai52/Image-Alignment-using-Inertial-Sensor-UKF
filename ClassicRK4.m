function X_New = ClassicRK4(f,X0,h)  
    F1 = h*f(X0);
    F2 = h*f(X0+F1/2);
    F3 = h*f(X0+F2/2);
    F4 = h*f(X0+F3);
    X_New = X0 + 1/6*(F1+2*F2+2*F3+F4);
end

