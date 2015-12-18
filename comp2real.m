function [Y_real,H_real,X_real]=comp2real(Y,H,X)
% complex value to real value
    H_real=[[real(H) -imag(H)];[imag(H) real(H)]];
    Y_real=[real(Y);imag(Y)];
    X_real=[real(X);imag(X)];
end
