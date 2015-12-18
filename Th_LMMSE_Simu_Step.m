function [S_out,S_in,MSE]= Th_LMMSE_Simu_Step(K,N,H,snRdB,modType,Q_StepSize,B_Bit1,B_Bit2,B_Bit3,S1,S2,S3)
    sigma2=10^(-snRdB/10);
    W=(randn(N,1)+1j*randn(N,1))*1/sqrt(2)*sqrt(sigma2);
    [X,M]=Source_Gen(K,modType);
    Y= H*X+W;
    YY=[real(Y);imag(Y)];
    
%     YY_hat1=Quan(YY,B_Bit1,Q_StepSize);
%     YY_hat2=Quan(YY,B_Bit2,Q_StepSize);
%     YY_hat = [YY_hat1(1:S1);YY_hat2(S1+1:N);YY_hat1(N+1:N+S1);YY_hat2(N+S1+1:2*N)];
    
    YY_hat1=Quan(YY,B_Bit1,Q_StepSize);
    YY_hat2=Quan(YY,B_Bit2,Q_StepSize);
    YY_hat3=Quan(YY,B_Bit3,Q_StepSize);
    YY_hat = [YY_hat1(1:S1);YY_hat2(S1+1:S1+S2);YY_hat3(S1+S2+1:N);...
              YY_hat1(N+1:N+S1);YY_hat2(N+S1+1:N+S1+S2);YY_hat3(N+S1+S2+1:2*N)];
    
    Y_hat=YY_hat(1:N)+1j*YY_hat(N+1:end);
    S_in=qamdemod(X,M,0);
    X_LMMSE=(H'*H+sigma2*eye(K))\H'*Y_hat; 
    MSE = norm([real(X);imag(X)]-[real(X_LMMSE);imag(X_LMMSE)],2)^2/(2*K);
    S_out=qamdemod(X_LMMSE,M,0);   % MMSE
end

    

