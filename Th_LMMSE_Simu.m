function [S_out,S_in,MSE,GMI,REF]= Th_LMMSE_Simu(K,N,H,snRdB,snrNo,modType,Q_StepSize,B_Bit1,B_Bit2,B_Bit3,S1,S2,S3)

    Q(1,:) = [2.638, 1.925, 1.519, 1.277, 1.131, 1.043, 0.990]; % optimal step size of 1-bit
    Q(2,:) = [0.992, 0.874, 0.801, 0.759, 0.735, 0.721, 0.713]; % optimal step size of 2-bit
    Q(3,:) = [0.583, 0.514, 0.475, 0.448, 0.432, 0.423, 0.419]; % optimal step size of 3-bit
    Q(4,:) = [0.187, 0.165, 0.152, 0.145, 0.145, 0.145, 0.145]; % optimal step size of 5-bit
    Q(5,:) = [0.103, 0.092, 0.084, 0.079, 0.077, 0.075, 0.075]; % optimal step size of 6-bit
    Q(6,:) = [0.039, 0.036, 0.034, 0.006, 0.004, 0.003, 0.003]; % optimal step size of 9-bit
    Q_StepSize1 = Q(1,snrNo); % chooose the optimal step for 1-bit or 2-bit
    Q_StepSize2 = Q(5,snrNo); % chooose the optimal step for 7-bit
    Q_StepSize3 = Q(3,snrNo); % chooose the optimal step for 9-bit
    sigma2=10^(-snRdB/10);
    W=(randn(N,1)+1j*randn(N,1))*1/sqrt(2)*sqrt(sigma2);
    [X,M]=Source_Gen(K,modType);
    Y= H*X+W;
    YY=[real(Y);imag(Y)];
    
%     Rrr = zeros(N,N);
%     Rrx = zeros(N,1);
%     delta = zeros(N,1);
%     delta(N-S2+1:N) = 1;
%     es = 1;
%     for n = 1 : N
%         for m = 1 : N
%             if n == m
%                 Rrr(n,m) = 1 + delta(n)*norm(H(n,:))^2*es + 1 - delta(n);
%             else
%                 Rrr(n,m) = H(n,:)*H(m,:)'*es*(delta(n)*delta(m) + delta(n)*(1 - delta(m))*sqrt(4/(pi*(norm(H(m,:))^2*es+1))) + ...
%                 (1 - delta(n))*delta(m)*sqrt(4/(pi*(norm(H(n,:))^2*es+1))) + (1 - delta(n))*(1 - delta(m))*4/pi*(...
%                 asin(real(H(n,:)*H(m,:)')*es/(sqrt(norm(H(n,:))^2*es+1)*sqrt(norm(H(m,:))^2*es+1))) + ...
%                 1j*asin(imag(H(n,:)*H(m,:)')*es/(sqrt(norm(H(n,:))^2*es+1)*sqrt(norm(H(m,:))^2*es+1)))));
%             end            
%         end       
%     end
%     save Rrr Rrr; 
%     GMI = 0;
%     for u = 1:K
%         for n = 1:N
%             Rrx(n) = H(n,u)*es*(delta(n) + (1 - delta(n))*sqrt(4/(pi*(norm(H(n,:))^2*es+1))));
%         end
%         k = real(Rrx'*(Rrr\Rrx))/es;
%         GMI = GMI + log(1 + k/(1-k))/log(2);
%     end
%     for n = 1 :N
%         sumH = 0;
%         for l = 1: K
%             sumH = sumH + H(n,l)*X(l,1);
%         end
%         r(n,1) = delta(n)*(sumH + W(n)) + (1 - delta(n))*sign(sumH + W(n));
%     end
%     cov(r,r)
%     GMI = GMI/K;
    es = 10.^(snRdB/10)/K;
    [kappa_multiuser_optimized]=rate_multiuser_optimized(N, H.', S2, es);
    GMI = log(1 + kappa_multiuser_optimized/(1-kappa_multiuser_optimized))/log(2);
    REF = real(log(det(eye(K)+es*H.'*conj(H)))/log(2)/K);

%       for n = 1:N
%           Rrx(n,1) = H(n,1)*es*(delta(n) + (1 - delta(n))*sqrt(4/(pi*(sum(abs(H(n,:)).^2)*es+1))));
%       end
%       k = real(1*Rrx'*inv(Rrr)*Rrx);
%       GMI = log(1 + k/(1-k))/log(2);
    
%     YY_hat1=Quan(YY,B_Bit1,Q_StepSize);
%     YY_hat2=Quan(YY,B_Bit2,Q_StepSize);
%     YY_hat = [YY_hat1(1:S1);YY_hat2(S1+1:N);YY_hat1(N+1:N+S1);YY_hat2(N+S1+1:2*N)];
    
    YY_hat1=Quan(YY,B_Bit1,Q_StepSize1);
    YY_hat2=Quan(YY,B_Bit2,Q_StepSize2);
    YY_hat3=Quan(YY,B_Bit3,Q_StepSize3);
    YY_hat = [YY_hat1(1:S1);YY_hat2(S1+1:S1+S2);YY_hat3(S1+S2+1:N);...
              YY_hat1(N+1:N+S1);YY_hat2(N+S1+1:N+S1+S2);YY_hat3(N+S1+S2+1:2*N)];
    
    Y_hat=YY_hat(1:N)+1j*YY_hat(N+1:end);
    S_in=qamdemod(X,M,0);
    X_LMMSE=(H'*H+sigma2*eye(K))\H'*Y_hat; 
    MSE = norm([real(X);imag(X)]-[real(X_LMMSE);imag(X_LMMSE)],2)^2/(2*K);
    S_out=qamdemod(X_LMMSE,M,0);   % MMSE
end

    

