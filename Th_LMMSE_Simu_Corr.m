function [S_out,S_in,MSE]= Th_LMMSE_Simu_Corr(K,N,H,snRdB,snrNo,modType,Q_StepSize,B_Bit1,B_Bit2,B_Bit3,S1,S2,S3)

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
    
%     H(:,K+1)=sum(abs(H').^2)';
%     H=sortrows(H,K+1);
%     H=H(:,1:K);
%     
%     y=zeros(1,K);
%     A=H(1:N-1,:);
%     x=H(N,:);
%     for i = 1:N-1
%         y = A(i,1:K);
%         corr = abs(corrcoef(x,y));
%         A(i,K+1) = corr(1,2);
%     end
%     A=sortrows(A,K+1);
%     A=A(:,1:K);
%     H(1:N-1,:)=A;
    max = 0;
    index = 0;
    H(:,K+1)=sum(abs(H').^2)';
    B = zeros(N - B_Bit2,K+1);
    A = H;
    cnt = 0;    
    while(size(A,1) > B_Bit2)
        max = 0;
        i = randi(size(A,1),1,1);
        if i > 1
            for j = 1:i-1
                x = A(i,1:K);
                y = A(j,1:K);
                temp = abs(corrcoef(x,y));
                corr = temp(1,2);
                if corr > max
                    max = corr;
                    index = j;
                end
            end
        end
        if i < size(A,1)
            for j = i+1:size(A,1)
                x = A(i,1:K);
                y = A(j,1:K);
                temp = abs(corrcoef(x,y));
                corr = temp(1,2);
                if corr > max
                    max = corr;
                    index = j;
                end
            end
        end
        if cnt < N - B_Bit2
            cnt = cnt + 1;
            if A(i,K+1) >= A(index,K+1)
                B(cnt,:) = A(index,:);
                A(index,:) = [];
            else
                B(cnt,:) = A(i,:);
                A(i,:) = [];
            end            
        end
    end
    H(1:N-B_Bit2,:) = B;
    H(N-B_Bit2+1:N,:) = A;
    H=H(:,1:K);
    
    Y= H*X+W;
%     save sortY Y;
    YY=[real(Y);imag(Y)];
%     save YY YY;
    
%     YY_hat1=Quan(YY,B_Bit1,Q_StepSize);
%     YY_hat2=Quan(YY,B_Bit2,Q_StepSize);
%     YY_hat = [YY_hat1(1:S1);YY_hat2(S1+1:N);YY_hat1(N+1:N+S1);YY_hat2(N+S1+1:2*N)];
    
    YY_hat1=Quan(YY,B_Bit1,Q_StepSize1);
    YY_hat2=Quan(YY,B_Bit2,Q_StepSize2);
    YY_hat3=Quan(YY,B_Bit3,Q_StepSize3);
    YY_hat = [YY_hat1(1:S1);YY_hat2(S1+1:S1+S2);YY_hat3(S1+S2+1:N);...
              YY_hat1(N+1:N+S1);YY_hat2(N+S1+1:N+S1+S2);YY_hat3(N+S1+S2+1:2*N)];
%     save YY_hat1 YY_hat1;
%     save YY_hat2 YY_hat2;
%     save YY_hat3 YY_hat3;
%     save YY_hat YY_hat;
    
    Y_hat=YY_hat(1:N)+1j*YY_hat(N+1:end);
    S_in=qamdemod(X,M,0);
    X_LMMSE=(H'*H+sigma2*eye(K))\H'*Y_hat; 
    MSE = norm([real(X);imag(X)]-[real(X_LMMSE);imag(X_LMMSE)],2)^2/(2*K);
    S_out=qamdemod(X_LMMSE,M,0);   % MMSE
end

    

