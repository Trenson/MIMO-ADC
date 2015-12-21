function [S_out,S_in,MSE] = Th_GAMP_Simu_Sort(K,N,H,snRdB,snrNo,modType,Q_StepSize,B_Bit1,B_Bit2,B_Bit3,S1,S2,S3)
% [mse_X Pe] = Th_GAMP_Simu(K,N,snRdB,modType,B,Delta)
% Simulation result for MIMO with perfect CSI --------------
% MIMO channel:
%   y = (1/sqrt(K)) H x +  sigma^2 w 
% H ~ complex Gaussian with zero mean and variance 1
% w ~ complex Gaussian with zero mean and variance 1
% sigma^2 = 10^(-snRdB/10)
% 
% Input: 
% K : number of users;
% N : number of antennas;
% snRdB : snr in dB; 
% modType = 'QPSK', '16QAM', 'Gaussian' modulation type;
% bind is the pairs of high resolution ADCs
% Outputs: 
% mse_X : mean square error of X
% Pe : Bit error rate
% 
% Example: 
% [mse_X,Pe] = Th_GAMP_Simu(10,20,10,'QPSK')
% 
% ------------------------------------------------------------------

%    initial parameters
%     H=(randn(N,K)+1j*randn(N,K))*1/sqrt(2*K);
 
    sigma2=10^(-snRdB/10);
    Q(1,:) = [2.636, 2.271, 2.057, 1.934, 1.865, 1.825, 1.803]; % optimal step size of 1-bit
    Q(2,:) = [1.058, 0.954, 0.886, 0.844, 0.818, 0.804, 0.795]; % optimal step size of 2-bit
    Q(3,:) = [0.596, 0.531, 0.489, 0.464, 0.449, 0.440, 0.436]; % optimal step size of 3-bit
    Q(4,:) = [0.104, 0.092, 0.085, 0.080, 0.077, 0.076, 0.075]; % optimal step size of 6-bit
    Q_StepSize1 = Q(1,snrNo); % chooose the optimal step for 1-bit or 2-bit
    Q_StepSize2 = Q(4,snrNo); % chooose the optimal step for 1-bit or 2-bit
    Q_StepSize3 = Q(3,snrNo); % chooose the optimal step for 7-bit
    
    W=(randn(N,1)+1j*randn(N,1))*1/sqrt(2)*sqrt(sigma2);
    [X,M]=Source_Gen(K,modType);
    %sum(H')
    %H = reshape(sort(reshape(H,1,[]),'descend'),N,K);
    %sort
    H(:,K+1)=sum(abs(H').^2)';
    H=sortrows(H,K+1);
    H=H(:,1:K);
    
    Y= H*X+W;
    [YY,HH,XX]=comp2real(Y,H,X);
    YY_hat1=Quan(YY,B_Bit1,Q_StepSize1);
    YY_hat2=Quan(YY,B_Bit2,Q_StepSize2);
    YY_hat3=Quan(YY,B_Bit3,Q_StepSize3);
    YY = [YY_hat1(1:S1);YY_hat2(S1+1:S1+S2);YY_hat3(S1+S2+1:N);...
        YY_hat1(N+1:N+S1);YY_hat2(N+S1+1:N+S1+S2);YY_hat3(N+S1+S2+1:2*N)];
%     YY=YY_hat(1:N)+1j*YY_hat(N+1:end);
    t=1;tol=10^(-3);
    iter_max=20;
    x=zeros(2*K,1);
    vx=ones(2*K,1);
    s=zeros(2*N,1);
    B=0; % mismatch at the receiver
    while (t<iter_max)
        vp_new=(HH.*HH)*vx;
        p_new=HH*x-vp_new.*s;
        [s_new,vs_new]=OutputNonlinear(vp_new,p_new,YY,sigma2,B,Q_StepSize);
        vr_new=1./((HH.*HH)'*vs_new);
        r_new=x+vr_new.*(HH'*s_new);
        [x_new,vx_new]=InputNonlinear(r_new,vr_new,M);
        if (norm(x_new-x,2)/norm(x)<tol)
            x=x_new;
            break;
        else
            x=x_new;
            vx=vx_new;
            s=s_new;
            t=t+1;
        end    
    end
    mse_X=norm(XX-x,2)/norm(XX);
    if t>30
        fprintf('---iteration times-----: %d\n',t)
    end
    X_GAMP=x(1:K)+1j*x(K+1:end);
    S_out=qamdemod(X_GAMP,M,0);  
    S_in=qamdemod(X,M,0);
    MSE = norm([real(X);imag(X)]-[real(X_GAMP);imag(X_GAMP)],2)^2/(2*K);

%     x_MMSE=(H'*H+sigma2*eye(K))\H'*Y; 
%     S_out=qamdemod(x_MMSE,M,0);   % MMSE

end

function [x_new,vx_new]=InputNonlinear(r_new,vr_new,M)
% Input nonlinear step
% x=E2(x|r,vr) 
% vx=V2(x|r,vr) 

    Am=1/sqrt(2/3*(M-1));
    s=(1:2:(sqrt(M)-1))'*Am;
    p_pos=zeros(length(r_new),sqrt(M)/2);
    p_neg=zeros(length(r_new),sqrt(M)/2);
    ss=[s;-s];
    for i=1:sqrt(M)/2
        temp1=0;temp2=0;
        for j=1:sqrt(M)    
            temp1=temp1+exp(((s(i)-r_new).^2-(ss(j)-r_new).^2)./(2*vr_new));
            temp2=temp2+exp(((-s(i)-r_new).^2-(ss(j)-r_new).^2)./(2*vr_new));       
        end
        p_pos(:,i)=1./temp1;
        p_neg(:,i)=1./temp2;
    end
    x_new=(p_pos-p_neg)*s;
    vx_new=(p_pos+p_neg)*(s.^2)-x_new.^2;

% 	p0=1./(1+exp(sqrt(2)*r_new./vr_new));
%     p1=1-p0;	
%     x_new=1/sqrt(2)*(p1-p0); 
%     vx_new=0.5-x_new.^2; 
end

function [s_new,vs_new]=OutputNonlinear(vp_new,p_new,YY,sigma2,B,Delta)
% Output nonlinear step
% s=E1(y|p,vp+sigma2,Q) 
% vs=V1(y|p,vp+sigma2,Q) 
    if B==0
        v=vp_new+0.5*sigma2;
        s_new=(YY-p_new)./v;
        vs_new=1./v;
    else
        v=vp_new+0.5*sigma2;
        [a,b]=DeQuan(YY,B,Delta);
        stdd=sqrt(v);
        M=g(a,b,p_new,stdd);
        step_size=10^(-5);
        dM=(g(a,b,p_new+step_size,stdd)-g(a,b,p_new,stdd))/step_size;
        N=p_new.*M+v.*dM;
        E=N./M;
        dM2=(g(a,b,p_new+step_size,stdd)+g(a,b,p_new-step_size,stdd)-2*g(a,b,p_new,stdd))/(step_size^2);
        L=2*p_new.*N-(p_new.^2-v).*M+v.^2.*dM2;
        V=L./M-E.^2;
        s_new=(E-p_new)./v;
        vs_new=(1-V./v)./v;
    end
end

function y=g(a,b,mu,sigma)
% Gaussian variables with mean mu and variance sigma2
% y is the cumulative density probabilty(cdf) on interval [m,n]
    y=normcdf(b,mu,sigma)-normcdf(a,mu,sigma);
end
    
    

