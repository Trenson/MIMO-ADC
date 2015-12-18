    function [snRdB,BER]=stepsize_choose
    clear;
    clc;
    N = 200; %number of antennas at the BS
    K = 50; %number of active users.

    modType = 'QPSK';
    deModType = 'LMMSE'; %'LMMSE';

    B_Bit1 = 9;
    B_Bit2 = 9; %7
    B_Bit3 = 5;
    
    S1 = N;
    S2 = 0;
    S3 = N-S1-S2;

    for snRdB = 0:2.5:15;
    Frame_Num = 10; %the number of simulation frame
    Frame_Len = 100; %length of each frame
    Q_StepSize = 0.001:0.001:0.07;
    BER=zeros(1,length(Q_StepSize));
    for isnr=1:length(Q_StepSize)
        s=[];
        s_out=[];
        for n_frame=1:Frame_Num %channel remains constant over each frame
           H=(randn(N,K)+1j*randn(N,K))*1/sqrt(2*K);
               % fprintf('-------running at frame %d/%d\n',n_frame,Frame_Num);
            for l_frame=1:Frame_Len %each frame contains Frame_Len symbols
                if strcmp(deModType,'GAMP')==1
                    [symbol,symbol_out]=Th_GAMP_Simu(K,N,H,snRdB,modType,Q_StepSize(isnr),B_Bit1,B_Bit2,B_Bit3,S1,S2,S3);
                else
                    [symbol,symbol_out]=Th_LMMSE_Simu_Step(K,N,H,snRdB,modType,Q_StepSize(isnr),B_Bit1,B_Bit2,B_Bit3,S1,S2,S3);
                end
                s=[s symbol];
                s_out=[s_out symbol_out];
            end
        end
        [~,BER(isnr)]=biterr(s,s_out);
    end
    hold on;
    [x,y]=min(BER);
    b=BER;
    [m,i]=min(b);
    b(i)=max(b);
    [m,i]=min(b);
    fprintf('-------ADC bit = %d ',B_Bit1);
    fprintf('-------SNR = %f ',snRdB);
    fprintf('-------Optimal step size = %f/%f\n',Q_StepSize(y),Q_StepSize(i));
    semilogy(Q_StepSize,BER,'-k');
    hold on;
    end
    %plot(snRdB,MSE,'b.')



