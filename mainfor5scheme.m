    function [snRdB,BER]=mainfor5scheme
    clear;
    clc;
    N = 200; %number of antennas at the BS
    K = 50; %number of active users.

    Q_StepSize = 0.5; %0.5
    modType = 'QPSK';
    deModType = 'LMMSE'; %'LMMSE';

    B_Bit1 = 1;
    B_Bit2 = 7; %7
    B_Bit3 = 5;
    
    S1 = N*0.8;
    S2 = N*0.2;
    S3 = N-S1-S2;

    Frame_Num = 1; %the number of simulation frame
    Frame_Len = 100; %length of each frame
    snRdB = 0:2.5:15;
    BER=zeros(1,length(snRdB));
    for isnr=1:length(snRdB)
        s=[];
        s_out=[];
        fprintf('[%d/%d] MIMO configure(N¡ÁK) : %3d¡Á%3d  running at SNR = %3d (dB)\n',isnr,length(snRdB),N,K,snRdB(isnr));
        MSE(isnr) = 0;
        MSE2(isnr) = 0;
        MSE3(isnr) = 0;
        MSE4(isnr) = 0;
        MSE5(isnr) = 0;
        s2=[];
        s_out2=[];
        s3=[];
        s_out3=[];
        s4=[];
        s_out4=[];
        s5=[];
        s_out5=[];
        for n_frame=1:Frame_Num %channel remains constant over each frame
           H=(randn(N,K)+1j*randn(N,K))*1/sqrt(2*K);
                fprintf('-------running at frame %d/%d\n',n_frame,Frame_Num);
            for l_frame=1:Frame_Len %each frame contains Frame_Len symbols
                if strcmp(deModType,'GAMP')==1
                    [symbol,symbol_out,MSE_temp]=Th_GAMP_Simu(K,N,H,snRdB(isnr),isnr,modType,Q_StepSize,B_Bit1,B_Bit2,B_Bit3,S1,S2,S3);
                    [symbol2,symbol_out2,MSE_temp2]=Th_GAMP_Simu_Sort(K,N,H,snRdB(isnr),isnr,modType,Q_StepSize,B_Bit1,B_Bit2,B_Bit3,S1,S2,S3);
                    [symbol3,symbol_out3,MSE_temp3]=Th_GAMP_Simu_Corr(K,N,H,snRdB(isnr),isnr,modType,Q_StepSize,B_Bit1,B_Bit2,B_Bit3,S1,S2,S3);
                    [symbol4,symbol_out4,MSE_temp4]=Th_GAMP_Simu_Value(K,N,H,snRdB(isnr),isnr,modType,Q_StepSize,B_Bit1,B_Bit2,B_Bit3,S1,S2,S3);
                    [symbol5,symbol_out5,MSE_temp5]=Th_GAMP_Simu_Det(K,N,H,snRdB(isnr),isnr,modType,Q_StepSize,B_Bit1,B_Bit2,B_Bit3,S1,S2,S3);
                else
                    [symbol,symbol_out,MSE_temp]=Th_LMMSE_Simu(K,N,H,snRdB(isnr),isnr,modType,Q_StepSize,B_Bit1,B_Bit2,B_Bit3,S1,S2,S3);
                    [symbol2,symbol_out2,MSE_temp2]=Th_LMMSE_Simu_Sort(K,N,H,snRdB(isnr),isnr,modType,Q_StepSize,B_Bit1,B_Bit2,B_Bit3,S1,S2,S3);
                    [symbol3,symbol_out3,MSE_temp3]=Th_LMMSE_Simu_Corr(K,N,H,snRdB(isnr),isnr,modType,Q_StepSize,B_Bit1,B_Bit2,B_Bit3,S1,S2,S3);
                    [symbol4,symbol_out4,MSE_temp4]=Th_LMMSE_Simu_Sort_De(K,N,H,snRdB(isnr),isnr,modType,Q_StepSize,B_Bit1,B_Bit2,B_Bit3,S1,S2,S3);
                    [symbol5,symbol_out5,MSE_temp5]=Th_LMMSE_Simu_Det(K,N,H,snRdB(isnr),isnr,modType,Q_StepSize,B_Bit1,B_Bit2,B_Bit3,S1,S2,S3);
                end
                MSE(isnr) = MSE(isnr)+MSE_temp;               
                s=[s symbol];
                s_out=[s_out symbol_out];
                MSE2(isnr) = MSE2(isnr)+MSE_temp2;
                s2=[s2 symbol2];
                s_out2=[s_out2 symbol_out2];
                MSE3(isnr) = MSE3(isnr)+MSE_temp3;               
                s3=[s3 symbol3];
                s_out3=[s_out3 symbol_out3];
                MSE4(isnr) = MSE4(isnr)+MSE_temp4;
                s4=[s4 symbol4];
                s_out4=[s_out4 symbol_out4];
                MSE5(isnr) = MSE5(isnr)+MSE_temp5;               
                s5=[s5 symbol5];
                s_out5=[s_out5 symbol_out5];
            end
        end
        [~,BER(isnr)]=biterr(s,s_out);
        [~,BER2(isnr)]=biterr(s2,s_out2);
        [~,BER3(isnr)]=biterr(s3,s_out3);
        [~,BER4(isnr)]=biterr(s4,s_out4);
        [~,BER5(isnr)]=biterr(s5,s_out5);
    end
    MSE = MSE/(Frame_Num*Frame_Len);
    MSE2 = MSE2/(Frame_Num*Frame_Len);
    MSE3 = MSE3/(Frame_Num*Frame_Len);
    MSE4 = MSE4/(Frame_Num*Frame_Len);
    MSE5 = MSE5/(Frame_Num*Frame_Len);
    hold on;
    BER
    BER2
    BER3
    BER4
    BER5
    semilogy(snRdB,BER,'-k');
    hold on;
    semilogy(snRdB,BER2,'-k');
    hold on;
    semilogy(snRdB,BER3,'-k');
    hold on;
    semilogy(snRdB,BER4,'-k');
    hold on;
    semilogy(snRdB,BER5,'-k');
    %plot(snRdB,MSE,'b.')



