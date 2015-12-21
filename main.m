    function [snRdB,BER]=main
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
    
    S1 = N*0.7;
    S2 = N*0.15;
    S3 = N-S1-S2;

    Frame_Num = 1000; %the number of simulation frame
    Frame_Len = 100; %length of each frame
    snRdB = 0:2.5:15;
    BER=zeros(1,length(snRdB));
    for isnr=1:length(snRdB)
        s=[];
        s_out=[];
        fprintf('[%d/%d] MIMO configure(N¡ÁK) : %3d¡Á%3d  running at SNR = %3d (dB)\n',isnr,length(snRdB),N,K,snRdB(isnr));
        MSE(isnr) = 0;
        MSE2(isnr) = 0;
        s2=[];
        s_out2=[];
        for n_frame=1:Frame_Num %channel remains constant over each frame
           H=(randn(N,K)+1j*randn(N,K))*1/sqrt(2*K);
                fprintf('-------running at frame %d/%d\n',n_frame,Frame_Num);
            for l_frame=1:Frame_Len %each frame contains Frame_Len symbols
                if strcmp(deModType,'GAMP')==1
                    [symbol,symbol_out,MSE_temp]=Th_GAMP_Simu(K,N,H,snRdB(isnr),isnr,modType,Q_StepSize,B_Bit1,B_Bit2,B_Bit3,S1,S2,S3);
                    [symbol2,symbol_out2,MSE_temp2]=Th_GAMP_Simu_Sort(K,N,H,snRdB(isnr),isnr,modType,Q_StepSize,B_Bit1,B_Bit2,B_Bit3,S1,S2,S3);
                else
                    [symbol,symbol_out,MSE_temp]=Th_LMMSE_Simu(K,N,H,snRdB(isnr),isnr,modType,Q_StepSize,B_Bit1,B_Bit2,B_Bit3,S1,S2,S3);
                    [symbol2,symbol_out2,MSE_temp2]=Th_LMMSE_Simu_Sort(K,N,H,snRdB(isnr),isnr,modType,Q_StepSize,B_Bit1,B_Bit2,B_Bit3,S1,S2,S3);
                end
                MSE(isnr) = MSE(isnr)+MSE_temp;               
                s=[s symbol];
                s_out=[s_out symbol_out];
                MSE2(isnr) = MSE2(isnr)+MSE_temp2;
                s2=[s2 symbol2];
                s_out2=[s_out2 symbol_out2];
            end
        end
        [~,BER(isnr)]=biterr(s,s_out);
        [~,BER2(isnr)]=biterr(s2,s_out2);
    end
    MSE = MSE/(Frame_Num*Frame_Len);
    MSE2 = MSE2/(Frame_Num*Frame_Len);
    hold on;
    BER
    BER2
    semilogy(snRdB,BER,'-k');
    hold on;
    semilogy(snRdB,BER2,'-k');
    %plot(snRdB,MSE,'b.')



