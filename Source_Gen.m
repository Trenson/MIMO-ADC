function [X,M]=Source_Gen(K,modType)
% Generate signal source --------------
% 
% Input: 
% K : number of users;
% modType = 'QPSK', '16QAM', 'Gaussian' modulation type;
% Outputs: 
% X : signal source

% Example: 
% X=Source_Gen(10,'QPSK')
% 
% ------------------------------------------------------------------
    switch modType    
        case 'QPSK'
            M=4;
        case '16QAM'
            M=16;
        case '64QAM' 
            M=64;            
    end
    Am=1/sqrt(2/3*(M-1)); %for M-QAM   
    symbol=randi(M,K,1)-1; %transmitted bits 
    x=qammod(symbol,M,0); %M-QAM modulation 
    X=x*Am; %the power normalization factor
end

