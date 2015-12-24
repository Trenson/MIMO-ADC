%%This function computes kappa of the typical user given
%%N: number of BS antennas
%%Channel: a realization of the channel matrix of dimension M*N
%%K: number of high-resolution ADC pairs
%%SNR: transmit signal-to-noise ratio, defined as \mathcal{E}_s/1
%%Norm-based ADC assignment

function [kappa_multiuser_optimized]=rate_multiuser_optimized(N,Channel,K,SNR)
    Delta=zeros(N,1);
    Delta(N-K+1:N) = 1;
    Rrr=zeros(N,N);
    RrX=zeros(N,1);
    for index=1:N
        if Delta(index)==1
            RrX(index)=Channel(1,index)*SNR;
            %high resolution
        elseif Delta(index)==0
            RrX(index)=Channel(1,index)*SNR*sqrt(4/pi/(norm(Channel(:,index))^2*SNR+1));
            %one-bit
        end
        %RrX
        for subindex=index:N
            if index==subindex
                Rrr(index,subindex)=1+Delta(index)*norm(Channel(:,index))^2*SNR+1-Delta(index);%diagonal elements
            elseif (Delta(index)==1)&&(Delta(subindex)==1)
                Rrr(index,subindex)=Channel(:,index).'*conj(Channel(:,subindex))*SNR;
                Rrr(subindex,index)=Rrr(index,subindex)';%Hermitian matrix
            elseif (Delta(index)==1)&&(Delta(subindex)==0)
                Rrr(index,subindex)=Channel(:,index).'*conj(Channel(:,subindex))*SNR*sqrt(4/pi/(norm(Channel(:,subindex))^2*SNR+1));
                Rrr(subindex,index)=Rrr(index,subindex)';
            elseif (Delta(index)==0)&&(Delta(subindex)==1)
                Rrr(index,subindex)=Channel(:,index).'*conj(Channel(:,subindex))*SNR*sqrt(4/pi/(norm(Channel(:,index))^2*SNR+1));
                Rrr(subindex,index)=Rrr(index,subindex)';
            elseif (Delta(index)==0)&&(Delta(subindex)==0)
                Rrr(index,subindex)=4/pi*asin(real(Channel(:,index).'*conj(Channel(:,subindex)))*SNR/sqrt(norm(Channel(:,index))^2*SNR+1)/sqrt(norm(Channel(:,subindex))^2*SNR+1))+1i*4/pi*asin(imag(Channel(:,index).'*conj(Channel(:,subindex)))*SNR/sqrt(norm(Channel(:,index))^2*SNR+1)/sqrt(norm(Channel(:,subindex))^2*SNR+1));
                Rrr(subindex,index)=Rrr(index,subindex)';
            end
        end
        %Rrr
    end
    kappa_multiuser_optimized=real(RrX'*(Rrr\RrX))/SNR;
end

