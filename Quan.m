function YY_hat=Quan(YY,B_Bit,Delta)
% The quantization output
% -inf=rb(0)<rb(1)<...<rb(2^B)=+inf
% input :
%       YY : signal (L*1 vector)
%       B :  quantization bit
%       Delta : quantization step size
%
% usage :  YY_hat=Quan(YY,3,0.5)
    if B_Bit==0 % unquantized
        YY_hat=YY;
    else
        b=0:2^B_Bit;
        rb=(b-2^(B_Bit-1))'*Delta;
        k=ceil((YY-rb(2))/Delta);
        k(k<=0)=0;
        k(k>=2^B_Bit)=2^B_Bit-1;
        YY_hat=rb(k+2)-Delta/2;    
    end
end