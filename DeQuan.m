function [m,n]=DeQuan(YY,B,Delta)
% The de-quantization function
% input :
%       YY : quantized signal (L*1 vector)
%       B :  quantization bit
%       Delta : quantization step size
% output:
%       YY falls into the interval [a,b]
% usage :  [a,b]=Quan(YY,3,0.5)
    
    b=(1:2^B-1)';
    rb=(-2^(B-1)+b)*Delta;
    lowbound=[-inf;rb];
    upbound=[rb;inf];

    k=ceil((YY-upbound(1))/Delta);
    k(k>=2^B-1)=2^B-1;
    k(k<=0)=0;
    m=lowbound(k+1);
    n=upbound(k+1);

