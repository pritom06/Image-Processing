function s = quantize_vector0(x,dither,delta1,n,K) 
x1=zeros(size(x));
for i = 1: n, 
    delta=delta1(i);
    x1(i) = myquantize(x(i)-K*dither(i), delta)+dither(i);
end
s=x1;
  
  
  
