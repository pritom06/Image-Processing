function sum=mydistance1(y,dither,delta1,n) 
sum1=0;
for i =1:n 
    delta=delta1(i);
    v=(y(i)-(myquantize(y(i)+dither(i),delta)-dither(i)));
    sum1=sum1+ v*v;
end
sum=sum1;  
 


