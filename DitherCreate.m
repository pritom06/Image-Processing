function DVec=DitherCreate(delta,delta2,L)
dither=zeros(2,L);
dither(1,:)=rand(1,L)*delta-delta2;
for i=1:L
    if dither(1,i)<0
       dither(2,i)=dither(1,i)+delta2;
    else
       dither(2,i)=dither(1,i)-delta2;
    end
end 
DVec=dither;