function val_o=myquantize(value,delta) 
q = round(value / delta);
if ((value - delta * q) <= (delta * (q + 1) - value))
    val_o=delta*q;
else
    val_o=delta*(q+1);
end    