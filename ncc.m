function nccV=ncc(wo,we) 
NCC1=and(wo,we);
NCCU=sum(sum(NCC1));
NCC2=and(wo,wo);
NCCD=sum(sum(NCC2));
nccV =NCCU/NCCD;