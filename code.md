# Improve-multi-population-genetic-algorithm-code
Improve multi-population genetic algorithm code
clear;
clc
NIND=40;               
NVAR=3;            
PRECI=10;             
GGAP=0.9;              
MP=2;            
FieldD=[rep(PRECI,[1,NVAR]);[100000,100000,100;2000000,60000000,2000000];rep([1;0;1;1],[1,NVAR])];  %译码矩阵
for i=1:MP
    Chrom{i}=crtbp(NIND, NVAR*PRECI);                  
end
pc=0.7+(0.9-0.7)*rand(MP,1);    
pm=0.001+(0.05-0.001)*rand(MP,1);  
gen=0; 
gen0=0; 
MAXGEN=15;  
maxY=0; 
for i=1:MP
    ObjV{i}=ObjectFunction(bs2rv(Chrom{i}, FieldD));
end
MaxObjV=-1*ones(MP,1);           
MaxChrom=zeros(MP,PRECI*NVAR); 
while gen0<=MAXGEN
    gen=gen+1;      
    for i=1:MP
        FitnV{i}=ranking(-ObjV{i});                      
        SelCh{i}=select('sus', Chrom{i}, FitnV{i},GGAP); 
        SelCh{i}=recombin('xovsp',SelCh{i}, pc(i));     
        SelCh{i}=mut(SelCh{i},pm(i));                  
        ObjVSel=ObjectFunction(bs2rv(SelCh{i}, FieldD)); 
        [Chrom{i},ObjV{i}]=reins(Chrom{i},SelCh{i},1,1,ObjV{i},ObjVSel);    
    end
    [Chrom,ObjV]=immigrantX(Chrom,ObjV);     
    MP=length(Chrom);  
for i=1:MP
    [MaxO,maxI]=max(ObjV{i}); 
    if MaxO>MaxObjV(i)
        MaxObjV(i)=MaxO;        
        MaxChrom(i,:)=Chrom{i}(maxI,:);  
    end
end                                              
    YY(gen)=max(MaxObjV);    
    if YY(gen)>maxY   
        maxY=YY(gen); 
        gen0=0;
    else
        gen0=gen0+1; 
    end
end
plot(1:gen,YY);
xlabel('Evolutionary algebra')
ylabel('Optimal solution change')
title('MPGA evolution process')
[Y,I]=max(MaxObjV);    
X=(bs2rv(MaxChrom(I,:), FieldD));   
disp(['The optimal value is：',num2str(Y)])
disp(['Corresponding independent variable value：',num2str(X)])
