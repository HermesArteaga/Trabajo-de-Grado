clear all
clc
Cadena = [load('Vars_LQR_d1'),load('Vars_LQR_d2'),load('Vars_LQR_d3'),load('Vars_LQR_i1'),load('Vars_LQR_i2')];
Prom=zeros(5,4);
for j=1:5
    Vars_s=Cadena(1,j);
    Vars=Vars_s.Vars_LQR;    
    k = find(Vars(:,4)<=0.349&Vars(:,4)>=-0.349);
    M = Vars(k(1,1):3001,1:6);
    Media(j,:)=mean(abs(M));
    Des_Est(j,:)=std(abs(M));
    t = M(:,1);    
    for i=2:6
        M1 = M(:,i);
        IND((i-1)+((j-1)*5),1)=trapz(t,M1.^2); %ISE
        IND((i-1)+((j-1)*5),2)=trapz(t,abs(M1)); %IAE
        IND((i-1)+((j-1)*5),3)=trapz(t,t.*(M1.^2)); %ITSE
        IND((i-1)+((j-1)*5),4)=trapz(t,t.*(abs(M1))); %ITAE
        Prom((i-1),1)=IND((i-1)+((j-1)*5),1)+Prom((i-1),1);        
        Prom((i-1),2)=IND((i-1)+((j-1)*5),2)+Prom((i-1),2);
        Prom((i-1),3)=IND((i-1)+((j-1)*5),3)+Prom((i-1),3);        
        Prom((i-1),4)=IND((i-1)+((j-1)*5),4)+Prom((i-1),4);
    end       
end
IND;
Prom_lqr=Prom./5;
disp('LQR');
Media
Des_Est;
%%
Cadena = [load('Vars_LQRI_i1'),load('Vars_LQRI_i2'),load('Vars_LQRI_i3'),load('Vars_LQRI_i4'),load('Vars_LQRI_i5')];
Prom=zeros(5,4);
for j=1:5
    Vars_s=Cadena(1,j);
    Vars=Vars_s.Vars_LQRI;
    k = find(Vars(:,4)<=0.349&Vars(:,4)>=-0.349);
    M = Vars(k(1,1):3001,1:6);
    Media(j,:)=mean(abs(M));
    Des_Est(j,:)=std(abs(M));
    t = M(:,1);
    for i=2:6
        M1 = M(:,i);
        IND((i-1)+((j-1)*5),1)=trapz(t,M1.^2); %ISE
        IND((i-1)+((j-1)*5),2)=trapz(t,abs(M1)); %IAE
        IND((i-1)+((j-1)*5),3)=trapz(t,t.*(M1.^2)); %ITSE
        IND((i-1)+((j-1)*5),4)=trapz(t,t.*(abs(M1))); %ITAE
        Prom((i-1),1)=IND((i-1)+((j-1)*5),1)+Prom((i-1),1);        
        Prom((i-1),2)=IND((i-1)+((j-1)*5),2)+Prom((i-1),2);
        Prom((i-1),3)=IND((i-1)+((j-1)*5),3)+Prom((i-1),3);        
        Prom((i-1),4)=IND((i-1)+((j-1)*5),4)+Prom((i-1),4);
    end       
end
IND;
Prom_lqri=Prom./5;
disp('LQR+i');
Media
Des_Est;
%%
Cadena = [load('Vars_Gain_d1'),load('Vars_Gain_d2'),load('Vars_Gain_i1'),load('Vars_Gain_i2'),load('Vars_Gain_i3')];
Prom=zeros(5,4);
for j=1:5
    Vars_s=Cadena(1,j);
    Vars=Vars_s.Vars_Gain;
    k = find(Vars(:,4)<=0.349&Vars(:,4)>=-0.349);
    M = Vars(k(1,1):3001,1:6);
    Media(j,:)=mean(M);
    Des_Est(j,:)=std(M);
    t = M(:,1);
    for i=2:6
        M1 = M(:,i);
        IND((i-1)+((j-1)*5),1)=trapz(t,M1.^2); %ISE
        IND((i-1)+((j-1)*5),2)=trapz(t,abs(M1)); %IAE
        IND((i-1)+((j-1)*5),3)=trapz(t,t.*(M1.^2)); %ITSE
        IND((i-1)+((j-1)*5),4)=trapz(t,t.*(abs(M1))); %ITAE
        Prom((i-1),1)=IND((i-1)+((j-1)*5),1)+Prom((i-1),1);        
        Prom((i-1),2)=IND((i-1)+((j-1)*5),2)+Prom((i-1),2);
        Prom((i-1),3)=IND((i-1)+((j-1)*5),3)+Prom((i-1),3);        
        Prom((i-1),4)=IND((i-1)+((j-1)*5),4)+Prom((i-1),4);
    end       
end
IND;
Prom_gain=Prom./5;
disp('Gain');
Media
Des_Est;
%%
Cadena = [load('Vars_Gaini_d1'),load('Vars_Gaini_d2'),load('Vars_Gaini_d3'),load('Vars_Gaini_i1'),load('Vars_Gaini_i2')];
Prom=zeros(5,4);
for j=1:5
    Vars_s=Cadena(1,j);
    Vars=Vars_s.Vars_Gain_i;
    k = find(Vars(:,4)<=0.349&Vars(:,4)>=-0.349);
    M = Vars(k(1,1):3001,1:6);
    Media(j,:)=mean(M);
    Des_Est(j,:)=std(M);
    t = M(:,1);
    for i=2:6
        M1 = M(:,i);
        IND((i-1)+((j-1)*5),1)=trapz(t,M1.^2); %ISE
        IND((i-1)+((j-1)*5),2)=trapz(t,abs(M1)); %IAE
        IND((i-1)+((j-1)*5),3)=trapz(t,t.*(M1.^2)); %ITSE
        IND((i-1)+((j-1)*5),4)=trapz(t,t.*(abs(M1))); %ITAE
        Prom((i-1),1)=IND((i-1)+((j-1)*5),1)+Prom((i-1),1);        
        Prom((i-1),2)=IND((i-1)+((j-1)*5),2)+Prom((i-1),2);
        Prom((i-1),3)=IND((i-1)+((j-1)*5),3)+Prom((i-1),3);        
        Prom((i-1),4)=IND((i-1)+((j-1)*5),4)+Prom((i-1),4);
    end       
end
IND;
Prom_gaini=Prom./5;
disp('Gain_i');
Media
Des_Est;
%%
Cadena = [load('Vars_Fuzzy_d1'),load('Vars_Fuzzy_d2'),load('Vars_Fuzzy_i1'),load('Vars_Fuzzy_i2'),load('Vars_Fuzzy_i3')];
Prom=zeros(5,4);
for j=1:5
    Vars_s=Cadena(1,j);
    Vars=Vars_s.Vars_Fuzzy;
    k = find(Vars(:,4)<=0.349&Vars(:,4)>=-0.349);
    M = Vars(k(1,1):3001,1:6);
    Media(j,:)=mean(M);
    Des_Est(j,:)=std(M);
    t = M(:,1);
    for i=2:6
        M1 = M(:,i);
        IND((i-1)+((j-1)*5),1)=trapz(t,M1.^2); %ISE
        IND((i-1)+((j-1)*5),2)=trapz(t,abs(M1)); %IAE
        IND((i-1)+((j-1)*5),3)=trapz(t,t.*(M1.^2)); %ITSE
        IND((i-1)+((j-1)*5),4)=trapz(t,t.*(abs(M1))); %ITAE
        Prom((i-1),1)=IND((i-1)+((j-1)*5),1)+Prom((i-1),1);        
        Prom((i-1),2)=IND((i-1)+((j-1)*5),2)+Prom((i-1),2);
        Prom((i-1),3)=IND((i-1)+((j-1)*5),3)+Prom((i-1),3);        
        Prom((i-1),4)=IND((i-1)+((j-1)*5),4)+Prom((i-1),4);
    end       
end
IND;
Prom_Fuzzy=Prom./5;
disp('Fuzzy_LQR');
Media
Des_Est;
%%
Cadena = [load('Vars_Fuzzyg_d1'),load('Vars_Fuzzyg_d2'),load('Vars_Fuzzyg_d3'),load('Vars_Fuzzyg_i1'),load('Vars_Fuzzyg_i2')];
Prom=zeros(5,4);
for j=1:5
    Vars_s=Cadena(1,j);
    Vars=Vars_s.Vars_Fuzzy_g;
    k = find(Vars(:,4)<=0.349&Vars(:,4)>=-0.349);
    M = Vars(k(1,1):3001,1:6);
    Media(j,:)=mean(M);
    Des_Est(j,:)=std(M);
    t = M(:,1);
    for i=2:6
        M1 = M(:,i);
        IND((i-1)+((j-1)*5),1)=trapz(t,M1.^2); %ISE
        IND((i-1)+((j-1)*5),2)=trapz(t,abs(M1)); %IAE
        IND((i-1)+((j-1)*5),3)=trapz(t,t.*(M1.^2)); %ITSE
        IND((i-1)+((j-1)*5),4)=trapz(t,t.*(abs(M1))); %ITAE
        Prom((i-1),1)=IND((i-1)+((j-1)*5),1)+Prom((i-1),1);        
        Prom((i-1),2)=IND((i-1)+((j-1)*5),2)+Prom((i-1),2);
        Prom((i-1),3)=IND((i-1)+((j-1)*5),3)+Prom((i-1),3);        
        Prom((i-1),4)=IND((i-1)+((j-1)*5),4)+Prom((i-1),4);
    end       
end
IND;
Prom_Fuzzy_g=Prom./5;
disp('Fuzzy Adaptativo');
Media
Des_Est;