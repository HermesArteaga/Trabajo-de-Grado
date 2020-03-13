
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Linealizacion Mediante Jacobiano
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculo de los ganacias K del LQR

clear all
clc
syms Q1 dQ1 Q2 dQ2 Vm;
J = 0.0023146;        %inercia pilar central
%J = 0.002765;        
mb = 0.245;           %masa del brazo
mp = 0.015;           %masa del pundulo
lb = 0.15;            %longitud pivote del brazo al pivote del pendulo
lp = 0.2;             %logitud total pendulo
g = 9.8;              %gravedad
Ke = 0.353;
Kt = 0.15;
R = 2.34;
A = J+((((1/3)*mb)+mp)*(lb^2));
B = (1/3)*mp*(lp^2);
C = (1/2)*mp*lb*lp;
D = (1/2)*mp*g*lp;
to = ((Ke/R)*Vm)-(((Ke*Kt)/R)*dQ1);     %Voltaje a Torque
den = (A*B)-(C^2)+(((B^2)+(C^2))*((sin(Q2))^2));
%% Espacio de Estados

X1 = dQ1;
X2 = ((B*C*(((sin(Q2))^2)-1)*sin(Q2)*(dQ1^2))...
        -(2*(B^2)*cos(Q2)*sin(Q2)*dQ1*dQ2)...
        +(B*C*sin(Q2)*(dQ2^2))...
        -(C*D*cos(Q2)*sin(Q2))...
        +(B*to))/den;
X3 = dQ2;
X4 = ((B*(A+(B*((sin(Q2))^2)))*cos(Q2)*sin(Q2)*(dQ1^2))...
       +(2*B*C*(1-(sin(Q2))^2)*sin(Q2)*dQ1*dQ2)...
       -((C^2)*cos(Q2)*sin(Q2)*dQ2)...
       +(D*(A+(B*(sin(Q2))^2))*sin(Q2))...
       -(C*cos(Q2)*to))/den;
%% Linealizacion 

F = [X1, X2, X3, X4];
X = [Q1 dQ1 Q2 dQ2];
%Jacobiano
Aa = jacobian(F,X);       
P = [0 0 0 0];
Af = subs(Aa,X,P);
Af = subs(Af,Vm,0);
Af=simplify(Af);
Af = double(Af);
Ba = jacobian(F,Vm);
Bf = subs(Ba,Vm,0);
Bf = subs(Bf,X,P);
Bf = double(Bf);
%% Calculo de las ganancias K del LQR
Qf=[1 0 0 0; 0 1.5 0 0; 0 0 4 0; 0 0 0 12];
Rf = 16;
%LQR Continuo
% Kcn = lqr(Af,Bf,Qf,Rf);         
% K1=Kcn(1)
% K2=Kcn(2)
% K3=Kcn(3)
% K4=Kcn(4)

%LQR Discreto
Kds = lqrd(Af,Bf,Qf,Rf,0.01)     
K1=Kds(1)
K2=Kds(2)
K3=Kds(3)
K4=Kds(4)




