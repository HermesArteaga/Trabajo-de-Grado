% tetha=angulo recorrido por el brazo positivo en sentido antihoario
% alfa=angulo recorrido por el pendulo positivo en sentido antihoario
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Linealizacion Mediante Jacobiano con ciclo FOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Mediante un ciclo FOR cambiamos los puntos donde deseamos 
%linealizar para poder generar las ganancias K en cada uno de los 
%puntos y asi obtener un controlador LQR con Gain Scheduling
%Creamos una matriz Cint la cual mediante un integrador nos ayudara  
%sacar al pendulo de la zona muerta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
syms Q1 dQ1 Q2 dQ2 Vm;
J = 0.0023146;        %inercia pilbr central
%J = 0.002765;
ma = 0.245;         %masa del brazo
mp = 0.015;         %masa del pundulo  
la = 0.15;          %longitud pivote del brazo al pivote del pendulo
lp = 0.2;           %logitud total pendulo
g = 9.8;            %gravedad
Ke = 0.353;
Kt = 0.15;
R = 2.34;
A = J+((((1/3)*ma)+mp)*(la^2));
B = (1/3)*mp*(lp^2);
C = (1/2)*mp*la*lp;
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
 %Puntos de Linealizacion
    %Creacion de la Matriz que contiene los puntos
    Rango = -20*(0.0175):0.0175:20*(0.0175);
    Rang = -20*(0.0175):0.0175:0;   %Matriz que contiene los puntos.
    t=size(Rang);
    x=t(1,2)*2;
    b=1;
        for i=Rang
        P = [0 0 i 0];
        Vmf = (R*((D*(A+(B*(sin(i))^2))*sin(i))+(C*D*cos(i)*sin(i))))...
              /(Ke*((C*cos(i))+B));     %Calculo del voltaje en cada punto
        Af = subs(Aa,X,P);
        Af = subs(Af,Vm,Vmf);
        Af=simplify(Af);
        Af = double(Af);
        Ba = jacobian(F,Vm);
        Bf = subs(Ba,Vm,Vmf);
        Bf = subs(Bf,X,P);
        Bf = double(Bf);        
        %% Calculo de las ganancias K del LQR
        %Matriz aumentada
        Cint = [1 0 0 0];
        Z = zeros(4,1);
        AA = [Af Z;-Cint 0];
        BA = [Bf; 0];
        Q = [0.05 0 0 0 0; 0 1 0 0 0; 0 0 4 0 0; 0 0 0 15 0; 0 0 0 0 0.2];
        R = 20;        
        %Kds = lqr(AA,BA,QA,RA);        %LQR Continuo 
        Kds = lqrd(AA,BA,Q,R,0.01);   %LQR Discreto
        
        %Creacion Matriz que contiene las ganacias K para diferetes puntos
        MKi(b,:)=Kds(1, :);
        MKi((x)-b,:)=Kds(1, :);
        b=b+1;
        
        end
        %Columnas con los diferentes valores de K
        M1=MKi(:,1);
        M2=MKi(:,2);
        M3=MKi(:,3);
        M4=MKi(:,4);
        M5=MKi(:,5);