clear all
clc
syms Q1 dQ1 Q2 dQ2 Vm;
J=0.0023146;
%J=0.0027;
ma=0.245;
mp=0.015;
la=0.15;
lp=0.2;
g=9.8;
Ke=0.353;
Kt=0.15;
R=2.34;
A = J+((((1/3)*ma)+mp)*(la^2));
B = (1/3)*mp*(lp^2);
C = (1/2)*mp*la*lp;
D = (1/2)*mp*g*lp;
to = ((Ke/R)*Vm)-(((Ke*Kt)/R)*dQ1);
den = (A*B)-(C^2)+(((B^2)+(C^2))*((sin(Q2))^2));

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
F = [X1, X2, X3, X4];
X = [Q1 dQ1 Q2 dQ2];
Aa = jacobian(F,X);
P = [0 0 0  0];
Af = subs(Aa,X,P);
Af = subs(Af,Vm,0);
Af=simplify(Af);
%%
Ba = jacobian(F,Vm);
Bf = subs(Ba,Vm,0);
Bf = subs(Bf,X,P);

Qf=[1 0 0 0; 0 1.5 0 0; 0 0 4 0; 0 0 0 12];
Rf =16;

Af = double(Af);
Bf = double(Bf);

Kcn = lqr(Af,Bf,Qf,Rf);
Kds = lqrd(Af,Bf,Qf,Rf,0.01)
K1=Kds(1);
K2=Kds(2);
K3=Kds(3);
K4=Kds(4);

%%
% Cint = [1 0 0 0];
% Z = zeros(4,1);
% AA = [Af Z;-Cint 0];
% BA = [Bf; 0];
% QA = [0.02 0 0 0 0; 0 1.5 0 0 0; 0 0 2 0 0; 0 0 0 12 0; 0 0 0 0 0.2];
% RA = 19;
% 
% Kds = lqrd(AA,BA,QA,RA,0.01)
% K1 = Kds(1);
% K2 = Kds(2);
% K3 = Kds(3);
% K4 = Kds(4);
% K5 = Kds(5);


