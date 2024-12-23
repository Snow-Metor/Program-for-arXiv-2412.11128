clc
clear

%% MODEL
%%%%  H=sum_j J*(cj*cj^\dagger+h.c.)
%%%%  L1=sqrt(gam1)*c1, L2=sqrt(gamL)*cL^\dagger, Lt,j=sqrt(gamtj)*nj
%%%% 增加 Lh1=sqrt(gamh)c1^dagger*cL, Lh2=sqrt(gamh)cL^dagger*c1,
%%%%      Ld1=sqrt(gamd)c1*cL, Ld2=sqrt(gamd)c1^dagger*cL^dagger
%%%%      Lhlj=sqrt(gamlr)*cj^dagger*cj+1,
%%%%      Lhgj=sqrt(gamlr)*cj+1^dagger*cj,

L=4;
J=1;

% jx=1;
% gam1=jx*rand();
% gamL=jx*rand();
% gamt1=jx*rand();
% gamtL=jx*rand();
% det=0.5;
% gam1=0.5*(1+det);
% gamL=0.5*(1-det);
gam1=3.2;
gamL=4.5;

gamt1=0;
gamtL=0;

gamh=0;
gamd=7.7/2;

% gamlr=1;

t0=0;
dt=1000000;
tf=4000000;

Initial=1;%1为全空态，2为全占据态，3为GHZ态，4为W态
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho_1=zeros(2^L);
rho_1(1,1)=1;


rho_2=zeros(2^L);
rho_2(end,end)=1;

rho_3=0.5*(rho_1+rho_2);
rho_3(1,end)=0.5;
rho_3(end,1)=0.5;

psi4=zeros(2^L,1);
psi4(2.^(0:(L-1))+1)=1;
rho_4=1/L*(psi4*psi4');

if Initial==1
    rho=rho_1;
elseif Initial==2
    rho=rho_2;
elseif Initial==3
    rho=rho_3;
elseif Initial==4
    rho=rho_4;
end

C_array = arrayfun(@(x) zeros(2^L), 1:L, 'UniformOutput', false);

for ss1=0:(2^L-1)
    for ss2=1:L
        bit_value=bitand(bitshift(ss1, -(ss2-1)), 1);
        if bit_value==1
            jo=jiou(ss1,ss2);
            C_array{ss2}(ss1-2^(ss2-1)+1,ss1+1)=jo;
        end
    end
end
%%% 注意此时表象为[0;1]对应占据粒子，[1;0]为不占粒子




ev_n=zeros(2^L);
for ss0=1:(L/2)
    ev_n=ev_n+C_array{ss0}'*C_array{ss0};
end

ev_n=2*ev_n/L;

exp_ev_n=kron(expm(ev_n),eye(2^L));





% n1nL=C_array{1}'*C_array{1}*C_array{L}'*C_array{L};
% 
% 
% 
% O1=kron(ev_n,eye(2^L));
% O2=kron(n1nL,eye(2^L));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % w1=1/sqrt(2)*kron((C_array{1}+C_array{1}'),eye(2^L));
% % % wL=1/sqrt(2)*kron((C_array{L}+C_array{L}'),eye(2^L));
% % % wz1=1i/sqrt(2)*kron((C_array{1}-C_array{1}'),eye(2^L));
% % % wzL=1i/sqrt(2)*kron((C_array{L}-C_array{L}'),eye(2^L));
% n1=kron(C_array{1}'*C_array{1},eye(2^L));
% nL=kron(C_array{L}'*C_array{L},eye(2^L));
% d1cL=kron(C_array{1}'*C_array{L},eye(2^L));
% c1dL=kron(C_array{1}*C_array{L}',eye(2^L));
% d1dL=kron(C_array{1}'*C_array{L}',eye(2^L));
% c1cL=kron(C_array{1}*C_array{L},eye(2^L));

j1=1;
j2=1;
j3=L;
j4=L;
o1=kron(C_array{j1}',eye(2^L));
o2=kron(C_array{j2},eye(2^L));
o3=kron(C_array{j3}',eye(2^L));
o4=kron(C_array{j4},eye(2^L));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




H=zeros(2^L);
for ss3=1:(L-1)
H=H+J*(C_array{ss3}'*C_array{ss3+1}+C_array{ss3+1}'*C_array{ss3});
end


L1=sqrt(gam1)*C_array{1};
L2=sqrt(gamL)*C_array{L}';
Lt_1=sqrt(gamt1)*C_array{1}'*C_array{1};
Lt_2=sqrt(gamtL)*C_array{L}'*C_array{L};

Lh_1=sqrt(gamh)*C_array{1}'*C_array{L};
Lh_2=sqrt(gamh)*C_array{L}'*C_array{1};

Ld_1=sqrt(gamd)*C_array{1}*C_array{L};
Ld_2=sqrt(gamd)*C_array{1}'*C_array{L}';


LL=-1i*(kron(H,eye(2^L))-kron(eye(2^L),H.'))...
    +2*kron(L1,conj(L1))+2*kron(L2,conj(L2))...
    +2*kron(Lt_1,conj(Lt_1))+2*kron(Lt_2,conj(Lt_2))...
    ...
    +2*kron(Lh_1,conj(Lh_1))+2*kron(Lh_2,conj(Lh_2))...
    +2*kron(Ld_1,conj(Ld_1))+2*kron(Ld_2,conj(Ld_2))...
    ...
    -kron(L1'*L1,eye(2^L))-kron(eye(2^L),L1.'*conj(L1))...
    -kron(L2'*L2,eye(2^L))-kron(eye(2^L),L2.'*conj(L2))...
    -kron(Lt_1'*Lt_1,eye(2^L))-kron(eye(2^L),Lt_1.'*conj(Lt_1))...
    -kron(Lt_2'*Lt_2,eye(2^L))-kron(eye(2^L),Lt_2.'*conj(Lt_2))...
    ...
    -kron(Lh_1'*Lh_1,eye(2^L))-kron(eye(2^L),Lh_1.'*conj(Lh_1))...
    -kron(Lh_2'*Lh_2,eye(2^L))-kron(eye(2^L),Lh_2.'*conj(Lh_2))...
    -kron(Ld_1'*Ld_1,eye(2^L))-kron(eye(2^L),Ld_1.'*conj(Ld_1))...
    -kron(Ld_2'*Ld_2,eye(2^L))-kron(eye(2^L),Ld_2.'*conj(Ld_2));

% 
% for sss=1:(L-1)
%     L_l=sqrt(gamlr)*C_array{sss}'*C_array{sss+1};
%     L_r=sqrt(gamlr)*C_array{sss+1}'*C_array{sss};
%     LL=LL+2*kron(L_l,conj(L_l))+2*kron(L_r,conj(L_r))...
%         -kron(L_l'*L_l,eye(2^L))-kron(eye(2^L),L_l.'*conj(L_l))...
%         -kron(L_r'*L_r,eye(2^L))-kron(eye(2^L),L_r.'*conj(L_r));
% end



[P,LAM]=eig(LL);
lam=sort(diag(LAM));
Pn=eye(4^L)/P;
t=t0:dt:tf;
elamt=exp(diag(LAM)*t);
lent=length(t);
IL=reshape(eye(2^L),1,[]);
% Eo1=zeros(lent,1);
% Eo2=zeros(lent,1);
% n1nL_wick=zeros(lent,1);
o12=zeros(lent,1);
o34=zeros(lent,1);
o1234=zeros(lent,1);
o1234_wick=zeros(lent,1);
expn=zeros(lent,1);

for ss=1:lent
    rhot=P*diag(elamt(:,ss))*Pn*reshape(rho.',[],1);
    % diag(elamt(:,ss))
    % rhot-reshape(rho.',[],1)
    % Eo1(ss)=IL*O1*rhot;
    % Eo2(ss)=IL*O2*rhot;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % n1nL_wick(ss)=0.25-1i/2*(IL*(w1*wz1+wL*wzL)*rhot)...
    %     +(IL*(w1*wL)*rhot)*(IL*(wz1*wzL)*rhot)...
    %     -(IL*(w1*wz1)*rhot)*(IL*(wL*wzL)*rhot)...
    %     +(IL*(w1*wzL)*rhot)*(IL*(wL*wz1)*rhot);



    % n1nL_wick(ss)=(IL*n1*rhot)*(IL*nL*rhot)...
    %     -(IL*d1dL*rhot)*(IL*c1cL*rhot)...
    %     +(IL*d1cL*rhot)*(IL*c1dL*rhot);

    o12(ss)=IL*o1*o2*rhot;
    o34(ss)=IL*o3*o4*rhot;
    o1234(ss)=IL*o1*o2*o3*o4*rhot;

    o1234_wick(ss)=(IL*o1*o2*rhot)*(IL*o3*o4*rhot)...
        -(IL*o1*o3*rhot)*(IL*o2*o4*rhot)...
        +(IL*o1*o4*rhot)*(IL*o2*o3*rhot);

    expn(ss)=IL*exp_ev_n*rhot;







    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end






% figure(3)
% hold on
% plot(t,real(o1234),'.','markersize',50,'LineWidth', 2.5)
% plot(t,real(o1234_wick),'o','markersize',8,'LineWidth', 3)
% plot(t,real(o1234),'.','markersize',60,'LineWidth', 6)
% plot(t,real(o1234_wick))
figure(2)
hold on
plot(t,real(o1234-o12.*o34),'.','markersize',60,'LineWidth', 6)
% plot(t,real(o1234_wick-o12.*o34))

% % figure(1)
% % hold on
% % plot(t,log(real(expn)),'.','markersize',60,'LineWidth', 2.5)

% % % 
% % % 
% % % figure(100)
% % % hold on
% % % plot(L,real(o1234(end))-real(o1234_wick(end)),'o','LineWidth',4,'MarkerSize',18)


figure(4)
hold on
plot(t,real(o12),'.','markersize',60,'LineWidth', 2.5)
hold on
plot(t,real(o34),'.','markersize',60,'LineWidth', 2.5)






function jo=jiou(a,j)
a=bitshift(a,-j);
num=0;% 计数二进制中 '1' 的个数
while a>0
    a=bitand(a,a-1);% 消除最右边的 '1'
    num=num+1;
end
jo=(-1)^num;%
end

