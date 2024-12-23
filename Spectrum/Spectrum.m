clc
clear


% global P2 P2n P4 P4n lam2 lam4

%% MODEL
%%%%%%%%   H=sum_j J*[sig_x(j)*sig_x(j+1)+sig_y(j)*sig_y(j+1)]
%%%%%%%%   L_1=sqrt(gam1)*c1, L_L=sqrt(gamL)*cL^\dagger,
%%%%%%%%   L_1=sqrt(gamt)*c1*cL, L_L=sqrt(gamt)*c1^dagger*cL^\dagger,

L=4;
t=1;
jx=1;
% gam1=jx*rand();
% gamL=jx*rand();
% gamt1=jx*rand();
% gamtL=jx*rand();
s=0;
gam=0.5;
gam1=gam*(1+s);
gamL=gam*(1-s);
gamt=0.5;



%% Initial state
Initial=1;%1为全空态，2为全占据态，3为GHZ态，4为W态

save('canshu.mat');

for ss=1:2:(2*L)
    F2p_F4p(ss);
    clear Fnp
    load('Fnp.mat')
    lam=eig(Fnp);
    hFig=figure(1);
    set(hFig, 'Position', [50, 50, 1600, 1300]);
    hold on
    plot(real(lam),imag(lam),'o','markersize',22,'LineWidth',5)
end


