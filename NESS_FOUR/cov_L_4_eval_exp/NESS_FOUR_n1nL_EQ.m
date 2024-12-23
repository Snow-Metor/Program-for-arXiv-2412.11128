%%% Dynamic_two diss

clc
clear


% global P2 P2n P4 P4n lam2 lam4

%% MODEL
%%%%%%%%   H=sum_j J*[sig_x(j)*sig_x(j+1)+sig_y(j)*sig_y(j+1)]
%%%%%%%%   L_1=sqrt(gam1)*c1, L_L=sqrt(gamL)*cL^\dagger,
%%%%%%%%   Lt1=sqrt(gamt)*c1cL, LtL=sqrt(gamt)*c1^dagger*cL^dagger


parpool(70)

L=4;
J=1;



bianx1=-2:0.005:1.5;
bianx2=-2:0.005:1.5;
len1=size(bianx1,2);
len2=size(bianx2,2);
cov=zeros(len1,len2);
ss1=0;
for sss1=bianx1
gam1=10^sss1;
    ss1=ss1+1;
    ss2=0;
    for sss2=-2:0.005:sss1
    gamL=10^sss2;

        gamt=0.5*(gam1+gamL);

        ss2=ss2+1;

        %% Initial state
        % Initial=2;%1为全空态，2为全占据态，3为GHZ态，4为W态

        save('canshu.mat');
        % Initial_state_T2p_T4p;
        F2p_F4p(2);
        F2p_F4p(4);

        G2p_G4p;


        clear
        load('canshu.mat')

        % load('T2p_T4p.mat')

        load('F2p.mat')
        load('F4p.mat')
        load('G2p_G4p.mat')



        j1=1;
        j2=L;
        j3=L+j1;
        j4=L+j2;

        xv1=[j1,j3];
        xv2=[j2,j4];
        xv3=[j1,j2,j3,j4];



        wei1=findwei(xv1,2*L);
        wei2=findwei(xv2,2*L);
        wei3=findwei(xv3,2*L);

        z2=size(F2p,1);
        z4=size(F4p,1);


        T2bar=-eye(z2)/F2p*G2p;
        T4bar=-eye(z4)/F4p*G4p*T2bar;


        n1nL=0.25-1i/2*(T2bar(wei1)+T2bar(wei2))+T4bar(wei3);
        n1=0.5-1i*T2bar(wei1);
        nL=0.5-1i*T2bar(wei2);

        cov(ss1,ss2)=n1nL-n1*nL;
        save('lin.mat','cov','J','L','bianx1','bianx2')
        % clear
        % load('lin.mat')
    end
end

% jie=[jie,jie((end-1):-1:1)];
% biandet=[biandet_L,-biandet_L((end-1):-1:1)];


% hFig=figure(1);
% set(hFig, 'Position', [50, 50, 1600, 1300]);
% hold on
% plot(biandet,jie,'.-','LineWidth',6,'MarkerSize',50)
% plot(biandet,zeros(1,size(jie,2)),':k','LineWidth',6,'MarkerSize',50)

save('cov.mat','cov','J','L','bianx1','bianx2')





function wei=findwei(xv,n)
len=size(xv,2);
com=nchoosek(1:n, len);
nxv=sum((com-kron(ones(size(com,1),1),xv)).^2,2);
wei=find(nxv==0);
end

