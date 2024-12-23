%%% Dynamic_two diss

clc
clear


% global P2 P2n P4 P4n lam2 lam4

%% MODEL
%%%%%%%%   H=sum_j J*[sig_x(j)*sig_x(j+1)+sig_y(j)*sig_y(j+1)]
%%%%%%%%   L_1=sqrt(gam1)*c1, L_L=sqrt(gamL)*cL^\dagger
%%%%%%%%   Lt1=sqrt(gamt)*c1cL, Lt2=sqrt(gamt)*c1^dagger*cL^dagger

L=4;
t=1;
det=0;
gam1=0.5*(1+det);
gamL=0.5*(1-det);
gamt=0.5;


t0=0;
dt=0.1;
tf=35;

%% Initial state
Initial=1;%1为全空态，2为全占据态，3为GHZ态，4为W态
% 
save('canshu.mat');
Initial_state_T2p_T4p;
F2p_F4p(2);
F2p_F4p(4);

G2p_G4p;


clear
load('canshu.mat')

load('T2p_T4p.mat')

load('F2p.mat')
load('F4p.mat')
load('G2p_G4p.mat')



[P2,LAM2]=eig(F2p);
P2n=eye(size(P2,1))/P2;
lam2=diag(LAM2);
F2pn=eye(size(P2,1))/F2p;

[P4,LAM4]=eig(F4p);
P4n=eye(size(P4,1))/P4;
lam4=diag(LAM4);
F4pn=eye(size(P4,1))/F4p;


X=lyap(-F4p,F2p,-G4p);

j1=1;
j2=L;
j3=L+j1;
j4=L+j2;

xv1=[j1,j3];
xv2=[j2,j4];
xv3=[j1,j2,j3,j4];
% xv1=[1,L+1];
% xv2=[L,2*L];
% xv3=[1,L,L+1,2*L];

wei1=findwei(xv1,2*L);
wei2=findwei(xv2,2*L);
wei3=findwei(xv3,2*L);


w12=findwei([j1,j2],2*L);
w34=findwei([j3,j4],2*L);
w13=findwei([j1,j3],2*L);
w24=findwei([j2,j4],2*L);
w14=findwei([j1,j4],2*L);
w23=findwei([j2,j3],2*L);
n1nL=[];
n1nL_wick=[];
n1nL_n1_nL=[];
n1nL_n1_nL_WICK=[];
for tt=t0:dt:tf
    % tt
    jie_n1nL=0.25...
        -1i/2*(T2p(wei1,tt,P2,P2n,lam2,T2p_0,F2pn,G2p)+T2p(wei2,tt,P2,P2n,lam2,T2p_0,F2pn,G2p))...
        +T4p(wei3,tt,P2,P2n,P4,P4n,lam2,lam4,F2pn,F4pn,G2p,G4p,X,T2p_0,T4p_0);
    n1nL=[n1nL,jie_n1nL];

    jie_wick=0.25...
        -1i/2*(T2p(wei1,tt,P2,P2n,lam2,T2p_0,F2pn,G2p)+T2p(wei2,tt,P2,P2n,lam2,T2p_0,F2pn,G2p))...
        +T2p(w12,tt,P2,P2n,lam2,T2p_0,F2pn,G2p)*T2p(w34,tt,P2,P2n,lam2,T2p_0,F2pn,G2p)...
        -T2p(w13,tt,P2,P2n,lam2,T2p_0,F2pn,G2p)*T2p(w24,tt,P2,P2n,lam2,T2p_0,F2pn,G2p)...
        +T2p(w14,tt,P2,P2n,lam2,T2p_0,F2pn,G2p)*T2p(w23,tt,P2,P2n,lam2,T2p_0,F2pn,G2p);
    n1nL_wick=[n1nL_wick,jie_wick];


    n1=0.5-1i*T2p(wei1,tt,P2,P2n,lam2,T2p_0,F2pn,G2p);
    nL=0.5-1i*T2p(wei2,tt,P2,P2n,lam2,T2p_0,F2pn,G2p);
    n1nL_n1_nL=[n1nL_n1_nL,jie_n1nL-n1*nL];%% <n1nL>-<n1><nL>
    n1nL_n1_nL_WICK=[n1nL_n1_nL_WICK,jie_wick-n1*nL];%% <n1nL>_wick-<n1><nL>




end


%% 绘图

hFig=figure(1);
set(hFig, 'Position', [50, 50, 1600, 1300]);

plot(t0:dt:tf,real(n1nL),'-','LineWidth',6,'MarkerSize',8)
hold on

% 获取当前坐标轴的默认颜色顺序
defaultColors = get(gca, 'ColorOrder');
% 获取当前绘制的对象数量
nCurves = numel(get(gca, 'Children')); % 获取当前绘图对象数量
% 计算下一个颜色的索引
cor=mod(nCurves, size(defaultColors,1))+1;
% 淡化颜色的方法（与白色混合）
fadeFactor=0.25; % 淡化比例
fadedColor1=(1-fadeFactor)*defaultColors(cor,:)+fadeFactor * [1 1 1]; % 淡化第cor个颜色

plot(t0:dt:tf,real(n1nL_wick),':','markersize',8,'LineWidth', 6,'Color',fadedColor1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hFig=figure(2);
hold on
set(hFig, 'Position', [50, 50, 1600, 1300]);
plot(t0:dt:tf,real(n1nL_n1_nL),'-','LineWidth',6,'MarkerSize',8)
% 获取当前坐标轴的默认颜色顺序
defaultColors = get(gca, 'ColorOrder');
% 获取当前绘制的对象数量
nCurves = numel(get(gca, 'Children')); % 获取当前绘图对象数量
% 计算下一个颜色的索引
cor=mod(nCurves, size(defaultColors,1))+1;
% 淡化颜色的方法（与白色混合）
fadeFactor=0.25; % 淡化比例
fadedColor1=(1-fadeFactor)*defaultColors(cor,:)+fadeFactor * [1 1 1]; % 淡化第cor个颜色

plot(t0:dt:tf,real(n1nL_n1_nL_WICK),':','LineWidth',6,'MarkerSize',8)

save('n1nL_and_wick.mat','n1nL','n1nL_wick','n1nL_n1_nL','n1nL_n1_nL_WICK','gam1','gamL','t0','dt','tf')









function wei=findwei(xv,n)
len=size(xv,2);
com=nchoosek(1:n, len);
nxv=sum((com-kron(ones(size(com,1),1),xv)).^2,2);
wei=find(nxv==0);
end

function jie=T2p(wei,t,P2,P2n,lam2,T2p_0,F2pn,G2p)
jie=P2(wei,:)*diag(exp(lam2*t))*P2n*(T2p_0+F2pn*G2p)...
    -F2pn(wei,:)*G2p;
end


function jie=T4p(wei,t,P2,P2n,P4,P4n,lam2,lam4,F2pn,F4pn,G2p,G4p,X,T2p_0,T4p_0)
jie=P4(wei,:)*diag(exp(lam4*t))*P4n*(T4p_0-X*(T2p_0+F2pn*G2p)-F4pn*G4p*F2pn*G2p)...
    +X(wei,:)*P2*diag(exp(lam2*t))*P2n*(T2p_0+F2pn*G2p)...
    +F4pn(wei,:)*G4p*F2pn*G2p;
end



