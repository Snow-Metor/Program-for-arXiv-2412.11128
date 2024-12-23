clc
clear

% parpool(80);



%% MODEL
%%%%  H=sum_j J*(cj*cj^\dagger+h.c.)
%%%%  L1=sqrt(gam1)*c1, L2=sqrt(gam2)*cL^\dagger


%% 参数
L=200;
J=1;
gam1=0.2;
gam2=0.1;
% gamt=0;
t0=0;
dt=2;
tf=5000;

Initial=2;%1为全空态，2为全占据态，3为GHZ态，4为W态

%% Initial state

if Initial==1
    T2_0=0.5*[eye(L),-1i*eye(L);1i*eye(L),eye(L)];
elseif Initial==2
    T2_0=0.5*[eye(L),1i*eye(L);-1i*eye(L),eye(L)];
else
    if Initial==3
        psijv=[0,(2^L-1);1/sqrt(2),1/sqrt(2)];
    else
        psijv=[2.^(0:(L-1));1/sqrt(L)*ones(1,L)];
    end
    com_idi=combvec(1:2*L, 1:2*L)';
    T2_0=zeros(2*L);
    for idi=1:size(com_idi, 1)
        i1=com_idi(idi,2);
        i2=com_idi(idi,1);
        wj1=wpsi(i2,L,psijv);
        wj2=wpsi(i1,L,wj1);

        jie=0;
        len=size(psijv,2);
        for ss=1:len
            xun=find(wj2(1,:)==psijv(1,ss));
            if isempty(xun)
                continue
            else
                jie=jie+psijv(2,ss)*sum(wj2(2,xun));
            end
        end
        T2lin=zeros(2*L);
        T2lin(i1,i2)=jie;
        T2_0=T2_0+T2lin;
    end
end





%% F1, FB

h=diag(J*ones(L-1,1),1)+diag(J*ones(L-1,1),-1);

Ml=zeros(L);
Ml(1,1)=gam1;

Mg=zeros(L);
Mg(end,end)=gam2;

X1=0.5*(Ml+Mg.'-1i*h);
X=[X1,zeros(L);zeros(L),conj(X1)];

Y=[zeros(L),Ml.';Mg.',zeros(L)];

Q=[eye(L),1i*eye(L);1i*eye(L),eye(L)];
Qn=eye(2*L)/Q;
Ix=[zeros(L),eye(L);eye(L),zeros(L)];

F1=-2*Qn*X.'*Q;
FB=2*1i*Qn*Y*Ix*Q;



%% Two Order

t=t0:dt:tf;

s1=L+1;
s2=1;
s3=2*L;
s4=L;

T2ss=lyap(F1,-1i*FB.');
[P,LAM]=eig(F1);
Pn=eye(2*L)/P;
lam=diag(LAM);


% Wr=1i/(2*L)*[zeros(L),-eye(L);...
%     eye(L),zeros(L)];

% Wr=1i/(2*L)*[zeros(L),-diag((-1).^(1:L));...
%     diag((-1).^(1:L)),zeros(L)];


Wr=2*1i/(2*L)*[zeros(L),-diag([ones(1,round(L/2)),zeros(1,round(L/2))]);...
    diag([ones(1,round(L/2)),zeros(1,round(L/2))]),zeros(L)];
exp2W=expm(2*Wr);

expn=zeros(size(t,2),1);
for idt=1:size(t,2)
    ts=t(idt);
    eF1t=P*diag(exp(ts*lam))*Pn;
    T2=eF1t*(T2_0-T2ss)*eF1t.'+T2ss;
    expn(idt)=sqrt(exp(1)*det(T2+(eye(2*L)-T2)*exp2W));
end




%% 绘图

hFig=figure(1);
hold on
set(hFig, 'Position', [50, 50, 1600, 1300]);
plot(t,log(real(expn)),'-','LineWidth',6,'MarkerSize',8)





%% wj*psi
%%%% 计算Wj*Psi所得的新 Psi
function wjpsi=wpsi(j,L,A)
len=size(A,2);
if j-L>0
    for ss=1:len
        bitj=bitget(A(1,ss),j-L);
        if bitj==0
            A(1,ss)=A(1,ss)+2^(j-L-1);
            A(2,ss)=-1i*A(2,ss)/sqrt(2);
        else
            A(1,ss)=A(1,ss)-2^(j-L-1);
            A(2,ss)=1i*A(2,ss)/sqrt(2);
        end
        A(2,ss)=jiou(A(1,ss),j-L)*A(2,ss);
    end
else
    for ss=1:len
        bitj=bitget(A(1,ss),j);
        if bitj==0
            A(1,ss)=A(1,ss)+2^(j-1);
        else
            A(1,ss)=A(1,ss)-2^(j-1);
        end
        A(2,ss)=jiou(A(1,ss),j)*A(2,ss)/sqrt(2);
    end
end
A(1,:)=round(A(1,:));
wjpsi=A;
end

%% 判断奇偶
function jo=jiou(a,j)
a=bitshift(a,-j);
num=0;% 计数二进制中 '1' 的个数
while a>0
    a=bitand(a,a-1);% 消除最右边的 '1'
    num=num+1;
end
jo=(-1)^num;%
end
