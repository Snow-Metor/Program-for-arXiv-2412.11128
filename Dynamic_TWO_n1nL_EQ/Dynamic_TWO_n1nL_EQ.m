clc
clear
% parpool(80);



%% MODEL
%%%%  H=sum_j J*(cj*cj^\dagger+h.c.)
%%%%  L1=sqrt(gam1)*c1, L2=sqrt(gam2)*cL^\dagger


%% 参数
L=4;
J=1;
gam1=0.2;
gam2=0.1;

t0=0;
dt=0.01;
tf=150;

Initial=1;%1为全空态，2为全占据态，3为GHZ态，4为W态

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
    parfor idi=1:size(com_idi, 1)
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


T4s1=zeros(1,size(t,2));
T2s=zeros(1,size(t,2));
parfor idt=1:size(t,2)
    ts=t(idt);
    eF1t=P*diag(exp(ts*lam))*Pn;
    T2=eF1t*(T2_0-T2ss)*eF1t.'+T2ss;
    T2p=eF1t*T2_0*eF1t.';
    A=T2.'-T2p.';
    B=T2+T2p;
    T4s1(idt)=-0.5*(-A(s2,s1)*B(s3,s4)+A(s3,s1)*B(s2,s4)...
                   -A(s3,s2)*B(s1,s4)-A(s4,s1)*B(s2,s3)...
                   +A(s4,s2)*B(s1,s3)-A(s4,s3)*B(s1,s2));
    T2s(idt)=T2(L+1,1)+T2(2*L,L);
end

%% Four Order
T4s2=zeros(1,size(t,2));
com_idx=combvec(1:2*L, 1:2*L)';% 得到所有 (j1, j2) 的组合
parfor idx=1:size(com_idx, 1)
    j1=com_idx(idx,2);
    j2=com_idx(idx,1);

    temp_sum=zeros(1,size(t,2));
    for idy=1:size(com_idx, 1)
        j3=com_idx(idy,2);
        j4=com_idx(idy,1);

        

        if Initial<=2
            T4_0=T2_0(j1,j2)*T2_0(j3,j4)...
                -T2_0(j1,j3)*T2_0(j2,j4)...
                +T2_0(j1,j4)*T2_0(j2,j3);
        else
            T4_0=T40(j1,j2,j3,j4,L,psijv);
        end
       




        temp_sum=temp_sum+...
              (P(s1,:).*Pn(:,j1).'*exp(lam*t))...
            .*(P(s2,:).*Pn(:,j2).'*exp(lam*t))...
            .*(P(s3,:).*Pn(:,j3).'*exp(lam*t))...
            .*(P(s4,:).*Pn(:,j4).'*exp(lam*t))*T4_0;
    end
    T4s2=T4s2+temp_sum;
end


T4s=T4s1+T4s2;%得到Ts1,s2,s3,s4的值


n1nL=0.25+1i/2*T2s-T4s;


% name=['L_',num2str(L),' J_',num2str(J),' gam1_',num2str(gam1),...
%     ' gam2_',num2str(gam2),' gamt_',num2str(gamt),' Initial_',num2str(Initial),'.mat'];
% save(name,'t','n1nL')
%% 绘图

hFig=figure(1);
set(hFig, 'Position', [50, 50, 1600, 1300]);
hold on
plot(t,real(n1nL),'-','LineWidth',6,'MarkerSize',8)


%% 计算Tj1j2j3j4_0
function T4_0=T40(j1,j2,j3,j4,L,psijv)


wj1=wpsi(j4,L,psijv);
wj2=wpsi(j3,L,wj1);
wj3=wpsi(j2,L,wj2);
wj4=wpsi(j1,L,wj3);
jie=0;
len=size(psijv,2);
for ss=1:len
    xun=find(wj4(1,:)==psijv(1,ss));
    if isempty(xun)
        continue
    else
        jie=jie+psijv(2,ss)*sum(wj4(2,xun));
    end
end
T4_0=jie;

end


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
jo=(-1)^num;
end
