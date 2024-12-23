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
dt=0.05;
tf=120;


Initial=3;%1为全空态，2为全占据态，3为GHZ态，4为W态

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





t=t0:dt:tf;

s1=L+1;
s2=1;
s3=2*L;
s4=L;

T2ss=lyap(F1,-1i*FB.');
[P,LAM]=eig(F1);
Pn=eye(2*L)/P;
lam=diag(LAM);


n1nL_wick=zeros(size(t,2),1);



%%%%%%%%%%%%%%%%%%%%%
T2_01=0.5*[eye(L),-1i*eye(L);1i*eye(L),eye(L)];
T2_02=0.5*[eye(L),1i*eye(L);-1i*eye(L),eye(L)];
n1nL_wick_xiu=zeros(size(t,2),1);
%%%%%%%%%%%%%%%%%%%



for idt=1:size(t,2)
    ts=t(idt);
    eF1t=P*diag(exp(ts*lam))*Pn;
    T2=eF1t*(T2_0-T2ss)*eF1t.'+T2ss;
    T4_wick=T2(L+1,1)*T2(2*L,L)...
        -T2(L+1,2*L)*T2(1,L)...
        +T2(L+1,L)*T2(1,2*L);
    n1nL_wick(idt)=0.25+0.5*1i*(T2(L+1,1)+T2(2*L,L))-T4_wick;


    %%%%%%%%%修改wick对于GHZ态
        T2_1=eF1t*(T2_01-T2ss)*eF1t.'+T2ss;
        T2_2=eF1t*(T2_02-T2ss)*eF1t.'+T2ss;
        T4_wick_1=T2_1(L+1,1)*T2_1(2*L,L)...
           -T2_1(L+1,2*L)*T2_1(1,L)...
           +T2_1(L+1,L)*T2_1(1,2*L);
        T4_wick_2=T2_2(L+1,1)*T2_2(2*L,L)...
           -T2_2(L+1,2*L)*T2_2(1,L)...
           +T2_2(L+1,L)*T2_2(1,2*L);
        n1nL_wick_xiu(idt)=0.25...
            +0.25*1i*(T2_1(L+1,1)+T2_1(2*L,L)+T2_2(L+1,1)+T2_2(2*L,L))...
            -0.5*(T4_wick_1+T4_wick_2);
    %%%%%%%%%%%%%%%%%%%%%%

    % % 
    % % %%%%%%%%%%修改wick对于W态
    % % lin=0;
    % % for s1=1:L
    % %     T200=0.5*[eye(L),1i*eye(L);-1i*eye(L),eye(L)];
    % %     T200(s1,L+s1)=-T200(s1,L+s1);
    % %     T200(L+s1,s1)=-T200(L+s1,s1);
    % %     T2_00=eF1t*(T200-T2ss)*eF1t.'+T2ss;
    % %     T4_wick_0=T2_00(L+1,1)*T2_00(2*L,L)...
    % %         -T2_00(L+1,2*L)*T2_00(1,L)...
    % %         +T2_00(L+1,L)*T2_00(1,2*L);
    % %     lin=lin+0.25+0.5*1i*(T2_00(L+1,1)+T2_00(2*L,L))-T4_wick_0;
    % % end
    % % n1nL_wick_xiu(idt)=1/L*lin;
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%






end




%% 绘图

hFig=figure(1);
set(hFig, 'Position', [50, 50, 1600, 1300]);
hold on
plot(t(1:end),real(n1nL_wick_xiu(1:end)),...
    '-','markersize',8,'LineWidth', 6)

plot(t(1:end),real(n1nL_wick_xiu(1:end)),...
    'o','markersize',20,'LineWidth', 4)

% 获取当前坐标轴的默认颜色顺序
defaultColors = get(gca, 'ColorOrder');
% 获取当前绘制的对象数量
nCurves = numel(get(gca, 'Children')); % 获取当前绘图对象数量
% 计算下一个颜色的索引
cor=mod(nCurves, size(defaultColors,1))+1;
% 淡化颜色的方法（与白色混合）
fadeFactor=0.25; % 淡化比例
fadedColor1=(1-fadeFactor)*defaultColors(cor,:)+fadeFactor * [1 1 1]; % 淡化第cor个颜色

% plot(t,real(n1nL_wick),'o','markersize',18,'LineWidth', 3)
plot(t,real(n1nL_wick),':',...
    'markersize',20,'LineWidth', 4,'Color',fadedColor1)






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
