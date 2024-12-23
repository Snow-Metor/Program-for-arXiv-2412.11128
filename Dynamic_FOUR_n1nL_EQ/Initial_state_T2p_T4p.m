%%%  T2p_0, T4p_0
clc
clear

load('canshu.mat')


%% T2p_0

com2=nchoosek(1:(2*L),2);
len2=size(com2,1);
T2p_0=zeros(len2,1);
if Initial==1
    T2_0=0.5*[eye(L),-1i*eye(L);1i*eye(L),eye(L)];
    for ss1=1:len2
        T2p_0(ss1)=T2_0(com2(ss1,1),com2(ss1,2));
    end
elseif Initial==2
    T2_0=0.5*[eye(L),1i*eye(L);-1i*eye(L),eye(L)];
    for ss1=1:len2
        T2p_0(ss1)=T2_0(com2(ss1,1),com2(ss1,2));
    end
else
    if Initial==3
        psijv=[0,(2^L-1);1/sqrt(2),1/sqrt(2)];
    else
        psijv=[2.^(0:(L-1));1/sqrt(L)*ones(1,L)];
    end
    
    T2_0=zeros(2*L);
    for idi=1:len2
        i1=com2(idi,1);
        i2=com2(idi,2);
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
        T2p_0(idi)=jie;
    end
end

%% T4p_0

com4=nchoosek(1:(2*L),4);
len4=size(com4,1);
T4p_0=zeros(len4,1);
if Initial<=2
    for ss3=1:len4
        j1=com4(ss3,1);
        j2=com4(ss3,2);
        j3=com4(ss3,3);
        j4=com4(ss3,4);
        T4p_0(ss3)=T2_0(j1,j2)*T2_0(j3,j4)...
            -T2_0(j1,j3)*T2_0(j2,j4)...
            +T2_0(j1,j4)*T2_0(j2,j3);
    end
else
    for ss3=1:len4
        j1=com4(ss3,1);
        j2=com4(ss3,2);
        j3=com4(ss3,3);
        j4=com4(ss3,4);
        T4p_0(ss3)=T40(j1,j2,j3,j4,L,psijv);
    end
end

save('T2p_T4p.mat','T2p_0','T4p_0')

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
jo=(-1)^num;%
end
