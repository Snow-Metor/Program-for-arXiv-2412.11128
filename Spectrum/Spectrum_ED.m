clc
clear

%% MODEL
%%%%  H=sum_j J*(cj*cj^\dagger+h.c.)
%%%%  L1=sqrt(gam1)*c1, L2=sqrt(gamL)*cL^\dagger,
%%%%  Lt1=sqrt(gamt)*c1cL, L2=sqrt(gamt)*c1^dagger*cL^\dagger,

L=4;
J=1;

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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





H=zeros(2^L);
for ss3=1:(L-1)
H=H+J*(C_array{ss3}'*C_array{ss3+1}+C_array{ss3+1}'*C_array{ss3});
end


L1=sqrt(gam1)*C_array{1};
L2=sqrt(gamL)*C_array{L}';
Lt_1=sqrt(gamt)*C_array{1}*C_array{L};
Lt_2=sqrt(gamt)*C_array{1}'*C_array{L}';


LL=-1i*(kron(H,eye(2^L))-kron(eye(2^L),H.'))...
    +2*kron(L1,conj(L1))+2*kron(L2,conj(L2))...
    +2*kron(Lt_1,conj(Lt_1))+2*kron(Lt_2,conj(Lt_2))...
    -kron(L1'*L1,eye(2^L))-kron(eye(2^L),L1.'*conj(L1))...
    -kron(L2'*L2,eye(2^L))-kron(eye(2^L),L2.'*conj(L2))...
    -kron(Lt_1'*Lt_1,eye(2^L))-kron(eye(2^L),Lt_1.'*conj(Lt_1))...
    -kron(Lt_2'*Lt_2,eye(2^L))-kron(eye(2^L),Lt_2.'*conj(Lt_2));


LAM=eig(LL);


hFig=figure(1);
set(hFig, 'Position', [50, 50, 1600, 1300]);
hold on
plot(real(LAM),imag(LAM),'.','markersize',50,'LineWidth', 2.5)












function jo=jiou(a,j)
a=bitshift(a,-j);
num=0;% 计数二进制中 '1' 的个数
while a>0
    a=bitand(a,a-1);% 消除最右边的 '1'
    num=num+1;
end
jo=(-1)^num;%
end

