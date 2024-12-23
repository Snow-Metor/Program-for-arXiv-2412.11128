

clc
clear


%% MODEL
%%%%%%%%   H=sum_j J*[sig_x(j)*sig_x(j+1)+sig_y(j)*sig_y(j+1)]
%%%%%%%%   L_1=sqrt(gam1)*c1, L_L=sqrt(gamL)*cL^\dagger,
%%%%%%%%   Lt1=sqrt(gamt)*c1cL, LtL=sqrt(gamt)*c1^dagger*cL^dagger
parpool(2)
L=4;

syms x y z J 


gaml=x-y;
gamg=x+y;
gamt=z;


F1=-diag([gaml,zeros(1,L-2),gamg,gaml,zeros(1,L-2),gamg]);
X0=diag(J*ones(1,L-1),1)+diag(J*ones(1,L-1),-1);
% X0=diag([J,Je,J],1)+diag([J,Je,J],-1);
F1=F1+[zeros(L),-X0;X0,zeros(L)];

FB=1i*diag([gaml,zeros(1,L-2),gamg,gaml,zeros(1,L-2),gamg])...
    +diag([-gaml,zeros(1,L-2),gamg],L)+diag([gaml,zeros(1,L-2),-gamg],-L);



%% U1
U1_real=sym(zeros(2*L));
U1_real(1,L)=1/2;
U1_real(L+1,2*L)=-1/2;
% U1_real(1,1)=1/2;
% U1_real(L+1,L+1)=1/2;

U1_real=U1_real*sqrt(gamt);


U1_imag=sym(zeros(2*L));
U1_imag(1,2*L)=-1/2;
U1_imag(L+1,L)=-1/2;
% U1_imag(1,L+1)=-1/2;
% U1_imag(L+1,1)=1/2;

U1_imag=U1_imag*sqrt(gamt);

U1a_real=U1_real-U1_real.';
U1a_imag=U1_imag-U1_imag.';
U1=U1_real+1i*U1_imag;
U1a=U1-U1.';

%% U2

U2_real=sym(zeros(2*L));
U2_real(1,L)=1/2;
U2_real(L+1,2*L)=-1/2;
% U2_real(L,L)=1/2;
% U2_real(2*L,2*L)=1/2;

U2_real=U2_real*sqrt(gamt);

U2_imag=sym(zeros(2*L));
U2_imag(1,2*L)=1/2;
U2_imag(L+1,L)=1/2;
% U2_imag(L,2*L)=-1/2;
% U2_imag(2*L,L)=1/2;

U2_imag=U2_imag*sqrt(gamt);


U2a_real=U2_real-U2_real.';
U2a_imag=U2_imag-U2_imag.';
U2=U2_real+1i*U2_imag;
U2a=U2-U2.';


U1ad=U1a_real.'-1i*U1a_imag.';
U2ad=U2a_real.'-1i*U2a_imag.';







save('canshu.mat');

%%

F2p_F4p(2);
F2p_F4p(4);
G2p_G4p;


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



n1=0.5-1i*T2bar(wei1);
nL=0.5-1i*T2bar(wei2);


ewei3=zeros(1,z4);
ewei3(wei3)=1;


save('G2p.mat','G2p')
save('G4p.mat','G4p')
save('T2bar.mat','T2bar')
save('ewei3.mat','ewei3')

save('lin.mat')







%%
[rows2, cols2, values2] = find(F2p);
fileID = fopen('F2p.txt', 'w');
fprintf(fileID, '{\n');
for i = 1:length(rows2)
    fprintf(fileID, ' {%d, %d, %s}', rows2(i), cols2(i), char(values2(i)));
    if i < length(rows2)
        fprintf(fileID, ',\n');
    else
        fprintf(fileID, '\n'); 
    end
end
fprintf(fileID, '}\n');
fclose(fileID);






%%
[rows2, cols2, values2] = find(G2p);
fileID = fopen('G2p.txt', 'w');
fprintf(fileID, '{\n');
for i = 1:length(rows2)
    fprintf(fileID, ' {%d, %d, %s}', rows2(i), cols2(i), char(values2(i)));
    if i < length(rows2)
        fprintf(fileID, ',\n');
    else
        fprintf(fileID, '\n'); 
    end
end
fprintf(fileID, '}\n');
fclose(fileID);


















%%
% 获取非零元素的行、列和值
[rows1, cols1, values1] = find(ewei3);

% 打开文件以写入
fileID = fopen('ewei3.txt', 'w');

% 写入开头的左大括号
fprintf(fileID, '{\n');

% 将行列值和对应的符号值写入文件
for i = 1:length(rows1)
    % 将符号值转为字符串
    fprintf(fileID, ' {%d, %d, %s}', rows1(i), cols1(i), char(values1(i)));
    if i < length(rows1)
        fprintf(fileID, ',\n');  % 如果不是最后一个元素，添加逗号
    else
        fprintf(fileID, '\n');  % 如果是最后一个元素，换行
    end
end

% 写入文件尾部的右大括号
fprintf(fileID, '}\n');

% 关闭文件
fclose(fileID);



%%
[rows2, cols2, values2] = find(G4p);
fileID = fopen('G4p.txt', 'w');
fprintf(fileID, '{\n');
for i = 1:length(rows2)
    fprintf(fileID, ' {%d, %d, %s}', rows2(i), cols2(i), char(values2(i)));
    if i < length(rows2)
        fprintf(fileID, ',\n');
    else
        fprintf(fileID, '\n'); 
    end
end
fprintf(fileID, '}\n');
fclose(fileID);


%%
[rows2, cols2, values2] = find(F4p);
fileID = fopen('F4p.txt', 'w');
fprintf(fileID, '{\n');
for i = 1:length(rows2)
    fprintf(fileID, ' {%d, %d, %s}', rows2(i), cols2(i), char(values2(i)));
    if i < length(rows2)
        fprintf(fileID, ',\n');
    else
        fprintf(fileID, '\n'); 
    end
end
fprintf(fileID, '}\n');
fclose(fileID);


%%
[rows2, cols2, values2] = find(T2bar);
fileID = fopen('T2bar.txt', 'w');
fprintf(fileID, '{\n');
for i = 1:length(rows2)
    fprintf(fileID, ' {%d, %d, %s}', rows2(i), cols2(i), char(values2(i)));
    if i < length(rows2)
        fprintf(fileID, ',\n');
    else
        fprintf(fileID, '\n'); 
    end
end
fprintf(fileID, '}\n');
fclose(fileID);




%%
% z4=70;
% 
niL=sym(zeros(1,2));
% niL=cell(1,z4);
K=21;

parfor jjj=1:2
    F4plin=F4p;
    F4plin(jjj+K-1,:)=[];
    F4plin(:,wei3)=[];
    niL(jjj)=simplify((-1)^(jjj+K-1+wei3)*det(F4plin));


end
save('niL.mat','niL')
% % % T4wei3=-niL*G4p*T2bar/det(F4p);
% % % 
% % % % T4wei3=-ewei3*eye(z4)/F4p*G4p*T2bar;
% % % % simplify(T4wei30-T4wei3)
% % % 
% % % n1nL=simplify(0.25-1i/2*(T2bar(wei1)+T2bar(wei2))+T4wei3);
% % % 
% % % cov=simplify(n1nL-n1*nL);
% % % 
% % % 
% % % 
% % % save('cov.mat','n1nL','n1','nL','cov')
% % % 
% % % save('all.mat')
% % % 


function wei=findwei(xv,n)
len=size(xv,2);
com=nchoosek(1:n, len);
nxv=sum((com-kron(ones(size(com,1),1),xv)).^2,2);
wei=find(nxv==0);
end

