

clc
clear
L=4;
syms x y z J 
load('G4p.mat')
load('T2bar.mat')
load('DET.mat')


Rvec=G4p*T2bar;
jie=0;
for ss=1:35
    load(['niL',num2str(2*ss-1),'.mat'])
    linR=Rvec((2*ss-1):(2*ss));
    jie=jie+niL*linR;
end
T4wei3=-jie/DET;




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




n1=0.5-1i*T2bar(wei1)
nL=0.5-1i*T2bar(wei2)

% 将符号变量转换为字符串
n1str = char(n1);

% 写入到 txt 文件
fileID = fopen('n1.txt', 'w'); % 打开文件
fprintf(fileID, '%s\n', n1str);             % 写入符号表达式
fclose(fileID);  

% 将符号变量转换为字符串
nLstr = char(nL);

% 写入到 txt 文件
fileID = fopen('nL.txt', 'w'); % 打开文件
fprintf(fileID, '%s\n', nLstr);             % 写入符号表达式
fclose(fileID);   

n1nL=simplify(0.25-1i/2*(T2bar(wei1)+T2bar(wei2))+T4wei3)

cov=simplify(n1nL-n1*nL)


% 将符号变量转换为字符串
n1nLstr = char(n1nL);

% 写入到 txt 文件
fileID = fopen('n1nL.txt', 'w'); % 打开文件
fprintf(fileID, '%s\n', n1nLstr);             % 写入符号表达式
fclose(fileID);                               % 关闭文件

% 将符号变量转换为字符串
covstr = char(cov);

% 写入到 txt 文件
fileID = fopen('cov.txt', 'w'); % 打开文件
fprintf(fileID, '%s\n', covstr);             % 写入符号表达式
fclose(fileID);                               % 关闭文件



save('cov.mat','n1nL','n1','nL','cov')








function wei=findwei(xv,n)
len=size(xv,2);
com=nchoosek(1:n, len);
nxv=sum((com-kron(ones(size(com,1),1),xv)).^2,2);
wei=find(nxv==0);
end