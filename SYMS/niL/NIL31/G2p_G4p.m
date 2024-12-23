
clc
clear

load('canshu.mat')



Fgam=sym(zeros(2*L));

for s1=1:2*L
    for s2=1:2*L
        Fgam(s1,s2)=...
            (U1a_real(s1,s2)*(0.5*(U1_real(s1,s1)+U1_real(s2,s2))-trace(U1_real)))...
            +(U1a_imag(s1,s2)*(0.5*(U1_imag(s1,s1)+U1_imag(s2,s2))-trace(U1_imag)))...
            +(U2a_real(s1,s2)*(0.5*(U2_real(s1,s1)+U2_real(s2,s2))-trace(U2_real)))...
            +(U2a_imag(s1,s2)*(0.5*(U2_imag(s1,s1)+U2_imag(s2,s2))-trace(U2_imag)));
    end
end


FNB2=Fgam+1i*FB;


com2=nchoosek(1:(2*L),2);
G2p=sym(zeros(size(com2,1),1));
for s3=1:size(com2,1)
    j1=com2(s3,1);
    j2=com2(s3,2);
    G2p(s3)=FNB2(j1,j2)-FNB2(j2,j1);
end

com4=nchoosek(1:2*L,4);
FNB4=2*Fgam+6*1i*FB;
G4p=sym(zeros(size(com4,1),size(com2,1)));
for i1=1:size(com4,1)
    xv1=com4(i1,:);
    for i2=1:4
        k1=xv1(i2);
        xv2=xv1;
        xv2(i2)=[];
        fu1=(-1)^(i2-1);
        for i3=1:3
            xv3=xv2;
            k2=xv3(i3);
            xv3(i3)=[];
            fu2=(-1)^(i3-1);
            if k1<k2
            nxv=[k1,k2];
            wei=findwei(nxv,2*L);
            G4p(i1,wei)=G4p(i1,wei)+fu1*fu2*(FNB4(xv3(1),xv3(2))-FNB4(xv3(2),xv3(1)));
            end
        end
    end
end
  
G2p=G2p/2;
G4p=G4p/4/3;
save('G2p_G4p.mat','G2p','G4p')

function wei=findwei(xv,n)
len=size(xv,2);
com=nchoosek(1:n, len);
nxv=sum((com-kron(ones(size(com,1),1),xv)).^2,2);
wei=find(nxv==0);
end

