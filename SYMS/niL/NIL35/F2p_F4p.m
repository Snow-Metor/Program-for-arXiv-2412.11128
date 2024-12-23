
function F2p_F4p(n)


load('canshu.mat')


%% Fnp

%%%%% Fn2

FS=F1-U1a*U1ad-U2a*U2ad...
    +trace(U1_real)*U1a_real+trace(U1_imag)*U1a_imag...
    +trace(U2_real)*U2a_real+trace(U2_imag)*U2a_imag;
Fn2=n/2*(kron(FS,eye(2*L))+kron(eye(2*L),FS))...
    +n*(n-1)*(kron(U1a_real,U1a_real)+kron(U1a_imag,U1a_imag)...
    +kron(U2a_real,U2a_real)+kron(U2a_imag,U2a_imag));


%%%% Fnp_2
com2=nchoosek(1:(2*L), 2);
Fnp=sym(zeros(size(com2, 1)));
for i1=1:size(com2, 1)
    for i2=1:size(com2, 1)
        j1=com2(i1,:);
        j2=com2(i2,:);
        s1=(j1(1)-1)*2*L+j1(2);
        s2=(j2(1)-1)*2*L+j2(2);
        s1p=(j1(2)-1)*2*L+j1(1);
        s2p=(j2(2)-1)*2*L+j2(1);
        Fnp(i1,i2)=Fn2(s1,s2)-Fn2(s1p,s2)-Fn2(s1,s2p)+Fn2(s1p,s2p);
    end
end

%%%% Fnp_n
if n>2
    for ss=3:n
        com=nchoosek(1:(2*L), ss);
        Fnpn=sym(zeros(size(com, 1)));
        for n1=1:size(com,1)
            % [num2str(ss),'/',num2str(n),';',num2str(n1),'/',num2str(size(com,1))]
            m1=com(n1,:);
            Faddn1=sym(zeros(1,size(com,1)));
            for k1=1:ss
                xv=1:(2*L);
                xv(m1(k1))=[];
                comz=nchoosek(xv,ss-1);
                for p1=1:size(comz,1)
                    nxvl=m1;
                    nxvl(k1)=[];
                    nxvr=comz(p1,:);
                    [nwl,nfl]=wei_fh(nxvl,2*L);
                    [nwr,nfr]=wei_fh(nxvr,2*L);
                    wxvr=[m1(k1),nxvr];
                    [wwr,wfr]=wei_fh(wxvr,2*L);
                    % Fnpn(n1,wwr)=Fnpn(n1,wwr)+(-1)^(k1-1)*wfr*Fnp(nwl,nwr);
                    
                    Faddn1(wwr)=Faddn1(wwr)+(-1)^(k1-1)*wfr*Fnp(nwl,nwr);
                    
                end
            end
            
            Fnpn(n1,:)=Fnpn(n1,:)+Faddn1;
            
        end
        Fnp=Fnpn;
    end
elseif n==1
    Fnp=FS;
end
Fnp=Fnp/factorial(n);

% name=['L_',num2str(L),' t_',num2str(t),' gam0_',num2str(gam0),...
%     ' gam1_',num2str(gam1),' gamL_',num2str(gamL),'.mat'];
% save(name,'n','Fnp')
if n==2
    F2p=Fnp;
    save('F2p.mat','F2p')
elseif n==4
    F4p=Fnp;
    save('F4p.mat','F4p')
end
end






function [wei,fh]=wei_fh(xv,n)
lincom=nchoosek(1:n,length(xv));
ss1=xv;
ss1(1)=[];
nn=sum((ss1-xv(1))<0);
fh=(-1)^nn;
if nn>0
    nxv=[xv(2:(nn+1)),xv(1),xv((nn+2):end)];
else
    nxv=xv;
end
len=size(lincom,1);

xun0=sum((lincom-kron(ones(len,1),nxv)).^2,2);
wei=find(xun0==0);
end


