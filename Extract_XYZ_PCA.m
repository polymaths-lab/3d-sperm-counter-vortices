function [D,C,B]=Extract_XYZ_PCA(Xdata,Ydata,Zdata)

% Input: 
% Xdata (Ydata, Zdata): np*nt. 
%
% Output:
% D: sorted eigenvalues of the decomposed PCA components.
% C: sorted eigenvectors of the decomposed PCA components.
% B: time coefficients of the sorted eigenvectors.

%% Pre-process the input data

[np nt] = size(Zdata);
% new array with space 3D->1D column
clear X

    data = zeros(3*np,nt);
    for i_nt=1:nt
      % reshape matrix for a column
      data(:,i_nt) = reshape([Xdata(:,i_nt),Ydata(:,i_nt),Zdata(:,i_nt)],[3*np,1]); 
    end  

% linear principal component analysis ====================================
psi = data';
nt=size(psi,1); % number of experimental observations (e.g. frames recorded) TIME
np3=size(psi,2); % number of features recorded with each observation SPACE


% mean of all n observations (represented as (n x m)-matrix)
psi0=repmat(mean(psi),nt,1); 

% mean-corrected (n x m)-observation matrix 
Delta=psi-psi0; 


%% PCA decomposition.

C=Delta'*Delta;
[V,D]=eig(C);
[D,ind]=sort(diag(D),'descend'); % sort eigenvalues (and corresponding eigenvectors) by decreasing magnitude

V=V(:,ind);
V1=V(:,1);  % np3*1 array 
V2=V(:,2);  
V3=V(:,3);  
V4=V(:,4);
V5=V(:,5);
V6=V(:,6);
V7=V(:,7);
V8=V(:,8);
V9=V(:,9);
V10=V(:,10);
V11=V(:,11);
V12=V(:,12);
V13=V(:,13);
V14=V(:,14);
V15=V(:,15);
V16=V(:,16);
V17=V(:,17);
V18=V(:,18);
V19=V(:,19);
V20=V(:,20);
V21=V(:,21);
V22=V(:,22);
V23=V(:,23);


c0=reshape(psi0(1,:),[np,3]);
c1=reshape(V1,[np,3]);
c2=reshape(V2,[np,3]);
c3=reshape(V3,[np,3]);
c4=reshape(V4,[np,3]);
c5=reshape(V5,[np,3]);
c6=reshape(V6,[np,3]);
c7=reshape(V7,[np,3]);
c8=reshape(V8,[np,3]);
c9=reshape(V9,[np,3]);
c10=reshape(V10,[np,3]);
c11=reshape(V11,[np,3]);
c12=reshape(V12,[np,3]);
c13=reshape(V13,[np,3]);
c14=reshape(V14,[np,3]);
c15=reshape(V15,[np,3]);
c16=reshape(V16,[np,3]);
c17=reshape(V17,[np,3]);
c18=reshape(V18,[np,3]);
c19=reshape(V19,[np,3]);
c20=reshape(V20,[np,3]);
c21=reshape(V21,[np,3]);
c22=reshape(V22,[np,3]);
c23=reshape(V23,[np,3]);
C=[c0 c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 c20 c21 c22];

%coefficient of the eigenvector
B1=V1\Delta';B1=B1(:); %[1,nt]=size(B1)
B2=V2\Delta';B2=B2(:); 
B3=V3\Delta';B3=B3(:); 
B4=V4\Delta';B4=B4(:); 
B5=V5\Delta';B5=B5(:); 
B6=V6\Delta';B6=B6(:); 
B7=V7\Delta';B7=B7(:); 
B8=V8\Delta';B8=B8(:); 
B9=V9\Delta';B9=B9(:); 
B10=V10\Delta';B10=B10(:); 
B11=V11\Delta';B11=B11(:); 
B12=V12\Delta';B12=B12(:); 
B13=V13\Delta';B13=B13(:); 
B14=V14\Delta';B14=B14(:); 
B15=V15\Delta';B15=B15(:); 
B16=V16\Delta';B16=B16(:); 
B17=V17\Delta';B17=B17(:); 
B18=V18\Delta';B18=B18(:); 
B19=V19\Delta';B19=B19(:); 
B20=V20\Delta';B20=B20(:); 
B21=V21\Delta';B21=B21(:); 
B22=V22\Delta';B22=B22(:); 
B23=V23\Delta';B23=B23(:); 
B=[B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22];


a=D(1);
c=D(2);


%sum(D(1:8))/sum(D)

%% PCA reconstruction.

for t=1:length(B2)
    Xpro0(:,:,t)=c0; %np*3  
    Xpro1(:,:,t)=c0+B1(t)*c1;
    Xpro2(:,:,t)=c0+B2(t)*c2;
    Xpro3(:,:,t)=c0+B3(t)*c3;
    Xpro12(:,:,t)=c0+B1(t)*c1+B2(t)*c2; 
    Xpro123(:,:,t)=c0+B1(t)*c1+B2(t)*c2+B3(t)*c3;
    
    xm0(:,t) = Xpro0(:,1,t);
    ym0(:,t) = Xpro0(:,2,t);
    zm0(:,t) = Xpro0(:,3,t); 
    
    xm1(:,t) = Xpro1(:,1,t);
    ym1(:,t) = Xpro1(:,2,t);
    zm1(:,t) = Xpro1(:,3,t);
    
    xm2(:,t) = Xpro2(:,1,t);
    ym2(:,t) = Xpro2(:,2,t);
    zm2(:,t) = Xpro2(:,3,t);
    

    xm3(:,t) = Xpro3(:,1,t);
    ym3(:,t) = Xpro3(:,2,t);
    zm3(:,t) = Xpro3(:,3,t);
    
    xm12(:,t) = Xpro12(:,1,t);
    ym12(:,t) = Xpro12(:,2,t);
    zm12(:,t) = Xpro12(:,3,t);
    
    
    xm123(:,t) = Xpro123(:,1,t);
    ym123(:,t) = Xpro123(:,2,t);
    zm123(:,t) = Xpro123(:,3,t);
    
end

