function [C1,eL2max,tvec]=gcqrk_PDE(sol,f,Kfun,RK,N,Tf,grad,c_u,node,elem,Dirichlet)

% approximate D_t^\alpha u+u=f(f is a given function) at prescribed general time mesh tn, n=1,...,N
%09-05-2022
%It uses gCQ based on the Euler method
% Written by Jing Guo Aug. 2023.
%%
% isBdNode = false(N,1); 
% if ~isempty(Dirichlet)
%     isBdNode(Dirichlet(:)) = true;
% else % Pure Neumann boundary condition
%     isBdNode(1) = true;
% end
% freeNode = find(~isBdNode);
%Approximation of convolution at first N steps
[s,A,c,b,eiginvA,V,invV]=setsolver(RK); %Buther tableau and eigenvalue decomposition of A^(-1)
I=speye(s); %estimated order of convergence plus one
tvec = ((0:N)*Tf^(1/grad)/N).^grad; %time nodes
hvec = tvec(2:end)-tvec(1:end-1); %temporal step

%Parameter for space
J= size(node,1);
Ij=speye(J);
IL=speye(s*J);
% L=lenghth(freenode);
C = zeros(J*s,N+1);
errL2 = zeros(1,N+1);
%% Assemble stiffness and fmass matrix
ve(:,:,3) = node(elem(:,2),:) - node(elem(:,1),:);
ve(:,:,1) = node(elem(:,3),:) - node(elem(:,2),:);
ve(:,:,2) = node(elem(:,1),:) - node(elem(:,3),:);
area = 0.5*abs(-ve(:,1,3).*ve(:,2,2)+ve(:,2,3).*ve(:,1,2));

[Stf,Mas] = assemblematrix(node,elem);

center = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
C(:,1)=repmat(accumarray(elem(:),repmat((sol(center(:,1),center(:,2),tvec(1))).*area/3,3,1),[J,1]),s,1);
errL2(1)=getL2error_j(node,elem,@(x,y) sol(x,y,tvec(1)),C(((s-1)*J+1):(s*J),1));

%% Assembling right side
% center = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
if RK==1
    circ_radi=3;%choice of the radius of the circle contour
else
    circ_radi=1;
end

lamda=max(abs(eiginvA));
switch circ_radi
    case 1
        m=min(real(eiginvA))/max(hvec(1:N));
        %  m=1/max(hvec(1:N));
        %         M=(2*s-1)*lamda/min(hvec(1:N));
        M=5*lamda/min(hvec(1:N));
    case 2
        m=sqrt(6)/max(hvec(1:N));
        M=sqrt(6)/min(hvec(1:N));
    case 3
        m=1/max(hvec(1:N));
        M=1/min(hvec(1:N));
end
% nq=max(50,ceil(N^2));
nq=N^2;
[Z,weights]=elliptic_quadrature(m,M,nq);

isBdNode = false(J,1); 
isBdNode(Dirichlet(:)) = true;
freeNode = find(~isBdNode);
K = Kfun(Z);
Omega0 = Kfun(eiginvA/hvec(1));
Omega0 = kron((A*hvec(1))\(V*diag(Omega0)*invV),Mas(freeNode,freeNode));
%% Dirichlet boundary conditions


%% Solve 
rhs=[];
for k=1:s
rhsk = accumarray(elem(:),repmat(f(center(:,1),center(:,2),tvec(1)+c(k)*hvec(1)).*area/3,3,1),[J,1]);
rhs=[rhs;rhsk];
end
% FreeV=kron((1:s)',freeNode);
freeV=freeNode;%generate free note vector for RK
for vk=1:s-1
  freeV=[freeV; vk*J+freeNode];
end

C(freeV,2)=(c_u*kron(I,Stf(freeNode,freeNode))+Omega0)\(rhs(freeV)+ Omega0*kron(ones(s,1),C((s-1)*J+freeNode,1)));
errL2(2)=getL2error_j(node,elem,@(x,y) sol(x,y,tvec(2)),C(((s-1)*J+1):(s*J),2));
U=zeros(s*J,length(Z));
% rhs=f(center(:,1),center(:,2),tvec(1)+c*hvec(1))+Omega0*C(:,end,1)*ones(s,1);%how to write this?? 01.08.2023
% rhs= accumarray(elem(:),repmat((f(center(:,1),center(:,2)tvec(1)+c*hvec(1)))...
%     .*area/3,3,1),[J,1]);



%Parameters for elliptic quadrature


for ll=3:N+1
    Omega0 = Kfun(eiginvA/hvec(ll-1));
    Omega0 =kron((A*hvec(ll-1))\(V*diag(Omega0)*invV),Mas(freeNode,freeNode));
    %integral along circular contour
    int = zeros(length(freeV),1);
     for mm=1:length(Z)
%         U(FreeV,mm)=kron(I-hvec(ll-2)*Z(mm)*A,Ij(freeNode,freeNode))\(kron(ones(s,1),U(s*freeNode,mm))+...
%             C(FreeV,ll-1)-kron(ones(s,1),C(s*freeNode,ll-2)));
%      U(FreeV,mm)=kron(I-hvec(ll-2)*Z(mm)*A,Mas(freeNode,freeNode))\(kron(I,Mas(freeNode,freeNode))*kron(ones(s,1),U(s*freeNode,mm))+...
%             kron(I,Mas(freeNode,freeNode))*(C(FreeV,ll-1)-kron(ones(s,1),C(s*freeNode,ll-2))));
%      int=int+kron(I-hvec(ll-1)*Z(mm)*A,Mas(freeNode,freeNode))\...
%             (weights(mm)*K(mm)* kron(I,Mas(freeNode,freeNode))*kron(ones(s,1),U(s*freeNode,mm)));
     U(freeV,mm)=kron(I-hvec(ll-2)*Z(mm)*A,Ij(freeNode,freeNode))\(kron(ones(s,1),U((s-1)*J+freeNode,mm))+...
            kron(I,Mas(freeNode,freeNode))*(C(freeV,ll-1)-kron(ones(s,1),C((s-1)*J+freeNode,ll-2))));
     int=int+kron(I-hvec(ll-1)*Z(mm)*A,Ij(freeNode,freeNode))\...
            (weights(mm)*K(mm)*kron(ones(s,1),U((s-1)*J+freeNode,mm)));
     end
rhs=[];
for k=1:s
rhs_k = accumarray(elem(:),repmat(f(center(:,1),center(:,2),tvec(ll-1)+c(k)*hvec(ll-1)).*area/3,3,1),[J,1]);
rhs=[rhs;rhs_k];
end
C(freeV,ll)=(Omega0+c_u*kron(I,Stf(freeNode,freeNode)))\(rhs(freeV)-int+Omega0*kron(ones(s,1),C((s-1)*J+freeNode,ll-1)));
errL2(ll)=getL2error_j(node,elem,@(x,y) sol(x,y,tvec(ll)),C(((s-1)*J+1):s*J,ll));
end
C1=C(((s-1)*J+1):(s*J),:);
% % %plot the stability region, poles and quadrature contour
% gamcirc = @(t,mid,r) mid+r*exp(1i*t);
% symbol=@(z) inv(A+z/(1-z)*ones(s,1)*b'); %Banjai 2010
% 
% nstab = 100;
% nodstab = -pi:2*pi/nstab:pi;
% nodstab = nodstab(1:end-1);
% zstab = gamcirc(nodstab,0,.99); %nodes on the circle
% F=[];
% for ll=1:length(nodstab)
%     freq = eig(symbol(zstab(ll)));
%     F=[F; freq];
% end
% figure(10); clf;
% plot(F/min(hvec(1:N)),'k.'); hold on
% for ll=1:s
%     plot(real(eiginvA(ll)./hvec(1:N)),imag(eiginvA(ll)./hvec(1:N)),'r*')%poles
% end
% plot(Z,'b.-'); hold off

%Quadrature for the memory term

ErrM=[];%err  matrix
eL2max=max(errL2);
%Error
% for k=1:N
% solvec =accumarray(elem(:),repmat((sol(center(:,1),center(:,2),tvec(k+1))).*area/3,3,1),[J,1]);
% err = C(s*(1:J),k+1)-solvec;
% ErrM=[ErrM, err];
% end

% [eamax]=max(errV);
end
