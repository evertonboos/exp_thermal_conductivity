function M = ort_coeff_matrix(K,mydata)
% Assembles coefficient matrix M given K, and respective variables.
% Orthotropic case.

%% Parameters
[n,N,t0,tt0,et,alpha,beta,D,x,y,P,hat_D,d0,dn,d0_e1,dn_en,X,Y] = mydata.cheby;
run('data.m')
run('data_adequation.m')

D1 = D(:,2:end-1);
D2 = D(2:end-1,:);

K11 = reshape(K(1:(n+1)^2),n+1,n+1);
K22 = reshape(K((n+1)^2+1:end),n+1,n+1);

%% Computing matrix
F = sparse((n+1)^2,(n+1)^2);
G = F;

for j = 1:n+1
    int = (j-1)*(n+1)+1:j*(n+1);
    A = l1*h1(y(j))*d0_e1 + D1*diag(K11(2:end-1,j))*D2 - l1*h2(y(j))*dn_en;
    B = l2*h3(x(j))*d0_e1 + D1*diag(K22(j,2:end-1))*D2 - l2*h4(x(j))*dn_en;
    F(int,int) = [A(1,:); D2*diag(K11(:,j))*D; A(end,:)];
    G(int,int) = [B(1,:); D2*diag(K22(j,:))*D; B(end,:)];
end

M = (1/l1^2)*F + (1/l2^2)*P*G*P';

end

