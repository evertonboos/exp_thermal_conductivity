function J = isox_Jacobian(KK,mydata)
% Jacobian in relation to 'K' (given as a vector or matrix), given the
% time steps 't0'. It's necessary to solve (n+1)^2 IVPs to compute de
% Jacobian.

%% Parameters
[n,N,t0,tt0,et,alpha,beta,D,x,y,P,hat_D,d0,dn,d0_e1,dn_en,X,Y] = mydata.cheby;
run('data.m')
run('data_adequation.m')

D1 = D(:,2:end-1);
D2 = D(2:end-1,:);
D1_null = [zeros(n+1,1), D1, zeros(n+1,1)]; % D1 with first and last columns of zeros (implementing issues)
D2_null = [zeros(1,n+1); D2; zeros(1,n+1)]; % D2 with first and last rows of zeros (implementing issues)

%% Assembling coefficient matrix M
% create a conductivity matrix from the input vector
KK = KK(:);
K = zeros(n+1,n+1);
for i = 1:n+1
    K(:,i) = KK(:);
end

% Computing matrix
F = sparse((n+1)^2,(n+1)^2);
G = F;

for j = 1:n+1
    int = (j-1)*(n+1)+1:j*(n+1);
    A = l1*h1(y(j))*d0_e1 + D1*diag(K(2:end-1,j))*D2 - l1*h2(y(j))*dn_en;
    B = l2*h3(x(j))*d0_e1 + D1*diag(K(j,2:end-1))*D2 - l2*h4(x(j))*dn_en;
    F(int,int) = [A(1,:); D2*diag(K(:,j))*D; A(end,:)];
    G(int,int) = [B(1,:); D2*diag(K(j,:))*D; B(end,:)];
end

M = (1/l1^2)*F + (1/l2^2)*P*G*P';

%% Compute solution T(K) (related to CPS discretization scheme)
TK = isox_generate_T(KK,T0,mydata);

%% Jacobian
J = zeros((n+1)^2,n+1,length(tt0));

C = c(X,Y);
C = C(:);
C = sparse(1:(n+1)^2, 1:(n+1)^2, C); % 'C' matrix, diagonal sparse

u = ones(n+1,1);
V = eye(n+1,n+1);
for jj = 1:n+1
    
    % preparing data
    v = V(:,jj);
    dK = v*u';
    PdK = reshape(P*dK(:),n+1,n+1);
    
    % Computing matrix
    F = sparse((n+1)^2,(n+1)^2);
    G = F;
    for j = 1:n+1
        int = (j-1)*(n+1)+1:j*(n+1);
        A = D1*diag(dK(2:end-1,j))*D2;
        B = D1*diag(PdK(j,2:end-1))*D2;
        F(int,int) = [A(1,:); D2*diag(dK(:,j))*D; A(end,:)];
        G(int,int) = [B(1,:); D2*diag(PdK(j,:))*D; B(end,:)];
    end
    dM1 = (1/l1^2)*F + (1/l2^2)*P*G*P';

    for k = 1:length(tt0)-1 % Crank-Nicolson's method
        % t_{i+1}
        h = tt0(k+1)-tt0(k); % \Delta t

        Mm = C - (h/2)*M;
        Mp = C + (h/2)*M;

        % right-hand side (source term W(t))
        W1a = dM1*TK(:,k); % W(t0(i)) (anterior) for K11
        W1p = dM1*TK(:,k+1); % W(t0(i+1)) (posterior) for K11

        % Crank-Nicolson's method
        J(:,jj,k+1) = Mm\( Mp*J(:,jj,k) + (h/2)*(W1a+W1p) ); % column 1 to (n+1)^2, related to K11

    end
end











% %for i = 1:n+1 % x
%     B = D1_null*D2_null;
%     hat_G = [B(1,:); D2*D; B(end,:)]; % matrix to replace at the derivative of M
%     
%     %hat_G = zeros(n+1,n+1); % matrix to replace at the derivative of M
%     
%     for j = 1:n+1 % y
%         %ell = (j-1) + (i-1)*(n+1) + 1;
%         %ell_22 = ell+(n+1)^2; % indexation for K22
%         
%         % derivative of M with respect to k_\ell        
%         A = D1_null(:,j)*D2_null(j,:);
%         F_block = [A(1,:); D2(:,j)*D(j,:); A(end,:)];
%         F = kron(speye(n+1),F_block);
%         
%         int = (j-1)*(n+1)+1:j*(n+1);
%         G = sparse((n+1)^2,(n+1)^2);
%         G(int,int) = hat_G;
% 
%         dM1 = (1/l1^2)*(F) + (1/l2^2)*(P*G*P'); % derivative matrix
% 
%         for k = 1:length(tt0)-1 % Crank-Nicolson's method
%             % t_{i+1}
%             h = tt0(k+1)-tt0(k); % \Delta t
% 
%             Mm = C - (h/2)*M;
%             Mp = C + (h/2)*M;
% 
%             % right-hand side (source term W(t))
%             W1a = dM1*TK(:,k); % W(t0(i)) (anterior) for K11
%             W1p = dM1*TK(:,k+1); % W(t0(i+1)) (posterior) for K11
% 
%             % Crank-Nicolson's method
%             J(:,j,k+1) = Mm\( Mp*J(:,j,k) + (h/2)*(W1a+W1p) ); % column 1 to (n+1)^2, related to K11
%            
%         end
%     end
%end

J = J(:,:,2:end);
% select time steps
J = J(:,:,(1:N)*(et+1));

end















