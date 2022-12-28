function J = ort_Jacobian(K,mydata)
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
M = ort_coeff_matrix(K,mydata);

%% Compute solution T(K) (related to CPS discretization scheme)
TK = ort_generate_T(K,T0,mydata);

%% Jacobian
J = zeros((n+1)^2,2*(n+1)^2,length(tt0));

C = c(X,Y);
C = C(:);
C = sparse(1:(n+1)^2, 1:(n+1)^2, C); % 'C' matrix, diagonal sparse

for i = 1:n+1 % x
    %hat_G = hat_D(:,i)*D(i,:); % matrix to replace at the derivative of M
    B = D1_null(:,i)*D2_null(i,:);
    hat_G = [B(1,:); D2(:,i)*D(i,:); B(end,:)];
    for j = 1:n+1 % y
        ell = (j-1) + (i-1)*(n+1) + 1;
        ell_22 = ell+(n+1)^2; % indexation for K22
        
        % derivative of M with respect to k_\ell
        F = sparse((n+1)^2,(n+1)^2);
        G = F;
        
        int = (i-1)*(n+1)+1:i*(n+1);
        A = D1_null(:,j)*D2_null(j,:);
        %F(int,int) = hat_D(:,j)*D(j,:);
        F(int,int) = [A(1,:); D2(:,j)*D(j,:); A(end,:)];
        int = (j-1)*(n+1)+1:j*(n+1);
        G(int,int) = hat_G;

        dM1 = (1/l1^2)*(F); % derivative matrix for K11
        dM2 = (1/l2^2)*(P*G*P'); % derivative matrix for K22

        for k = 1:length(tt0)-1 % Crank-Nicolson's method
            % t_{i+1}
            h = tt0(k+1)-tt0(k); % \Delta t

            Mm = C - (h/2)*M;
            Mp = C + (h/2)*M;

            % right-hand side (source term W(t))
            W1a = dM1*TK(:,k); % W(t0(i)) (anterior) for K11
            W1p = dM1*TK(:,k+1); % W(t0(i+1)) (posterior) for K11
            W2a = dM2*TK(:,k); % W(t0(i)) (anterior) for K22
            W2p = dM2*TK(:,k+1); % W(t0(i+1)) (posterior) for K22

            % Crank-Nicolson's method
            J(:,ell,k+1) = Mm\( Mp*J(:,ell,k) + (h/2)*(W1a+W1p) ); % column 1 to (n+1)^2, related to K11
            J(:,ell_22,k+1) = Mm\( Mp*J(:,ell_22,k) + (h/2)*(W2a+W2p) ); % column (n+1)^2+1 to 2*(n+1)^2, related to K22

        end
    end
end

J = J(:,:,2:end);
% select time steps
J = J(:,:,(1:N)*(et+1));

end















