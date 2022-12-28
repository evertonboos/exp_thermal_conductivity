function TK = ort_generate_T(K,T0,mydata)
% calculate temperature T at time step N given a vector K=[K11(:); K22(:)]
% containing the heat conductivity parameter.

%% Parameters
[n,N,t0,tt0,et,alpha,beta,D,x,y,P,hat_D,d0,dn,d0_e1,dn_en,X,Y] = mydata.cheby;
run('data.m')
run('data_adequation.m')

%% Assembling coefficient matrix
M = ort_coeff_matrix(K,mydata);

%% IVP
TK = zeros((n+1)^2,length(tt0));
TK(:,1) = T0;

C = c(X,Y);
%C = flip(C);
C = C(:);
C = sparse(1:(n+1)^2, 1:(n+1)^2, C); % 'C' matrix, diagonal sparse

for i = 1:length(tt0)-1 % Crank-Nicolson's method
    
    % t_{i+1}
    h = tt0(i+1)-tt0(i); % \Delta t
    
    Mm = C - (h/2)*M;
    Mp = C + (h/2)*M;
    %Mm = speye((n+1)^2) - (h/2)*M;
    %Mp = speye((n+1)^2) + (h/2)*M;
    
    % right-hand side (vector S(t))
    Sa = ort_coeff_vector(tt0(i),mydata); % S(t_i) (anterior)
    Sp = ort_coeff_vector(tt0(i+1),mydata); % S(t_{i+1}) (posterior)
    
    % Crank-Nicolson's method
    TK(:,i+1) = Mm\( Mp*TK(:,i) + (h/2)*(Sa+Sp) );
end

%TK = TK(2:end);
% select time steps
%TK = TK(:,(1:N)*(et+1));

end

