%% Solution Section 3 article BRANDWAJN PEVA 93
% 3 classes C=3
% Load as submitted load (not attained)
function [S1, S2, S3, U1, U2, U3] = Engset_3classes (B, N1, N2, N3, b1, b2, b3, lambda1, lambda2, lambda3, mu1, mu2, mu3)

rho1 = lambda1/mu1;
rho2 = lambda2/mu2;
rho3 = lambda3/mu3;

% compute steady-state
p = zeros (N1+1,N2+1,N3+1); % +1 bc starts at 1

for i = 1:N1+1
    for j = 1:N2+1
        for l = 1:N3+1
            % check if current state is feasible
            if (i-1)*b1 + (j-1)*b2 + (l-1)*b3 <= B
                % check if product non nul
                block1 = realmin^(1/4); % 1 modified 2021/03/21
                if i > 1 % check if product is non null
                    for k = 1:i-1
                        block1 = block1 * (N1-k+1)*rho1/k;
                    end
                end
                block2 = realmin^(1/4); % 1 modified 2021/03/21
                if j > 1 % check if product is non null
                    for k = 1:j-1
                        block2 = block2 * (N2-k+1)*rho2/k;
                    end
                end
                block3 = realmin^(1/4); % 1 modified 2021/03/21
                if l > 1 % check if product is non null
                    for k = 1:l-1
                        block3 = block3 * (N3-k+1)*rho3/k;
                    end
                end                
                p(i,j,l) = block1*block2*block3;
            else
                1; % just for checking
            end
        end
    end
end

p = p./sum(sum(sum(p)));

if isnan(sum(sum(sum(p))))
    disp('Error: vector p not computable');
end

% Number of classes 1, 2 and 3 in service
U1 = 0;
for i = 2:N1+1
    U1 = U1 + (i-1)*sum(sum(p(i,:,:)));
end

U2 = 0;
for j = 2:N2+1
    U2 = U2 + (j-1)*sum(sum(p(:,j,:)));
end

U3 = 0;
for l = 2:N3+1
    U3 = U3 + (l-1)*sum(sum(p(:,:,l)));
end

% Proba of access for classes 1, 2 and 3
S1 = U1./(rho1*(N1-U1));
S2 = U2./(rho2*(N2-U2));
S3 = U3./(rho3*(N3-U3));

