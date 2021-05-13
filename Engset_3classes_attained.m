%% Solution Section 3 article BRANDWAJN PEVA 93
% 3 classes C=3
% Load as attained load (not submitted)
function [S1, S2, S3, lambda1, lambda2, lambda3, N1, N2, N3, X1, X2, X3] = Engset_3classes_attained (B, N1, N2, N3, b1, b2, b3, X1, X2, X3, mu1, mu2, mu3)

EPS = 10^-4; % accuracy for iteration/dichotomy
ITMAX = 50; % max nbr of iterations % 10^3, modified 27 March 2021
ITLOCALMAX = 500; % max nbr of iterations % 10^3, modified 27 March 2021 % 50 not enough
ALPHA = 0.1; % speed of dichotomy (min : 0+, max : 1)
LAMBDA_MAX = 10^2; % 10^6; modified 2021/03/21, o/w Inf in Engset_submit()

% check feasability of attained throughput
lambda1_min = 10^-6;
lambda2_min = 10^-6;
lambda3_min = 10^-6;
lambda1_max = LAMBDA_MAX; %10^6; modified 2021/03/21, o/w Inf in Engset_submit()
lambda2_max = LAMBDA_MAX; %10^6;
lambda3_max = LAMBDA_MAX; %10^6;
N1_max = 256;
N2_max = 256;
N3_max = 256;

[S1_max, S2_max, S3_max, U1_max, U2_max, U3_max] = Engset_3classes_submit(B,N1,N2,N3,b1,b2,b3,lambda1_max,lambda2_max,lambda3_max,mu1,mu2,mu3);
X1_max = (N1-U1_max)*lambda1_max * S1_max; % 2021/03/21 added (N1-U1_max)*
X2_max = (N2-U2_max)*lambda2_max * S2_max;
X3_max = (N3-U3_max)*lambda3_max * S3_max;

while (X1_max < X1 || X2_max < X2 || X3_max < X3) && (N1 < N1_max && N2 < N2_max && N3 < N3_max)
    if X1_max < X1
        %lambda1_max = lambda1_max * 10;
        % better to increase N1
        N1 = N1+2;
    end
    if X2_max < X2        
        %lambda2_max = lambda2_max * 10;
        N2 = N2+2;
    end
    if X3_max < X3
        %lambda3_max = lambda3_max * 10;
        N3 = N3+2;
    end
    [S1_max, S2_max, S3_max, U1_max, U2_max, U3_max] = Engset_3classes_submit(B,N1,N2,N3,b1,b2,b3,lambda1_max,lambda2_max,lambda3_max,mu1,mu2,mu3);
    X1_max = (N1-U1_max)*lambda1_max * S1_max; % 2021/03/21 added (N1-U1_max)*
    X2_max = (N2-U2_max)*lambda2_max * S2_max; 
    X3_max = (N3-U3_max)*lambda3_max * S3_max; 
end

if false % condition should be when too many states for Engset_3classes_submit 
    disp('Unattainable level of X1, X2 and X3');
    S1 = 0;
    S2  = 0;
    S3  = 0;
    lambda1  = 0;
    lambda2 = 0;
    lambda3 = 0;
    %return
end

% iterate on the two classes
delta = 1;
it = 1;
lambda1_cur = lambda1_min;
lambda2_cur = lambda2_min;
lambda3_cur = lambda3_min;
while delta > 100*EPS && it < ITMAX % accuracy cannot be the same as in dichotomy
    
    % keep classes 2 and 3 constant, find rate of class 1
    eps1 = 1;
    lambda1_min = 10^-6;
    lambda1_max = LAMBDA_MAX; %10^6; modified 2021/03/21, o/w Inf in Engset_submit()
    % check feasability
    [S1_max, S2_max, S3_max, U1_max, U2_max, U3_max] = Engset_3classes_submit(B,N1,N2,N3,b1,b2,b3,lambda1_max,lambda2_cur,lambda3_cur,mu1,mu2,mu3);
    X1_max = N1*lambda1_max * S1_max;
    if X1_max < X1
        disp('Unattainable level of X1, X2 and X3');
        S1 = 0;
        S2  = 0;
        S3  = 0;
        lambda1  = 0;
        lambda2 = 0;
        lambda3 = 0;
        return
    end
    itlocal = 0;
    while eps1 > EPS && itlocal < ITLOCALMAX
        lambda1_cur = (lambda1_min + lambda1_max)/2;
        
        [S1_cur, S2_cur, S3_cur, U1_cur, U2_cur, U3_cur] = Engset_3classes_submit(B,N1,N2,N3,b1,b2,b3,lambda1_cur,lambda2_cur,lambda3_cur,mu1,mu2,mu3);
        
        if (N1-U1_cur)*lambda1_cur*S1_cur < X1
            lambda1_min = ALPHA*lambda1_cur+(1-ALPHA)*lambda1_min; % TMP, should be lambda1_cur
        else
            lambda1_max = ALPHA*lambda1_cur+(1-ALPHA)*lambda1_max; % TMP, should be lambda1_cur
        end
        lambda1_cur = (lambda1_min + lambda1_max)/2;
        eps1 = abs(X1 - (N1-U1_cur)*lambda1_cur*S1_cur);
        eps1tab(it) = eps1; % TMP
        lambda1tab(it) = lambda1_cur; % TMP
        itlocal = itlocal + 1;
    end
    
    % keep classes 1 and 3 constant, find rate of class 2
    eps2 = 1;
    lambda2_min = 10^-6;
    lambda2_max = LAMBDA_MAX; %10^6; modified 2021/03/21, o/w Inf in Engset_submit()
    % check feasability
    [S1_max, S2_max, S3_max, U1_max, U2_max, U3_max] = Engset_3classes_submit(B,N1,N2,N3,b1,b2,b3,lambda1_cur,lambda2_max,lambda3_cur,mu1,mu2,mu3);
    X2_max = N2*lambda2_max * S2_max;
    if X2_max < X2
        disp('Unattainable level of X1, X2 and X3');
        S1 = 0;
        S2  = 0;
        S3 = 0;
        lambda1  = 0;
        lambda2 = 0;
        lambda3 = 0;
        return
    end
    itlocal = 0;    
    while eps2 > EPS && itlocal < ITLOCALMAX
        lambda2_cur = (lambda2_min + lambda2_max)/2;
        
        [S1_cur, S2_cur, S3_cur, U1_cur, U2_cur, U3_cur] = Engset_3classes_submit(B,N1,N2,N3,b1,b2,b3,lambda1_cur,lambda2_cur,lambda3_cur,mu1,mu2,mu3);
        
        if (N2-U2_cur)*lambda2_cur*S2_cur < X2
            lambda2_min = ALPHA*lambda2_cur+(1-ALPHA)*lambda2_min; % TMP, should be lambda2_cur
        else
            lambda2_max = ALPHA*lambda2_cur+(1-ALPHA)*lambda2_max; % TMP, should be lambda2_cur
        end
        lambda2_cur = (lambda2_min + lambda2_max)/2;
        eps2 = abs(X2 - (N2-U2_cur)*lambda2_cur*S2_cur);
        eps2tab(it) = eps2;  % TMP
        lambda2tab(it) = lambda2_cur;   % TMP     
        itlocal = itlocal + 1;        
    end
    
	% keep classes 1 and 2 constant, find rate of class 3
    eps3 = 1;
    lambda3_min = 10^-6;
    lambda3_max = LAMBDA_MAX; %10^6; modified 2021/03/21, o/w Inf in Engset_submit()
    % check feasability
    [S1_max, S2_max, S3_max, U1_max, U2_max, U3_max] = Engset_3classes_submit(B,N1,N2,N3,b1,b2,b3,lambda1_cur,lambda2_cur,lambda3_max,mu1,mu2,mu3);
    X3_max = N3*lambda3_max * S3_max;
    if X3_max < X3
        disp('Unattainable level of X1, X2 and X3');
        S1 = 0;
        S2  = 0;
        S3 = 0;
        lambda1  = 0;
        lambda2 = 0;
        lambda3 = 0;
        return
    end
    itlocal = 0;    
    while eps3 > EPS && itlocal < ITLOCALMAX
        lambda3_cur = (lambda3_min + lambda3_max)/2;
        
        [S1_cur, S2_cur, S3_cur, U1_cur, U2_cur, U3_cur] = Engset_3classes_submit(B,N1,N2,N3,b1,b2,b3,lambda1_cur,lambda2_cur,lambda3_cur,mu1,mu2,mu3);
        
        if (N3-U3_cur)*lambda3_cur*S3_cur < X3
            lambda3_min = ALPHA*lambda3_cur+(1-ALPHA)*lambda3_min; % TMP, should be lambda2_cur
        else
            lambda3_max = ALPHA*lambda3_cur+(1-ALPHA)*lambda3_max; % TMP, should be lambda2_cur
        end
        lambda3_cur = (lambda3_min + lambda3_max)/2;
        eps3 = abs(X3 - (N3-U3_cur)*lambda3_cur*S3_cur);
        eps3tab(it) = eps3;  % TMP
        lambda3tab(it) = lambda3_cur;   % TMP     
        itlocal = itlocal + 1;        
    end
    
    % update values of eps1 and compute stopping criterion
    eps1 = abs(X1 - (N1-U1_cur)*lambda1_cur*S1_cur);
    eps2 = abs(X2 - (N2-U2_cur)*lambda2_cur*S2_cur);    
    delta = max([eps1,eps2,eps3]);
    it = it + 1;
end

if it == ITMAX
   disp('Convergence not found'); 
end

% figure();
% plot(1:length(eps1tab),eps1tab,'-r',1:length(eps1tab),eps2tab,'-g') % TMP

lambda1 = lambda1_cur;
lambda2 = lambda2_cur;
lambda3 = lambda3_cur;
X1 = (N1-U1_cur)*lambda1_cur*S1_cur;
X2 = (N2-U2_cur)*lambda2_cur*S2_cur;
X3 = (N3-U3_cur)*lambda3_cur*S3_cur;
[S1, S2,S3, U1, U2, U3] = Engset_3classes_submit(B,N1,N2,N3,b1,b2,b3,lambda1,lambda2,lambda3,mu1,mu2,mu3);