clear all;
close all;

%% Example for Mix 1 with 2 PMs 
%% Consider 19 levels of loads and 13 penalty costs for Disaggregated architecture
C = 2;
B = 128;
b1 = 2;
b2 = 4;
b3 = 8;
wload = [5 6.25 7.5 8.75 10 11.25 12.5 15 17.5 18.75 20 20.625 21.25 22.5 23.75 25 26.25 27.5 28.75];
mu1 = 1;
mu2 = 1;
mu3 = 1;
p1 = 1/2;
p2 = 1/4;
p3 = 1/4;
penalty = [1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.1 2.2];
pen1 = 1;
pen2 = 2;
pen3 = 3;

% added March 29 2021 to unify parameters
N1_init = round(B/12);
N2_init = round(B/12);
N3_init = round(B/12);
B = 64*C;

% Disaggregated architecture
for j = 1:length(penalty)
    imax(j) = length(wload);
    for i = 1:length(wload)
        [S1(j,i), S2(j,i), S3(j,i), lambda1(j,i), lambda2(j,i), lambda3(j,i),...
            N1(j,i), N2(j,i), N3(j,i), X1(j,i), X2(j,i), X3(j,i)] ...
            = Engset_3classes_attained (...
            B, N1_init, N2_init, N3_init, b1, b2, b3, ...
            wload(i)*p1, wload(i)*p2, wload(i)*p3, ...
            mu1/penalty(j), mu2/penalty(j), mu3/penalty(j));
        if S1(j,i) < 0.8 || S2(j,i) < 0.8 || S3(j,i) < 0.8
            imax(j) = i;
            break;
        end
    end
end

% Classic architecture
imaxc = length(wload);
for i = 1:length(wload)
    [S1c(i), S2c(i), S3c(i), lambda1c(i), lambda2c(i), lambda3c(i),...
        N1c(i), N2c(i), N3c(i), X1c(i), X2c(i), X3c(i)] ...
        = Engset_3classes_attained (...
        B/C, N1_init, N2_init, N3_init, b1, b2, b3, ...
        wload(i)*p1/C, wload(i)*p2/C, wload(i)*p3/C, ...
        mu1, mu2, mu3);
    if S1c(i) < 0.8 || S2c(i) < 0.8 || S3c(i) < 0.8
        imaxc = i;
        break;
    end
end

%% Compute performance based on class 3 requirement
target = 0.99;

% Disaggregated architecture
for j = 1:length(penalty)
    [SS3, index] = unique(S3(j,1:imax(j))); % need to remove duplicate values for interpolation
    Res3(j) = interp1(SS3,X3(j,index),target,'linear');
    Res1(j) = p1/p3*Res3(j);
    Res2(j) = p2/p3*Res3(j);
end

% Classic architecture
[SS3c, index] = unique(S3c(1:imaxc));
Res3c = C*interp1(SS3c,X3c(index),target,'linear','extrap'); % extrap to extrapolate beyond ranged interval of wload, adde C* April 1st 2021
Res1c = p1/p3*Res3c;
Res2c = p2/p3*Res3c;


disp(['Max Thr for an acceptance rate of ' num2str(target)]);
First_col = [9 9 'Control '];
for j = 1:length(penalty)
    First_col = [First_col ' Pen' num2str(j) ':' num2str(penalty(j))];
end
disp(First_col)
disp(['Class 1 ' 9 num2str(Res1c) 9 num2str(Res1)]);
disp(['Class 2 ' 9 num2str(Res2c) 9 num2str(Res2)]);
disp(['Class 3 ' 9 num2str(Res3c) 9 num2str(Res3)]);
disp(['All Classes ' 9 num2str(Res1c + Res2c + Res3c) 9 num2str(Res1 + Res2 + Res3)]);




