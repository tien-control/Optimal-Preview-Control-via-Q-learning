%% On policy
clc; clear; close all;
%% Simple system
N = 4;
% (1)
A{1} = [-0.6 0.4; -1.2 1.2]; B{1} = [1.1; 0.6]; D{1} = [1; 1];
Q{1} = eye(2); R{1} = eye(1);
% (2)
A{2} = [-0.8 0.4; -1.2 -0.2]; B{2} = [1.3; 0.7]; D{2} = [-1; 1];
Q{2} = eye(2); R{2} = eye(1);
% (3)
A{3} = [0.5 -0.7; 0.3 1.2]; B{3} = [0.5; -0.5]; D{3} = [1; -1];
Q{3} = eye(2); R{3} = eye(1);
% (4)
A{4} = [-0.6 1.4; -1.2 -0.4]; B{4} = [1.1; 0.7]; D{4} = [-1; 1];
Q{4} = eye(2); R{4} = eye(1);
% x_0
x_0 = [0.5 0.5]';
%% Parameter
a = 5; t_end = 100;
n = size(A{1},2); m = size(B{1},2); l = size(D{1},2);
%% Signals
% Disturbance
for i = 1:t_end
        if i >= 30 && i<= 70    d{i} = 0.07*sin((1/15)*pi*i);
        else                    d{i} = 0;
        end
end
%% System trans
Psin = []; Lambdan = []; Qsn = []; Rsn = [];
Ad = zeros((a+N+1)*l); Ad(1:end-l,l+1:end) = eye((a+N)*l);
% Trans1
for i = 1:N
    G{i} = [D{i} zeros(n,(a+N)*l)];
    Psi{i} = [A{i} G{i}; zeros((a+N+1)*l,n) Ad];
    Lambda{i} = [B{i}; zeros((a+N+1)*l,m)];
    Qs{i} = blkdiag(Q{i}, zeros((a+N+1)*l));
    Rs{i} = R{i};
end
    X_0 = [x_0; zeros((a+N+1)*l,1)];
% Trans2
M0 = eye(n+(a+N+1)*l);
for i = 1:N
    M0 = Psi{i}*M0;
    M1 = eye(n+(a+N+1)*l); M2 = [];
    for j = 1:N
        if j < i
            M2 = [M2; zeros(n+(a+N+1)*l,m)];
        elseif j == i
            M2 = [M2; Lambda{i}];
        else
            M1 = Psi{j}*M1;
            M2 = [M2; M1*Lambda{i}];
        end
    end
    Psin = [Psin; M0]; Lambdan = [Lambdan, M2];
    Qsn = blkdiag(Qsn,Qs{i}); Rsn = blkdiag(Rsn,Rs{i});
end
Psin1 =  Psin(1:(N-1)*(n+(a+N+1)*l),:); Psin2 = Psin((N-1)*(n+(a+N+1)*l)+1:N*(n+(a+N+1)*l),:);
Psin = [zeros(N*(n+(a+N+1)*l),(N-1)*(n+(a+N+1)*l)),Psin];
Lambdan1 = Lambdan(1:(N-1)*(n+(a+N+1)*l),:); Lambdan2 = Lambdan((N-1)*(n+(a+N+1)*l)+1:N*(n+(a+N+1)*l),:);
Qsn1 = Qsn(1:(N-1)*(n+(a+N+1)*l),1:(N-1)*(n+(a+N+1)*l)); Qsn2 = Qsn((N-1)*(n+(a+N+1)*l)+1:N*(n+(a+N+1)*l),(N-1)*(n+(a+N+1)*l)+1:N*(n+(a+N+1)*l));
Rsn1= Rsn(1:(N-1)*m,1:(N-1)*m); Rsn2 = Rsn((N-1)*m+1:N*m,(N-1)*m+1:N*m);
    Xn_0 = [zeros((N-1)*(n+(a+N+1)*l),1); X_0];
%% On policy
% start controller
%[K,P,~] = dlqr(Psin,Lambdan,Qsn,Rsn);
An = []; Bn = []; Qn = []; Rn = [];
M0 = eye(n); %(M0,M1,M2) are medium matrix
% Process follow column
for i = 1:N
    M0 = A{i}*M0;
    M1 = eye(n); M2 = [];
    for j = 1:N
        if j < i
            M2 = [M2; zeros(n,m)];
        elseif j == i
            M2 = [M2; B{i}];
        else
            M1 = A{j}*M1;
            M2 = [M2; M1*B{i}];
        end
    end
    An = [An; M0]; Bn = [Bn, M2];
    Qn = blkdiag(Qn,Q{i}); Rn = blkdiag(Rn,R{i});
end
An = [zeros(N*n,(N-1)*n),An];
% Optimal controller
[K,P,~] = dlqr(An,Bn,Qn,Rn);
K = 0.5*K(:,size(K,2)*(N-1)/N+1:end);
K = [K, zeros(m*N,(a+N+1)*l)];

%% Thuat toan Q_learning_on_policy
% Khoi tao bien
xM = Xn_0;
uuM = [];
xxM = [];

L = K;

L_s = [];
Norm_L(1)=norm(L);
for j = 1:6
    Z = []; Y = [];
    L_s(:,:,j) = L;
    L=[zeros(m*N,(n+(a+N+1)*l)*(N-1)),L]; %Create full L
 
    for i = 1:800
        xold=xM;
        uM = -L*xold + 0.01*randn(N*m,1);
        uuM = [uuM uM];
        xM = Psin*xold + Lambdan*uM + 0.001*randn((n+(a+N+1)*l)*N,1);
        xxM = [xxM xM];

        xold1 = xold(1:end-(n+(a+N+1)*l));
        xold2 = xold(end-(n+(a+N+1)*l)+1:end);
        % Lay x1, x2
        xM1 = xM(1:end-(n+(a+N+1)*l));
        xM2 = xM(end-(n+(a+N+1)*l)+1:end);
        
        pk = [xold2;uM];
        phi_1 = vecv(pk);
        Jn = xold2'*Qsn2*xold2 + uM'*Rsn*uM + xM1'*Qsn1*xM1; %r bi thay doi
        uM = -L*xM;
        
        pk1=[xM2;uM];
        phi_2=vecv(pk1);
        
        phi=phi_1-phi_2;
        
        Z=[Z;phi'];
        Y=[Y;Jn];


    end
    rank(Z'*Z)
    theta = inv(Z'*Z)*Z'*Y;
    H = vecs_inv(theta);
    H22 = H(end-m*N+1:end,end-m*N+1:end);
    H12 = H(1:(n+(a+N+1)*l),end-m*N+1:end);
    rank(H)
    L = inv(H22)*H12';

    Norm_L(j+1) = norm(L);
    Norm_H(j+1) = norm(H);
    
    %if (abs(Norm_L(j+1) - Norm_L(j)) <= 0.05)
    %    break;
    %end
end
    L_s(:,:,j+1) = L;
L=[zeros(m*N,(n+(a+N+1)*l)*(N-1)),L];
for i=1:500
    xM=Psin*xM - Lambdan*L*xM;
    xxM=[xxM xM];
    uuM=[uuM -L*xM];
end

%% Ve do thi
% Qua trinh Q_learning
figure(1)
%subplot(1,2,1);
plot(0:length(Norm_L)-1,Norm_L,'Color', [0 0.4470 0.7410],'LineStyle', '-', 'LineWidth', 1,'Marker', '*');
ylabel('$||K||$', 'FontSize', 10, 'Interpreter', 'latex');
%{
subplot(1,2,2);
plot(0:length(Norm_H)-1,Norm_H,'Color', [1, 0, 0],'LineStyle', '-', 'LineWidth', 1,'Marker', '*');
ylabel('$||H||$', 'FontSize', 10, 'Interpreter', 'latex');
%}
%{
figure(2)
subplot(2,1,1);
plot(0:length(xxM)-1,xxM(1:2,:));
xlabel('Time step ($k$)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
legend('$x_1$','$x_2$', 'Interpreter', 'latex', 'FontSize', 10);
subplot(2,1,2);
plot(0:length(uuM)-1,uuM(1,:));
xlabel('Time step ($k$)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$u$', 'Interpreter', 'latex', 'FontSize', 12);
legend('$u$', 'Interpreter', 'latex', 'FontSize', 10);
%}

%% simulation of controller
% Optimal controller
for o = 1:4
    x = {}; u = {};
    if o<=3 
        K = L_s(:,:,o);
    else
        K = L_s(:,:,size(L_s,3));
    end
    % Simulation
    x{1} = x_0;
    for i = 1:t_end-1
        k = mod(i,N)+1; % State matrix
        % Find dot control signal
        ud = zeros(l,1);
        for j = 0:a
            if i+j <= t_end
                    ud = ud + K(k,n+j*l+1:n+(j+1)*l)*d{i+j};
            else    ud = ud;
            end
        end
        u{i} = - (K(k,1:n)*x{i} + ud);
        % State value
        x{i+1} = A{k}*x{i} + B{k}*u{i} + D{k}*d{i};
    end
    xx{o} = cell2mat(x);
    uu{o} = cell2mat(u);
end

% Tinh trang cac bien trang thai qua cac bo dieu khien
figure(2)
% Bo dieu khien 1
subplot(4, 2, 1);
hold on
plot(0:size(xx{1},2)-1,xx{1}(1,:))
plot(0:size(xx{1},2)-1,xx{1}(2,:))
ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 10);
legend('$x_1$','$x_2$','show', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 6);

subplot(4, 2, 2);
hold on
plot(0:size(uu{1},2)-1,uu{1}(1,:))
ylabel('$u$', 'Interpreter', 'latex', 'FontSize', 10);
%legend('$u$','Interpreter', 'latex', 'FontSize', 10);

% Bo dieu khien 2
subplot(4, 2, 3);
hold on
plot(0:size(xx{2},2)-1,xx{2}(1,:))
plot(0:size(xx{2},2)-1,xx{2}(2,:))
ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 10);
legend('$x_1$','$x_2$','show', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 6);

subplot(4, 2, 4);
hold on
plot(0:size(uu{2},2)-1,uu{2}(1,:))
ylabel('$u$', 'Interpreter', 'latex', 'FontSize', 10);
%legend('$u$','Interpreter', 'latex', 'FontSize', 10);

% Bo dieu khien 3
subplot(4, 2, 5);
hold on
plot(0:size(xx{3},2)-1,xx{3}(1,:))
plot(0:size(xx{3},2)-1,xx{3}(2,:))
ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 10);
legend('$x_1$','$x_2$','show', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 6);

subplot(4, 2, 6);
hold on
plot(0:size(uu{3},2)-1,uu{3}(1,:))
ylabel('$u$', 'Interpreter', 'latex', 'FontSize', 10);
%legend('$u$','Interpreter', 'latex', 'FontSize', 10);


% Bo dieu khien cuoi
subplot(4, 2, 7);
hold on
plot(0:size(xx{4},2)-1,xx{4}(1,:))
plot(0:size(xx{4},2)-1,xx{4}(2,:))
ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 10);
legend('$x_1$','$x_2$','show', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 6);

subplot(4, 2, 8);
hold on
plot(0:size(uu{4},2)-1,uu{4}(1,:))
ylabel('$u$', 'Interpreter', 'latex', 'FontSize', 10);
%legend('$u$','Interpreter', 'latex', 'FontSize', 10);

%% Funtion
function y = vec2matrix(v,n,m)
    y=reshape(v,n,m);
end

function y = vecs_inv(v)
    s=length(v);
    n=0.5*(-1+sqrt(1+8*s));
    O=tril(ones(n));
    O(O==1)=v;
    y=0.5*(O'+O);
end

function y = vecv(n)
    %m la ma tran vuong
    %C1:
    m=n*n';
    p=[];
    s=size(m,1);
    for i=1:s
        for j=i:s
            p=[p;m(i,j)];
        end
    end
    y=p;
end