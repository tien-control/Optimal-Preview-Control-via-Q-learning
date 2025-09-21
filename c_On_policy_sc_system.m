%% On policy
clc; clear; close all;
%% Spacecraft
J = [250 150 100];
Jw = [100 100 100];
Q0 = diag([0.001, 0.001,0.001,0.001,0.001,0.001,0.02,0.02,0.02]);
R0 = diag([10^3,10^3,10^3,10^2,10^2,10^2]);

im = 57;
N = 5;
w0 = 0.0011;
ts = 58.6352;
td = [0.01; 0.01; 0.01];
muyf = 7.9*10^15;
r = 657000;

q_0 = [0.1,0.1,0.1];
w_0 = [0.00001,0.00001,0.00001]; % bien trang thai
Ohm_0 = [0.00001,0.00001,0.00001];
x_0 = [w_0 Ohm_0 q_0]';

A = cell(1,N);
B = cell(1,N);
D = cell(1,N);

for k = 1:N
    t = k*ts;

    A{k} = zeros(9,9);
    A{k}(1,3) = w0*(J(1)-J(2)+J(3))/(-J(1));
    A{k}(1,6) = w0*Jw(3)/(-J(1));
    A{k}(1,7) = 8*(w0^2)*(J(3)-J(2))/J(1);
    A{k}(2,8) = 6*(w0^2)*(J(3)-J(1))/J(2);
    A{k}(3,1) = w0*(J(1)-J(2)+J(3))/J(3);
    A{k}(3,4) = w0*Jw(1)/J(3);
    A{k}(3,9) = 2*(w0^2)*(J(1)-J(2))/J(3);
    A{k}(7,1) = 0.5;
    A{k}(8,2) = 0.5;
    A{k}(9,3) = 0.5;
    A{k} = (eye(9)+ts*A{k}); % roi rac he thong

    B{k} = zeros(9,6);
    B{k}(1,1) = J(1)^(-1);
    B{k}(1,5) = (2*muyf/(J(1)*r^3))*sin(im*pi/180)*sin(w0*t);
    B{k}(1,6) = (muyf/(J(1)*r^3))*cos(im*pi/180);
    B{k}(2,2) = J(2)^(-1);
    B{k}(2,4) = -(2*muyf/(J(2)*r^3))*sin(im*pi/180)*sin(w0*t);
    B{k}(2,6) = (muyf/(J(2)*r^3))*sin(im*pi/180)*cos(w0*t);
    B{k}(3,3) = J(3)^(-1);
    B{k}(3,4) = -(muyf/(J(3)*r^3))*cos(im*pi/180);
    B{k}(3,5) = -(muyf/(J(3)*r^3))*sin(im*pi/180)*cos(w0*t);
    B{k}(4,1) = -(Jw(1))^(-1);
    B{k}(5,2) = -(Jw(2))^(-1);
    B{k}(6,3) = -(Jw(3))^(-1);
    B{k} = ts*B{k}; % roi rac he thong

    D{k} = zeros(9,3);
    D{k}(1,1) = J(1)^(-1);
    D{k}(2,2) = J(2)^(-1);
    D{k}(3,3) = J(3)^(-1);
    D{k} = ts*D{k}; % roi rac he thong

    Q{k} = Q0;
    R{k} = R0;
end

%% Parameter
a = 3; t_end = 100;
n = size(A{1},2); m = size(B{1},2); l = size(D{1},2);
%% Signals
% Disturbance
for i = 1:t_end
        if i >= 300 && i<= 600      d{i} = td*sin((1/15)*pi*i);
        else                      d{i} = [0;0;0];
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
    Rs{i} = blkdiag(R{i});
end
    X_0 = [x_0; zeros((a+N+1)*l,1)];
    noise = [zeros(n,1); ones((a+N+1)*l,1)];
    Noise = [];
% Trans2
M0 = eye(n+(a+N+1)*l);
for i = 1:N
    M0 = Psi{i}*M0;
    M1 = eye(n+(a+N+1)*l); M2 = [];
    for j = 1:N
        if j < i
            M2 = [M2; zeros(n+(a+N+1)*l,m)];
        elseif j == i
            M2 = [M2; Lambda{j}];
        else
            M1 = Psi{j}*M1;
            M2 = [M2; M1*Lambda{i}];
        end
    end
    Psin = [Psin; M0]; Lambdan = [Lambdan, M2];
    Qsn = blkdiag(Qsn,Qs{i}); Rsn = blkdiag(Rsn,Rs{i});
    Noise = [Noise; noise];

end
Psin1 =  Psin(1:(N-1)*(n+(a+N+1)*l),:); Psin2 = Psin((N-1)*(n+(a+N+1)*l)+1:N*(n+(a+N+1)*l),:);
Psin = [zeros(N*(n+(a+N+1)*l),(N-1)*(n+(a+N+1)*l)),Psin];
Lambdan1 = Lambdan(1:(N-1)*(n+(a+N+1)*l),:); Lambdan2 = Lambdan((N-1)*(n+(a+N+1)*l)+1:N*(n+(a+N+1)*l),:);
Qsn1 = Qsn(1:(N-1)*(n+(a+N+1)*l),1:(N-1)*(n+(a+N+1)*l)); Qsn2 = Qsn((N-1)*(n+(a+N+1)*l)+1:N*(n+(a+N+1)*l),(N-1)*(n+(a+N+1)*l)+1:N*(n+(a+N+1)*l));
Rsn1= Rsn(1:(N-1)*m,1:(N-1)*m); Rsn2 = Rsn((N-1)*m+1:N*m,(N-1)*m+1:N*m);
    Xn_0 = [zeros((N-1)*(n+(a+N+1)*l),1); X_0];
%% On policy
% start controller
%[K,P,~] = dlqr(Psin,Lambdan,100*Qsn,100*Rsn);
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
[K,P,~] = dlqr(An,Bn,Qn,0.1*Rn);
K = 0.6*K(:,size(K,2)*(N-1)/N+1:end);
K = [K, zeros(m*N,(a+N+1)*l)];
%K = 0.7*K(:,size(K,2)*(N-1)/N+1:end);

%% Thuat toan Q_learning_on_policy
% Khoi tao bien
xM = Xn_0;
uuM = [];
xxM = [];

L = K;

L_s = [];
Norm_L(1)=norm(L);
for j = 1:12
    Z = []; Y = [];
    L_s(:,:,j) = L;
    L=[zeros(m*N,(n+(a+N+1)*l)*(N-1)),L]; %Create full L
 
    for i = 1:2500
        xold=xM;
        uM = -L*xold + 0.5*randn(N*m,1);
        uuM = [uuM uM];
        xM = Psin*xold + Lambdan*uM + 0.001*diag(randn((n+(a+N+1)*l)*N,1))*Noise;
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
    theta = inv(Z'*Z)*Z'*Y;
    H = vecs_inv(theta);
    H22 = H(end-m*N+1:end,end-m*N+1:end);
    H12 = H(1:(n+(a+N+1)*l),end-m*N+1:end);
    rank(H)
    L = inv(H22)*H12';

    Norm_L(j+1) = norm(L)
    Norm_H(j+1) = norm(H);
    
   % if (abs(Norm_L(j+1) - Norm_L(j)) <= 0.05)
    %    break;
    %end
end

L=[zeros(m*N,(n+(a+N+1)*l)*(N-1)),L];
for i=1:10
    xM=Psin*xM - Lambdan*L*xM;
    xxM=[xxM xM];
    uuM=[uuM -L*xM];
end


%% Ve do thi
% Qua trinh Q_learning
figure(1)
%subplot(2,1,1);
plot(0:length(Norm_L)-1,Norm_L,0:length(Norm_L)-1,Norm_L,'Color', [0 0.4470 0.7410],'LineStyle', '-', 'LineWidth', 1,'Marker', '*');
ylabel('$||K||$', 'FontSize', 10, 'Interpreter', 'latex');
%{
subplot(2,1,2);
plot(0:length(Norm_H)-1,Norm_H,0:length(Norm_H)-1,Norm_H,'Color', [1, 0, 0],'LineStyle', '-', 'LineWidth', 1,'Marker', '*');
ylabel('$||H||$', 'FontSize', 10, 'Interpreter', 'latex');
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
        ud = zeros(m,1);
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
%legend('$x_1$','$x_2$', 'Interpreter', 'latex', 'FontSize', 10);

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
%legend('$x_1$','$x_2$', 'Interpreter', 'latex', 'FontSize', 10);

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
%legend('$x_1$','$x_2$', 'Interpreter', 'latex', 'FontSize', 10);

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
%legend('$x_1$','$x_2$', 'Interpreter', 'latex', 'FontSize', 10);

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