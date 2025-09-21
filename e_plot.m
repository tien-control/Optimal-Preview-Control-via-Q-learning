t_end = 1000
%% Signals
% Disturbance
for i = 1:t_end
        if i >= 200 && i<= 400     d{i} = 0.01*td*(sin((1/120)*pi*i)+cos(i*pi/40));
        else                      d{i} = [0;0;0];
        end
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
                    ud = ud + K((k-1)*m+1:k*m,n+j*l+1:n+(j+1)*l)*d{i+j};
            else    ud = ud;
            end
        end
        u{i} = - (K((k-1)*m+1:k*m,1:n)*x{i} + ud);
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
plot(0:size(xx{1},2)-1,xx{1}(1:3,:))
%plot(0:size(xx{1},2)-1,xx{1}(2,:))
ylabel('$\omega$', 'Interpreter', 'latex', 'FontSize', 10);
legend('$\omega_1$','$\omega_2$','$\omega_3$','show', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 6);

subplot(4, 2, 2);
hold on
plot(0:size(uu{1},2)-1,uu{1}(:,:))
ylabel('$u$', 'Interpreter', 'latex', 'FontSize', 10);
%legend('$u_1$','$u_2$','$u_3$','$u_4$','$u_5$','$u_6$','show', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 6);

% Bo dieu khien 2
subplot(4, 2, 3);
hold on
plot(0:size(xx{2},2)-1,xx{2}(1:3,:))
%plot(0:size(xx{2},2)-1,xx{2}(2,:))
ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 10);
legend('$\omega_1$','$\omega_2$','$\omega_3$','show', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 6);

subplot(4, 2, 4);
hold on
plot(0:size(uu{2},2)-1,uu{2}(:,:))
ylabel('$u$', 'Interpreter', 'latex', 'FontSize', 10);
%legend('$u_1$','$u_2$','$u_3$','$u_4$','$u_5$','$u_6$','show', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 6);

% Bo dieu khien 3
subplot(4, 2, 5);
hold on
plot(0:size(xx{3},2)-1,xx{3}(1:3,:))
%plot(0:size(xx{3},2)-1,xx{3}(2,:))
ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 10);
legend('$\omega_1$','$\omega_2$','$\omega_3$','show', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 6);

subplot(4, 2, 6);
hold on
plot(0:size(uu{3},2)-1,uu{3}(:,:))
ylabel('$u$', 'Interpreter', 'latex', 'FontSize', 10);
%legend('$u_1$','$u_2$','$u_3$','$u_4$','$u_5$','$u_6$','show', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 6);


% Bo dieu khien cuoi
subplot(4, 2, 7);
hold on
plot(0:size(xx{4},2)-1,xx{4}(1:3,:))
%plot(0:size(xx{4},2)-1,xx{4}(2,:))
ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 10);
legend('$\omega_1$','$\omega_2$','$\omega_3$','show', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 6);

subplot(4, 2, 8);
hold on
plot(0:size(uu{4},2)-1,uu{4}(:,:))
ylabel('$u$', 'Interpreter', 'latex', 'FontSize', 10);
%legend('$u_1$','$u_2$','$u_3$','$u_4$','$u_5$','$u_6$','show', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 6);

figure(3)
% Bo dieu khien 1
subplot(4, 2, 1);
hold on
plot(0:size(xx{1},2)-1,xx{1}(4:6,:))
%plot(0:size(xx{1},2)-1,xx{1}(2,:))
ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 10);
legend('$\Omega_1$','$\Omega_2$','$\Omega_3$','show', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 6);

subplot(4, 2, 2);
hold on
plot(0:size(xx{1},2)-1,xx{1}(7:9,:))
%plot(0:size(xx{1},2)-1,xx{1}(2,:))
ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 10);
legend('$q_1$','$q_2$','$q_3$','show', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 6);

% Bo dieu khien 2
subplot(4, 2, 3);
hold on
plot(0:size(xx{2},2)-1,xx{2}(4:6,:))
%plot(0:size(xx{2},2)-1,xx{2}(2,:))
ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 10);
legend('$\Omega_1$','$\Omega_2$','$\Omega_3$','show', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 6);

subplot(4, 2, 4);
hold on
plot(0:size(xx{2},2)-1,xx{2}(7:9,:))
%plot(0:size(xx{2},2)-1,xx{2}(2,:))
ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 10);
legend('$q_1$','$q_2$','$q_3$','show', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 6);

% Bo dieu khien 3
subplot(4, 2, 5);
hold on
plot(0:size(xx{3},2)-1,xx{3}(4:6,:))
%plot(0:size(xx{3},2)-1,xx{3}(2,:))
ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 10);
legend('$\Omega_1$','$\Omega_2$','$\Omega_3$','show', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 6);

subplot(4, 2, 6);
hold on
plot(0:size(xx{3},2)-1,xx{3}(7:9,:))
%plot(0:size(xx{3},2)-1,xx{3}(2,:))
ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 10);
legend('$q_1$','$q_2$','$q_3$','show', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 6);


% Bo dieu khien cuoi
subplot(4, 2, 7);
hold on
plot(0:size(xx{4},2)-1,xx{4}(4:6,:))
%plot(0:size(xx{4},2)-1,xx{4}(2,:))
ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 10);
legend('$\Omega_1$','$\Omega_2$','$\Omega_3$','show', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 6);

subplot(4, 2, 8);
hold on
plot(0:size(xx{4},2)-1,xx{4}(7:9,:))
%plot(0:size(xx{4},2)-1,xx{4}(2,:))
ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 10);
legend('$q_1$','$q_2$','$q_3$','show', 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 6);