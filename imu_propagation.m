addpath(genpath('/home/zhangtao/workings/projects/psins'));

S = size(tout,1);
T = 0.001;
I3 = eye(3);
Z3 = zeros(3,3);
g = [ 0;0;-9.81];

if(size(pos,2) < size(pos,1))
    pos = pos.';
    vel = vel.';
    eul = eul.';
    acc = acc.';
    gyr = gyr.';
end


















gyro_noise = 0.003;
acc_noise = 0.03;

acc_meas = acc + acc_noise * randn(size(S));
gyro_meas = w + gyro_noise * randn(size(S));

x_kf = zeros(15, S);
cov = zeros(15, 15, S);

pos_estimate = zeros(3,S);
vel_estimate = zeros(3,S);
rot_estimate = zeros(3,3,S);
ba_estimate = zeros(3,S);
bg_estimate = zeros(3,S);

rot_estimate(:,:,1) = I3;

Q = diag([0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001]);

pos_var = zeros(3,S);
vel_var = zeros(3,S);
rot_var = zeros(3,S);
ba_var = zeros(3,S);
bg_var = zeros(3,S);


for i = 2:S
    posk = pos_estimate(:,i-1);
    velk = vel_estimate(:,i-1);
    rotk = rot_estimate(:,:,i-1);
    bak = ba_estimate(:,i-1);
    bgk = bg_estimate(:,i-1);
    
    acck = acc_meas(:,i-1) - bak;
    gyrok = gyro_meas(:,i-1) - bgk;
    
    pos_pred = posk + velk * T+ 0.5 * (acck - rotk.' * g) * T^2;
    vel_pred = velk + (acck - rotk.' * g) * T;
    rot_pred = rotk * hatmap(gyrok * T);
    ba_pred = bak;
    bg_pred = bgk;
    
    pos_estimate(:,i) = pos_pred;
    vel_estimate(:,i) = vel_pred;
    rot_estimate(:,:,i) = rot_pred;
    ba_estimate(:,i) = ba_pred;
    bg_estimate(:,i) = bg_pred;
    
    PHI = [I3 I3*T Z3 Z3 Z3 ; Z3 I3 -hatmap(rot_pred*acck)*T -rot_pred*T Z3; Z3 Z3 I3 Z3 -rot_pred*T; Z3 Z3 Z3 I3 Z3; Z3 Z3 Z3 Z3 I3];
    
    P_pred = PHI * cov(:,:,i-1) * PHI.' + Q * T;
    
    cov(:,:,i) = P_pred;

    pos_var(:,i) = [cov(1,1,i);cov(2,2,i);cov(3,3,i)];
    vel_var(:,i) = [cov(4,4,i);cov(5,5,i);cov(6,6,i)];
    rot_var(:,i) = [cov(7,7,i);cov(8,8,i);cov(9,9,i)];
    ba_var(:,i) = [cov(10,10,i);cov(11,11,i);cov(12,12,i)];
    bg_var(:,i) = [cov(13,13,i);cov(14,14,i);cov(15,15,i)];
end

figure;
ax_pos_dev = subplot(3,1,1);
plot(tout, sqrt(pos_var(1,:)), tout, sqrt(pos_var(2,:)), tout, sqrt(pos_var(3,:)));xlabel('time(s)');title('pos estimation deviation');ylabel('m');grid on;legend('x','y','z');
ax_vel_dev = subplot(3,1,2);
plot(tout, sqrt(vel_var(1,:)), tout, sqrt(vel_var(2,:)), tout, sqrt(vel_var(3,:)));xlabel('time(s)');title('velo estimation deviation');ylabel('m/s');grid on;legend('x','y','z');
ax_euler_dev = subplot(3,1,3);
plot(tout, sqrt(rot_var(1,:)), tout, sqrt(rot_var(2,:)), tout, sqrt(rot_var(3,:)));xlabel('time(s)');title('euler estimation deviation');ylabel('rad');grid on;legend('roll','pitch','yaw');
linkaxes([ax_pos_dev ax_vel_dev ax_euler_dev],'x')

figure;
ax_ba_dev = subplot(2,1,1);
plot(tout, sqrt(ba_var(1,:)), tout, sqrt(ba_var(2,:)), tout, sqrt(ba_var(3,:)));xlabel('time(s)');title('ba estimation deviation');ylabel('m/ss');grid on;legend('x','y','z');
ax_bg_dev = subplot(2,1,2);
plot(tout, sqrt(bg_var(1,:)), tout, sqrt(bg_var(2,:)), tout, sqrt(bg_var(3,:)));xlabel('time(s)');title('bg estimation deviation');ylabel('rad/s');grid on;legend('x','y','z');
linkaxes([ax_ba_dev ax_bg_dev],'x')

figure;
ax_pos_dev = subplot(3,1,1);
plot(tout, pos_var(1,:), tout, pos_var(2,:), tout, pos_var(3,:));xlabel('time(s)');title('pos estimation var');ylabel('m');grid on;legend('x','y','z');
ax_vel_dev = subplot(3,1,2);
plot(tout, vel_var(1,:), tout, vel_var(2,:), tout, vel_var(3,:));xlabel('time(s)');title('velo estimation var');ylabel('m/s');grid on;legend('x','y','z');
ax_euler_dev = subplot(3,1,3);
plot(tout, rot_var(1,:), tout, rot_var(2,:), tout, rot_var(3,:));xlabel('time(s)');title('euler estimation var');ylabel('rad');grid on;legend('roll','pitch','yaw');
linkaxes([ax_pos_dev ax_vel_dev ax_euler_dev],'x')

figure;
ax_ba_dev = subplot(2,1,1);
plot(tout, ba_var(1,:), tout, ba_var(2,:), tout, ba_var(3,:));xlabel('time(s)');title('ba estimation var');ylabel('m/ss');grid on;legend('x','y','z');
ax_bg_dev = subplot(2,1,2);
plot(tout, bg_var(1,:), tout, bg_var(2,:), tout, bg_var(3,:));xlabel('time(s)');title('bg estimation var');ylabel('rad/s');grid on;legend('x','y','z');
linkaxes([ax_ba_dev ax_bg_dev],'x')    
    
    
    
    
