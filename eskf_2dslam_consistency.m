%% filter configuration
closeloop = false;
FEJ = false;

%% transpose simulation signal
if(size(pos,1) > size(pos,2))
    pos = pos';
    velo = velo';
end

%% constants
N = 1;
landmark_meas_noise = 0.1;
velo_noise = 0.1;
w_noise = 0.1;
% 
% pos_init_error = [0.1;0.1];
% theta_init_error = 0.01;
% landmark_init_error = 0.1 * randn(2*N,1);

S = size(tout,1);
T = 0.01;
cross2 = [0 -1;
          1  0];
I2 = eye(2);

rng(1);

%% landmark set
landmarks = zeros(2*N,1);
for i = 1:N
    landmarks(2*i-1:2*i,1) = diag([60 30])*rand(2,1);
end

%% measurement generation
landmark_meas = zeros(2*N,S);
obsvalid = zeros(S,1);
for i = 1:S
        R_theta = theta2R(theta(i));
        k = 1;
        for j = 1:N
            if(abs(landmarks(2*j-1:2*j,1) - pos(:,i)) < 10)
                landmark_meas(2*k-1:2*k,i) = R_theta' * (landmarks(2*j-1:2*j,1) - pos(:,i)) + landmark_meas_noise * randn(2,1);
                landmark_meas(2*N+1) = landmark_meas(2*N+1) + 1;
                k = k +1;
            end
        end
        if(mod(i,10)==0)
            obsvalid(i) = 1;
        end
end
velo_meas = velo + velo_noise * randn(size(velo));
w_meas = w + w_noise * randn(size(w));


%% kalman filter initial state
x_kf = zeros(3+2*N, S);
cov_kf = zeros(3+2*N, 3+2*N, S);
pos_estimate = zeros(2, S);
theta_estimate = zeros(1, S);
landmark_estimate = zeros(2*N, S);
uncertanty_estimate = zeros( 3+2*N, S);

pos_estimate(:,1) = [0;0];
theta_estimate(1) = 0;
landmark_estimate(:,1) = landmark_meas(:,1);

x_kf(:,1) = [pos_estimate(:,1);
             theta_estimate(1);
             landmark_estimate(:,1)];
         
% cov_kf(1:2,1:2,1) = 0.01 * I2;
% cov_kf(3,3,1) = 0.01;
cov_kf(1:2,1:2,1) = 0 * I2;
cov_kf(3,3,1) = 0;

for i = 1:N
    cov_kf(2*i+2:2*i+3,2*i+2:2*i+3,1) = 0.01 * I2;
end

for i = 1:3+2*N
    uncertanty_estimate(i,1) = cov_kf(i,i,1);
end

if(FEJ)
    pos_pred_latch = pos_estimate(:,1);
    lm_FEJ = landmark_estimate(:,1);
end

%% kalman filter loop !

for k = 2:S

    %% rate value to delta value
    dpos = T * velo_meas(:,k);
    dtheta = T * w_meas(k);

    %% extract state vectors
    %   |--------------- x -------------|
    %   | pos(2) theta(1) landmark(2*N) |
    pos_k = x_kf(1:2,k-1);
    theta_k = x_kf(3,k-1);
    lm_k = x_kf(4:(3+2*N),k-1);
    R_theta_k = theta2R(theta_k);

    %% nominal prediction
    pos_k1_pred = pos_k + R_theta_k * dpos;
    theta_k1_pred = theta_k + dtheta;
    lm_k1_pred = lm_k;
    R_theta_k1_pred = theta2R(theta_k1_pred);

    %% covariance propagation
    if(FEJ)
        PHI_pos_theta = [I2 cross2 * (pos_k1_pred - pos_pred_latch);
                         0   0   1];
        pos_pred_latch = pos_k1_pred;
    else
        PHI_pos_theta = [I2 cross2 * R_theta_k * dpos;
                         0   0   1];
    end                 
    
    PHI_lm = eye(2*N);
    for i = 1:N
        PHI_lm(2*i-1:2*i, 2*i-1:2*i) = I2;
    end

    PHI = blkdiag(PHI_pos_theta, PHI_lm);

    Q_pos_theta = diag(T^2*[velo_noise^2 velo_noise^2 w_noise^2]);
    Q_lm = zeros(2*N, 2*N);
    Q = blkdiag(Q_pos_theta, Q_lm);

    cov_k1_pred = PHI * cov_kf(:,:,k-1) * PHI' + Q;

    if(closeloop && obsvalid(k))
        %% measurement processing
        % |----------- z -----------|
        % | pos_body_landmark (2*N) |
        R = landmark_meas_noise^2 * eye(2*N);

        predict_obs = zeros(2*N,1);
        for i = 1:N
            predict_obs(2*i-1:2*i,1) = R_theta_k1_pred' * (lm_k1_pred(2*i-1:2*i,1) - pos_k1_pred);
        end

        z = landmark_meas(:,k) - predict_obs;

        %% construct jacobian H
        H = zeros(2*N, 3+2*N);
        for i= 1:N % the ith measurement
            H_pos = - R_theta_k1_pred';
            if(FEJ)
                H_theta = - R_theta_k1_pred' * cross2 * (lm_FEJ(2*i-1:2*i) - pos_k1_pred);
            else    
                H_theta = - R_theta_k1_pred' * cross2 * (lm_k1_pred(2*i-1:2*i) - pos_k1_pred);
            end
            H_lm = R_theta_k1_pred';

            H(2*i-1:2*i , 1:2) = H_pos;
            H(2*i-1:2*i , 3) = H_theta;
            H(2*i-1:2*i , 2*i+2:2*i+3) = H_lm;
        end

        %% kalman filter 
        K = cov_k1_pred * H' * inv(H * cov_k1_pred * H' + R);
        delta_x = K * z;
        Ik = eye(size(K,1));
        cov_kf(:,:,k) = (Ik - K*H) * cov_k1_pred;   

        xk1_pred = [pos_k1_pred; 
               theta_k1_pred;
               lm_k1_pred];

        x_kf(:,k) = xk1_pred + delta_x;
    else
        %% only prediction and propagation
         x_kf(:,k) = [pos_k1_pred; 
               theta_k1_pred;
               lm_k1_pred];
        cov_kf(:,:,k) = cov_k1_pred;
    end
    
    pos_estimate(:,k) = x_kf(1:2,k);
    theta_estimate(k) = x_kf(3,k);
    landmark_estimate(:,k) = x_kf(4:(3+2*N),k);    
    
    for i = 1:3+2*N
        uncertanty_estimate(i,k) = cov_kf(i,i,k);
    end

end

angle_error = zeros(size(theta));
for i = 1:S
    angle_error(i) = angle_error_pi(theta(i), theta_estimate(i));
end

figure;
f311 = subplot(311);
plot(tout, pos(1,:), 'r', tout, pos_estimate(1,:), 'g');title('x pos(m)');legend('gt','estimate');grid on;
f312 = subplot(312);
plot(tout, pos(2,:), 'r', tout, pos_estimate(2,:), 'g');title('y pos(m)');legend('gt','estimate');grid on;
f313 = subplot(313);
plot(tout, theta, 'r', tout, theta_estimate, 'g');title('theta(rad)');legend('gt','estimate');grid on;
linkaxes([f311,f312,f313],'x');

figure;
g311 = subplot(311);
plot(tout, pos(1,:)-pos_estimate(1,:), 'r', tout, 3*sqrt(uncertanty_estimate(1,:)), 'g-.', tout, -3*sqrt(uncertanty_estimate(1,:)), 'g-.');title('x pos error/uncertainty(m)');legend('error','+3 sigma','-3 sigma');grid on;
g312 = subplot(312);
plot(tout, pos(2,:)-pos_estimate(2,:), 'r', tout, 3*sqrt(uncertanty_estimate(2,:)), 'g-.', tout, -3*sqrt(uncertanty_estimate(2,:)), 'g-.');title('y pos error/uncertainty(m)');legend('error','+3 sigma','-3 sigma');grid on;
g313 = subplot(313);
plot(tout, angle_error, 'r', tout, 3*sqrt(uncertanty_estimate(3,:)), 'g-.', tout, -3*sqrt(uncertanty_estimate(3,:)), 'g-.');title('theta error/uncertainty(rad)');legend('error','+3 sigma','-3 sigma');grid on;
linkaxes([g311,g312,g313],'x');

figure;
plot(pos(1,:), pos(2,:), 'r-.', pos_estimate(1,:), pos_estimate(2,:), 'g-.');title('estimate visualization');xlabel('x(m)');ylabel('y(m)');
hold on;
plot(landmarks(1:2:end), landmarks(2:2:end), 'ro', landmark_estimate(1:2:end,end), landmark_estimate(2:2:end,end), 'b+');legend('pos gt','pos estimate','ldmk gt','ldmk estimate');
for i = 1:N
    hold on;
    plot(landmark_estimate(2*i-1,:), landmark_estimate(2*i,:), 'cyan');
end
axis equal;