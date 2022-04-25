if(size(gyr,1) > size(gyr,2))
    gyr = gyr.';
end

rng(1);

gyro_noise_scale = 0.001;
gyro_noise = randn(size(gyr,1),size(gyr,2)) * gyro_noise_scale;
gyr = gyr + gyro_noise;

len = size(gyr,2);
T = 0.01;
ts = (0:1:len-1)*T;

%% trival mid-point integration
Rwb_estimate = Rbw;
error_angle = zeros(3,len);
q_estimate = zeros(4,len);
q_estimate(:,1) = [1;0;0;0];

intmode = 2; % 0-euler 1-midpoint 2-rk4


for i = 2:len
    if(intmode == 0)
        dphi = gyr(:,i) * T;
        dq = rv2q(dphi);
        q_estimate(:,i) = LQPM(q_estimate(:,i-1)) * dq;
        normalize(q_estimate(:,i),'norm');    
    elseif(intmode == 1)
        dphi = (gyr(:,i) + gyr(:,i-1)) * 0.5 * T;            
        dq = rv2q(dphi);
        q_estimate(:,i) = LQPM(q_estimate(:,i-1)) * dq;
        normalize(q_estimate(:,i),'norm');  
    else
        
        %k1 = 0.5 * q0 ** w0;        
        k1 = 0.5 * LQPM(q_estimate(:,i-1)) * [0;gyr(:,i-1)];
        %k2 = 0.5 * q0 ** (k1*0.5T) ** (w0+w1)/2
        k2 = 0.5 * LQPM((q_estimate(:,i-1)+k1*0.5*T)) * [0;(gyr(:,i-1)+gyr(:,i))*0.5];
        %k3 = 0.5 * q0 ** (k2*0.5T) ** (w0+w1)/2
        k3 = 0.5 * LQPM((q_estimate(:,i-1)+k2*0.5*T)) * [0;(gyr(:,i-1)+gyr(:,i))*0.5];
        %k4 = 0.5 * q0 ** (k3*T) ** w1
        k4 = 0.5 * LQPM((q_estimate(:,i-1)+k3*T)) * [0;gyr(:,i)];
        %q1 = q0 ** T/6 * (k1 + 2k2 + 2k3 + k4)      

        q_estimate(:,i) = q_estimate(:,i-1) + (k1+2*k2+2*k3+k4) * T /6;
        normalize(q_estimate(:,i),'norm');  
    end
        
        
    %dphi = (gyr(:,i) + gyr(:,i-1)) * 0.5 * T + skew_sym(gyr(:,i-1)) * gyr(:,i) * T^2 /12 ;
    

    
    Rwb_estimate(:,:,i) = q2R(q_estimate(:,i));
    
    error_angle(:,i) = reverse_skew_sym(Rbw(:,:,i) * Rwb_estimate(:,:,i));
end

figure;
plot(ts, error_angle(1,:), 'r', ts, error_angle(2,:), 'g', ts, error_angle(3,:), 'b');
legend('x', 'y', 'z');


    
    