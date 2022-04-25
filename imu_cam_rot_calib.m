%% data generation
Ric = [  0.9810603, -0.1804681, -0.0703706;
         0.1729874,  0.9797412, -0.1009085;
         0.0871557,  0.0868241,  0.9924039 ];

MIS = [1 0.01 -0.01;
       0.02 1 -0.12;
       -0.003 0.013 1];
   
SCL = diag([0.998 0.979 1.022]);

BIA = [-0.02; 0.032;-0.019];

noise  = 0;
N = 1000;

wc_sample = randn(3,N);
figure;
plot3(wc_sample(1,:),wc_sample(2,:),wc_sample(3,:),'.');

%wc_sample = normalize(wc_sample,1, 'norm');

figure;
plot3(wc_sample(1,:),wc_sample(2,:),wc_sample(3,:),'.');

wi_sample = Ric * err_model(wc_sample, MIS, SCL, BIA) + noise * randn(3,N);

%%
syms k11 k12 k13 k21 k22 k23 k31 k32 k33 b1 b2 b3 w1 w2 w3
K = [k11 k12 k13;
     k21 k22 k23;
     k31 k32 k33];

Kw = K * [w1;w2;w3] + [b1;b2;b3];


%% solve

A = zeros(3*N, 12);
b = zeros(3*N, 1);

z13 = zeros(1,3);
I3 = eye(3);

for i = 1:N
    w123 = [wc_sample(1,i) wc_sample(2,i) wc_sample(3,i)];
    A(i*3+1:i*3+3, :) = [[w123   z13   z13  ;
                          z13   w123   z13  ;
                          z13   z13    w123 ] I3];
    b(i*3+1:i*3+3) = [wi_sample(1,i);wi_sample(2,i);wi_sample(3,i)];
end

x = A\b;

Kic = [x(1) x(2) x(3);
       x(4) x(5) x(6);
       x(7) x(8) x(9)];
bic = [x(10);x(11);x(12)];


dK = Kic - Ric * MIS * SCL;
db = bic - Ric * MIS * SCL * BIA;

norm(dK)/norm(Kic)
norm(db)/norm(bic)


%% forward error model
function w_meas = err_model(w_true, MIS, SCL, BIA)
    w_meas = MIS * SCL * w_true + BIA;
end
    