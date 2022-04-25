
%% data generation
% data1 2 cams:
Rw1 = eye(3);
tw1 = zeros(3,1);

Rw2 = eye(3);
tw2 = [1;0;0];

Rw1 = [  0.9810603, -0.1804681, -0.0703706;
   0.1729874,  0.9797412, -0.1009085;
   0.0871557,  0.0868241,  0.9924039 ];
tw1 = [-2;1;1];

Rw2 = [  0.9305476,  0.3372476,  0.1426366;
  -0.3386916,  0.9407814, -0.0147762;
  -0.1391731, -0.0345599,  0.9896649 ];
tw2 = [8;3;13];

p = [0;0;30];

p1 = Rw1.' *( p - tw1);
p2 = Rw2.' *( p - tw2);

obs_c1 = obs(Rw1,tw1,p);
obs_c2 = obs(Rw2,tw2,p);

noise_c1 = 0.001 * randn(2,1);
noise_c2 = 0.001 * randn(2,1);

obs_c1 = obs_c1 + noise_c1;
obs_c2 = obs_c2 + noise_c2;


% naive two obs least square method using xyz

R12 = Rw1.' * Rw2;
t12 = Rw1.' * (tw2 - tw1);

m = R12 * [obs_c2(1);obs_c2(2);1];

A = [[obs_c1;1] -m];
b = t12;

x_lsq = inv(A.'*A)* A.' * b;

err1 = p1(3) - x_lsq(1);
err2 = p2(3) - x_lsq(2);


% inverse depth method using alpha beta rho


function uv = obs(Rwc,twc,p)
    Rcw = Rwc.';
    tcw = - Rcw * twc;
    
    pc = Rcw * p + tcw;
    
    %% pc = Rwc.'*(p - twc);
    
    uv = [pc(1)/pc(3);pc(2)/pc(3)];
end