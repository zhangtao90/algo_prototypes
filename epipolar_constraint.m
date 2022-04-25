%% data generation
% data1 2 cams:
% Rw1 = eye(3);
% tw1 = zeros(3,1);
% 
% Rw2 = eye(3);
% tw2 = [1;0;0];

N = 30;

Rw1 = [  0.9810603, -0.1804681, -0.0703706;
   0.1729874,  0.9797412, -0.1009085;
   0.0871557,  0.0868241,  0.9924039 ];
tw1 = [-2;1;1];

Rw2 = [  0.9305476,  0.3372476,  0.1426366;
  -0.3386916,  0.9407814, -0.0147762;
  -0.1391731, -0.0345599,  0.9896649 ];
tw2 = [8;2;1];

R12 = Rw1.' * Rw2;
t12 = Rw1.' * (tw2 - tw1);

R1w = Rw1.';
t1w = - R1w * tw1;
R2w = Rw2.';
t2w = - R2w * tw2;

ldmks = rand(3,N);
ldmks = ldmks + [0;0;10];
ldmks = diag([6,8,1]) * ldmks;

figure;
plot3(ldmks(1,:),ldmks(2,:),ldmks(3,:),'b.');
hold on;plot3(tw1(1),tw1(2),tw1(3),'r*');
hold on;plot3(tw2(1),tw2(2),tw2(3),'g*');
axis equal;
grid on;
xlabel('x');ylabel('y');zlabel('z');

pc1 = R1w * ldmks + t1w;
obs1 = pc1(1:2,:)./pc1(3,:);

pc2 = R2w * ldmks + t2w;
obs2 = pc2(1:2,:)./pc2(3,:);

%% essential matrix solving
A = zeros(8,9);
for i = 1:8
    u1 = obs1(1,i);
    v1 = obs1(2,i);
    u2 = obs2(1,i);
    v2 = obs2(2,i);
    
    A(i,:) = [u1*u2, u1*v2, u1, v1*u2, v1*v2, v1, u2, v2, 1];
end

[uA,sA,vA] = svd(A);

em = vA(:,9);

E12 = [em(1), em(2), em(3); 
     em(4), em(5), em(6); 
     em(7), em(8), em(9)];
 
[uE,sE,vE] = svd(E12);

if(det(uE)<0)
    uE = -uE;
end
if(det(vE)<0)
    vE = -vE;
end

W = [0,-1,0; 1,0,0; 0,0,1];
R2 = uE * W * vE.';
t2 = uE(:,3);

M = skew_sym(R2.'*t2);
X1 = M * [obs1(:,1);1];
X2 = M * R2.' * [obs2(:,1);1];

if(X1(3) * X2(3) < 0)
   R2 = uE * W' * vE.'; 
   M = skew_sym(R2.'*t2);
   X1 = M * [obs1(:,1);1];
end

if(X1(3) < 0)
    t2 = -t2;
end




 
% Rz_halfpi_pos = Exp([0 0 pi/2]);
% Rz_halfpi_neg = Exp([0 0 -pi/2]);
% 
% t_candi1_ss = uE * Rz_halfpi_pos * sE * uE.';
% t_candi2_ss = uE * Rz_halfpi_neg * sE * uE.';
% 
% t_candi1 = reverse_skew_sym(t_candi1_ss);
% t_candi2 = reverse_skew_sym(t_candi2_ss);
% 
% R_candi1 = uE * Rz_halfpi_pos.' * vE.';
% R_candi2 = uE * Rz_halfpi_neg.' * vE.';




