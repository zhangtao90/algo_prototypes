% left quaternion-product matrices (HAMILTON CONVENTION)
% q1 * q2 -> [q1]L * q2
% input : [qw qx qy qz]
% output : qL

function R = q2R(q)
    qw = q(1);
    qv = [q(2);q(3);q(4)];
    
    R = (qw^2 - qv.' * qv) * eye(3) + 2 * qv * qv.' + 2 * qw * skew_sym(qv);
end