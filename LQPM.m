% left quaternion-product matrices (HAMILTON CONVENTION)
% q1 * q2 -> [q1]L * q2
% input : [qw qx qy qz]
% output : qL

function lqpm = LQPM(q)
    lqpm = [q(1) -q(2) -q(3) -q(4);
            q(2)  q(1) -q(4)  q(3);
            q(3)  q(4)  q(1) -q(2);
            q(4) -q(3)  q(2)  q(1)];
end