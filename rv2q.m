% left quaternion-product matrices (HAMILTON CONVENTION)
% q1 * q2 -> [q1]L * q2
% input : [qw qx qy qz]
% output : qL

function q = rv2q(rv)
    norm_rv = norm(rv);
    
    if(norm_rv < 1e-4)
        q = [1;0.5*rv];
    else
        q = [cos(norm_rv * 0.5); normalize(rv, 'norm') * sin(norm_rv * 0.5)];
    end
        
    q = normalize(q, 'norm');
    
end