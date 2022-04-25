function R = Exp(rv)

phi = norm(rv);
u = normalize(rv, 'norm');

R = eye(3) + sin(phi) * skew_sym(u) + (1-cos(phi)) * skew_sym(u) * skew_sym(u);