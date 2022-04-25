function dpsi = angle_error_pi(psi, psi0)

psi_2pi = mod(psi, 2*pi);
psi0_2pi = mod(psi0, 2*pi);

dpsi_all = (psi_2pi - psi0_2pi) + [- 2*pi, 0, 2*pi] ;
[minval,minind] = min(abs(dpsi_all));

dpsi = dpsi_all(minind);