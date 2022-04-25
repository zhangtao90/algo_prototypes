psi = 4*pi*2*(rand(10000,1)-0.5);
psi0 = 4*pi*2*(rand(10000,1)-0.5);

raw_dif = psi - psi0;

proc_dif = zeros(size(psi));
for i = 1:size(psi,1)
    proc_dif(i) = angle_error_pi(psi(i), psi0(i));
end


%% test if -pi <-> pi limit ok
off_limit_cnt = uint16(0);
for i = 1:size(proc_dif,1)
    if(proc_dif(i)<-pi || proc_dif(i) > pi)
        off_limit_cnt = off_limit_cnt + uint16(1);
    end
end

disp('off limit count is :');
disp(off_limit_cnt);

%% test if dif is correct
dif_wrong_cnt = uint16(0);
difdif = abs(mod(raw_dif - proc_dif, 2*pi));

for i = 1:size(difdif,1)
    if(abs(difdif(i))> 1e-5 && abs(difdif(i)-2*pi)>1e-5)
        dif_wrong_cnt = dif_wrong_cnt + uint16(1);
    end
end

disp('dif_wrong_cnt is :');
disp(dif_wrong_cnt);