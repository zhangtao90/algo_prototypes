dmin = 0.03;



iz = dmin * randn(100000,1);

izp = [];

for i = 1:size(iz,1)
    if(iz(i)>0)
        izp = [izp iz(i)];
    end
end

figure;
histogram(izp, 'Normalization', 'pdf');

z = 1./izp;

figure;
histogram(z,[-1:0.1:300], 'Normalization', 'pdf');