
dof_tbl = [1 2 5 10];

figure;
for i = 1:size(dof_tbl,2)
    v = randn(1e6, dof_tbl(i));
    X2 = sum(v.*v,2);
    subplot(2,2,i);   
    histogram(X2,1e3,'Normalization','pdf');
    title(dof_tbl(i));
end
    
tbl = [2 4 4 6 6 6];
edges = [1 3 5 7];
figure;
subplot(3,1,1);
histogram(tbl, 'BinWidth', 2);
subplot(3,1,2);
histogram(tbl,'Normalization','probability');
subplot(3,1,3);
histogram(tbl,'Normalization','pdf');