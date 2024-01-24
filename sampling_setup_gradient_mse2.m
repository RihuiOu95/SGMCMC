Ls = [2, 12];
%L = 2;
%Bs = [5, 10];
B = 5;
Ts = [1e4, 5e4, 1e5];
STDs = [0, 1, 2, 3];
n1 = length(Ls);
n2 = length(Ts); 
n3 = length(STDs);
cell_clu = cell(n1,n2,n3);
cell_unif = cell(n1,n2,n3);
for i = 1:n1
    for j = 1:n2
        for k = 1:n3
            L = Ls(i);
            T = Ts(j);
            STD = STDs(k);
            [res1, res2] = sampling_setup_gradient_mse(T, L, B, STD);
            cell1{i, j, k} = res1;
            cell2{i, j, k} = res2;
        end
    end
end

            