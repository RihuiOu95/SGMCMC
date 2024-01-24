resha_z = reshape(z,[5,347865/5])';
uni_z = unique(resha_z,'rows');
prob_vec = zeros([size(uni_z,1),1]);
rare_state = 4;
weights = weights_alex.trans{rare_state,rare_state};
for i = 1:size(uni_z,1)
    target = uni_z(i,:);
    dif_mat = sum(abs(resha_z - target),2);
    index = find(dif_mat<1e-6,1);
    prob_vec(i) = weights(index);
end
prob_mat = [uni_z ,prob_vec];
prob_tab = sortrows(prob_mat,6,'descend');
T = array2table(prob_tab);
T.Properties.VariableNames(1:6) = {'z_first','z_second','z_third','z_forth','z_fifth','weights'};