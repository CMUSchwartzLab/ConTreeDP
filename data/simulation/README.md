configuration:

collapse: 
if freq_sum_simu[simulate_cp_tree[node]] - freq_sum_simu[node] <= 0.5/k:
    threshold = 0.9
elif freq_sum_simu[simulate_cp_tree[node]] - freq_sum_simu[node] <= 3/k:
    threshold = 0.5
else:
    threshold = 0.1

move:
prob_move = 0.1
