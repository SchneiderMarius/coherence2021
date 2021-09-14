function WilsonCowan(W, input, alpha, beta,init_state, par,response_fn,folder)

[spike_times,spike_ids] = ...
    full_sim_batch(response_fn, W, input, alpha, beta, ...
    init_state, par.total_secs);
SPt = zeros(sum(par.cellnum),par.total_secs*1000);

for cnt1 = 1 : sum(par.cellnum) 
    sp = ceil(spike_times(spike_ids==cnt1));
    SPt(cnt1,sp(sp<=3000)) = 1;
end
SPt = sparse(SPt);

save(fullfile(folder),'SPt','par')
end