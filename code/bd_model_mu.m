function [ mu] = gen_model_mu(p,tsamp, Ninit, Nsim)
bdmodel_out= @(p)run_bdmodel(p, tsamp, Ninit, Nsim);
[~, ~, Cmoms, ~, ~, ~] = bdmodel_out(p);
mu_exp = Cmoms(:,1,:);
mu_exp_long = reshape(mu_exp,size(mu_exp,1)*size(mu_exp,2)*size(mu_exp,3),1);
mu = mu_exp_long;
end