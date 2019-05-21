%% MCMC Search algorithm
% 1. Uniformly sample domain and calculate negLL at each point
% 2. Take lowest negLL as initial guess 
% 3. Random walk to search for better neg LL

% set domain of b, d, and A
bvec = [ 0: 0.005: 0.05];
dvec = [ 0: 0.001: 0.01];

% make a 3D mesh
[B,D] = meshgrid(bvec, dvec); % 3D grid of parameters
Bflat = reshape(B,1,[]); % vector of bs
Dflat = reshape(D,1,[]); % vector of ds

% run loop through vector of parameters and calculate negLL

    for i = 1:length(Bflat)
        % set params
        pguess(i,1:2) = horzcat(Bflat(i),Dflat(i));
        negLL(i)= objfun_J(pfxform(pguess(i,1:2)));
        pguess(i,3) = negLL(i);   
    end

NEGLL = reshape(negLL, size(B));
[Jinit,imin] = min(pguess(:,3));
pinit = pguess(imin,1:2)

figure;
hold off;
surf(bvec,dvec,NEGLL(:,:));
hold on
plot3(pinit(1), pinit(2), Jinit, 'r*', 'LineWidth',8)
xlabel( 'b')
ylabel('d')
zlabel('negLL')
title('Initial parameter space search')
%%
J_curr = Jinit; % set current negLL value to inital value to initialize
count_accepted = 0; % start with no changes accepted

% Initialize a vector to store b and d params in
nruns = 10000;
store_params = zeros(nruns, 3); % store J, b,  d,for each run

% Outside the loop: initialize temperature and iterations
T0= 20;
k = 1:1:nruns;
T= T0*exp(-4*k./(nruns)); % sets up gradual cooling


% check that temperature annealing is working
figure;
hold off
set(gca,'LineWidth',1.5,'FontSize',12);
plot(k, T, 'LineWidth', 3)
xlabel('Markov step')
ylabel('Temperature')
title('Temperature Annealing')

bstep = 0.001; % step size for searching birth param
dstep = 0.001; % step size for searchign death param

% set initial guess as b1 and d1
b1=pinit(1);
d1=pinit(2);


store_acc_params = [];
store_acc_params(1,1) = Jinit;
store_acc_params(1,2) = b1;
store_acc_params(1,3) = d1;
theta1 = pinit;
    for k = 1:nruns
        r = rand;
        % set up to only change one parameter at a time, use rand to choose the
        % parameter to update
        if r<0.5
            % update b
            b2 = b1 + bstep*(2*rand-1);
            d2 = d1;
        end
        if r>0.5
            % update d
            d2 = d1 + dstep*(2*rand-1);
            b2 = b1;
        end


    % Constrain search region to domain
        if b2<0 
            b2=0;
        end
        if b2>.05
            b2 = .05;
        end
        if d2<0
            d2=0;
        end
        if d2>.01
            d2 = .01;
        end


        % find the neg LL of the new params
        theta2 = horzcat(b2,d2);
        J_new = objfun_J( pfxform(theta2));
        % store the  negLL and searched parameters
        store_params(k,1) = J_new;
        store_params(k,2) =b2;
        store_params(k,3) = d2;


        prob(k) = exp((J_curr-J_new)./T(k));
        % if Jcurr > Jnew, Jnew is better--the numerator in the exponent is positive &
        % prob>1--> change will always be accepted
        % if Jcurr < Jnew, numerator is negative, & prob <1 --> chance change
        % will be accepted

        if rand < prob(k) % if true, accept the change!
            b1 = b2;
            d1 = d2;

            theta1 = theta2;
            J_curr = J_new;
            count_accepted = count_accepted +1;
            % decrease search step size
            if r<0.5
                % update b step
                bstep = 0.999*bstep;
            end
            if r>0.5
                % update d step
                 dstep = 0.999*dstep;
            end

            params(1,1) = J_curr;
            params(1,2)= b1;
            params(1,3)= d1;
            store_acc_params= vertcat(store_acc_params, params);
        else
            % increase search step size
            if r<0.5
                % update b step
                bstep = 1.001*bstep;
            end
            if r>0.5
                % update d step
                 dstep = 1.001*dstep;
            end
        end
      end
%%
params_best_MC = store_acc_params(end,2:3)
[lowest_LL,ind] = min(store_acc_params(:,1));
lowestLLparams = store_acc_params(ind, 2:3)
negLLMC = lowest_LL;

likelihoodval = exp(-store_acc_params(:,1));
likelihoodvalall = exp(-store_params(:,1));
pointsize = 30;
figure;
scatter(store_acc_params(:,2), store_acc_params(:,3), pointsize, likelihoodval,'filled');
colorbar
xlabel('b')
ylabel('d')
title(' b versus d colored by likelihood for accepted parameters')

