% Gibbs sampler for Dirichlet Process Mixtures of Gaussians
% The base distribution G0 is Normal inverse Wishart of parameters given by
% hyperG0
% Reference: Algorithm 1 of Neal

function c_st = gibbsDPM_algo1(P, hyperG0, alpha, niter, doPlot)

    if doPlot
        figure('name','Gibbs sampling for DPM');
        colormap('default');
        cmap = colormap;
    end

    [p, n] = size(P);
    c_st = zeros(n, niter/2);

    if (hyperG0.prior == 'NIW')
        theta_Sigma = zeros(p, p, n);
    elseif (hyperG0.prior == 'NIG')
        p=p-1; % the y-coordinate doesn't count as dimension
        theta_Sigma = zeros(n);
    end
    theta_mu = zeros(p, n);

    % Initialisation
    hyper = update_SS(P(:, 1), hyperG0);
    %[theta_mu(:, 1), theta_Sigma(:,:, 1)] = normalinvwishrnd(hyper);
    R = sample_pdf(hyper);
    theta_mu(:, 1) = R.mu;
    if (hyperG0.prior == 'NIW')
        theta_Sigma(:, :, 1) = R.Sigma;
    else
        theta_Sigma(1) = R.Sigma;
    end

    for k=2:n
        if (hyperG0.prior == 'NIW')
            [theta_mu(:, k), theta_Sigma(:, :, k)] = sample_theta(alpha,...
                P(:,k), hyperG0, theta_mu(:, 1:k-1), theta_Sigma(:, :, 1:k-1));
        else
            [theta_mu(:, k), theta_Sigma(k)] = sample_theta(alpha,...
                P(:,k), hyperG0, theta_mu(:, 1:k-1), theta_Sigma(1:k-1));
        end
    end

    % Iterations
    for i=2:niter
        for k=1:n
            % Sample theta_k | theta_{-k}, y_k
            ind_notk = (1:n)~=k;
            if (hyperG0.prior == 'NIW')
                [theta_mu(:, k), theta_Sigma(:, :, k)] = sample_theta(alpha,...
                    P(:,k), hyperG0, theta_mu(:, ind_notk), theta_Sigma(:, :, ind_notk));
            else
                [theta_mu(:, k), theta_Sigma(k)] = sample_theta(alpha,...
                    P(:,k), hyperG0, theta_mu(:, ind_notk), theta_Sigma(ind_notk));
            end
            if doPlot==1
                some_plot(y, theta_mu, theta_Sigma, k, i, n, cmap);
            end
        end
        % get the partition
        [~, ~, c] = unique(theta_mu(1, :));

        print_clusters=true;
        if (print_clusters)
            fprintf('Iteration %d/%d\n', i, niter);
            fprintf('%d clusters\n\n', length(unique(c)));
        end

        if doPlot==2
            some_plot(P(:,:), theta_mu, theta_Sigma, k, i, n, cmap);
        end

        if i>niter/2
           c_st(:, i-niter/2) = c;
        end
    end
end

%%%% Subfunctions

% each parameter theta has two fields, mu and Sigma
% we want to sample theta_k | theta_not_k, y_k
% we sample it from sum_j_not_k q_kj delta(theta_j) + r_k H_k
% here is H_k the posterior distribution for theta based on G_0 and the observation y_k with likelihood F(y_k,theta)
% q_kj = b F(y_k,theta_j)
% r_k = b alpha integral F(y_k,theta) dG_0(theta)
% b = sum of all j not k over q_kj + r_k = 1

% z is y_k
% n0 is H_k, the posterior for theta based on G0 and y_k, or is n0 r_k?
% n is sum j not k F(y_k,theta)

% Escobar calls the dispersion factor A_0 instead of alpha
%         writes the normal as phi(Y_i - X_i), which is often more generally noted as the likelihood F(y_i,theta_i)
%

% TODO: replace theta_mu_notk and theta_Sigma_notk by variable theta.mu_notk and theta.Sigma_notk
function [theta_mu_k, theta_Sigma_k] = sample_theta(alpha, z, hyperG0, theta_mu_notk, theta_Sigma_notk)
    % number of items (grows to number of data points with initialization and stays at N-1)
    p = size(theta_mu_notk, 2);

    % q_ij, calculate b F(y_i,theta_j) for all items in the sum (except j)
    % equation 3.3 in Neal for all j unequal to i (sum in eq. 3.2)
    % n = exp(loggausspdf(repmat(z, 1, p)', theta_mu_notk', theta_Sigma_notk));
    % I transposed sigma here!!
    n = likelihoods(hyperG0, repmat(z, 1, p), theta_mu_notk', theta_Sigma_notk');
    % prediction is A(Y_i) = A_0 integral phi(Y_i-X) G_0(dX) with Y_i the new observation z
    n0 = pred(z, hyperG0);
    % the denominator, A(Y_i) + sum_l phi(Y_i - X_l) with A(Y_i) = A_0 integral phi(Y_i-X) G_0(dX)
    const = alpha*n0 +sum(n);

    % (sum(n) + alpha*n0)/const = 1

    % we only need to calculate bottom equation in 3
    p0 = alpha*n0 / const; % p0 is the r_i H_i part
    %fprintf('1-p0=%i\n', 1-p0);
    % pick uniform random number between 0 and 1
    u = rand;
    if u<p0
	% Sample new value
	% H_i, h(X_i|Y_i) = A_0/A(Y_i) phi(Y_i-X) g_0(X)
	% but better to first update the hyper parameters given new observation z

	% THIS IS CALCULATED MANY TIMES! SHOULD BE STORED IN TABLE
	hyper = update_SS(z, hyperG0);
	% and then get X_i given these hyper parameters
	%[theta_mu_k, theta_Sigma_k] = normalinvwishrnd(hyper);
	R = sample_pdf(hyper);
	theta_mu_k = R.mu;
	theta_Sigma_k = R.Sigma;
    else
	% Sample old value
	% Remaining part is uniform as well
	u1 = u - p0;
	% find normal distribution in n with probabiliy as in Neal equation 3 at top
	% this is a normal trick to pick something out of an array of probability masses
	% make the array cumulative and pick then the first item beyond the randomly picked value
	%cs=cumsum(n/const);
	%	lastcsshouldbe1minusp0=cs(length(cs)) % should be 1-p0
	%fprintf('last item in array is indeed 1-p0=%i\n', cs(length(cs)));
	ind = find(cumsum(n/const)>=u1, 1);
	theta_mu_k = theta_mu_notk(:, ind);
	if (hyperG0.prior == 'NIW')
		theta_Sigma_k = theta_Sigma_notk(:, :, ind);
	else
		theta_Sigma_k = theta_Sigma_notk(ind);
	end
    end
end

function some_plot(zt, theta_mu, theta_Sigma, k, i, n, cmap)
    print_clusters=false;
    hold off
    z=zt(2:end,:);
    if (print_clusters)
       printf('--------------- New plot --------------\n');
    end
    % ugly... picks just one dimension and sorts... but yeah, it's only for visualization, so who cares? :-)
    ind = unique(theta_mu(1, :));
    cnt = hist(theta_mu(1,:),ind);
    [scnt,si]=sort(cnt,'descend');
    %si=si(1:length(ind)/3); % show only 1/3 largest clusters
    c = zeros(n, 1);
    U_mu = zeros(2, length(ind));
    for j=1:length(ind)
	ind2 = find(theta_mu(1, :)==ind(j));
	c(ind2) = j;
	U_mu(:, j) = theta_mu(:, ind2(1));
    end

    N=10;
    maxN=si(1:min(N,length(si))); % print N largest

    if (print_clusters)
        for j=maxN
            cluster=j
            ind2 = find(theta_mu(1, :)==ind(j));
            cluster_line_coefficients=theta_mu(:,ind2(1))
            sigma=theta_Sigma(ind2(1))
            nmbr=cnt(j)
       end
   end
    % print also one when tere is a single cluster
%    if (length(ind) == 1)
%	j=1
%	cluster_line_coefficients=theta_mu(:,j)
%	sigma=theta_Sigma(j)
 %   end

    %for j=1:length(ind)
    for j=maxN
	plot(z(1,c==j),z(2,c==j),'.','color',cmap(mod(15*j,63)+1,:), 'markersize', 15);
	hold on
	plot(U_mu(1,j),U_mu(2,j),'.','color',cmap(mod(15*j,63)+1,:), 'markersize', 30);
	plot(U_mu(1,j),U_mu(2,j),'ok', 'linewidth', 2, 'markersize', 10);
    end
    plot(z(1,k),z(2,k),'or', 'linewidth', 3)
    title(['i=' num2str(i) ',  k=' num2str(k) ', Nb of clusters: ' num2str(length(ind))]);
    xlabel('X')
    ylabel('Y')
    title(['i=' num2str(i) ',  k=' num2str(k) ', Nb of clusters: ' num2str(length(ind))]);
    xlabel('X')
    ylabel('Y')
    xlim([-1 1]*20);
    ylim([-1 1]*20);
    pause(.01)
end
