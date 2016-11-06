% Gibbs sampler for Dirichlet Process Mixtures of Gaussians
% The base distribution G0 is Normal inverse Wishart of parameters given by
% hyperG0
% Reference: Algorithm 2 of Neal

% In the paper we have used this one, not the collapsed CRP.

function c_st = gibbsDPM_algo2(y, hyperG0, alpha, niter, doPlot)

    if doPlot
        figure('name','Gibbs sampling for DPM');
        colormap('default')
        cmap = colormap;
    end

    n = size(y,2);
    c_st = zeros(n, niter/2);

    % U_SS is a structure array where U_SS(k) contains the sufficient
    % statistics associated to cluster k
    if (hyperG0.prior == 'NIW')
        U_SS = struct('prior', 'NIW', 'mu', cell(n, 1), 'kappa', cell(n, 1), ...
        'nu', cell(n, 1), 'lambda', cell(n, 1));
    else
        U_SS = struct('prior', 'NIG', 'mu', cell(n, 1), 'a', cell(n, 1), 'b', cell(n, 1), 'Lambda', cell(n, 1));
    end

    % Just reserve enough tables/clusters
    m = zeros(1,200);
    c = zeros(n, 1);

    % Initialisation
    for k=1:n
        c(k) = ceil(30*rand); % Sample new allocation uniform
        m(c(k)) = m(c(k)) + 1;
        if m(c(k))>1
            U_SS(c(k)) = update_SS(y(:,k), U_SS(c(k)));
        else
            U_SS(c(k)) = update_SS(y(:,k), hyperG0);
        end
    end

    % Iterate over the tables (unique indices in customer allocation array c)
    % to get the first samples for U_R.
    ind = unique(c);
    for j=1:length(ind)
        R = sample_pdf(U_SS(ind(j)));
        U_R(ind(j)) = R;
    end

    % Iterations
    for i=2:niter
        % Update cluster assignments c
        for k=1:n
            % Remove data k from the partition
            m(c(k)) = m(c(k)) - 1;
	
            U_SS(c(k)) = downdate_SS(y, k, U_SS(c(k)));

            % Sample allocation of data k
            c(k) = sample_c(m, alpha, y(:,k), hyperG0, U_R);
            m(c(k)) = m(c(k)) + 1;
            if m(c(k))>1
                % There were already data items at this table
                % Update the sufficient statistics with the new data y(:,k)
                U_SS(c(k)) = update_SS(y, k, U_SS(c(k)));
            else
                % It is the first item at this table, calculate the sufficient
                % statistics using the base distribution G0 and the first data
                % item
                U_SS(c(k)) = update_SS(y, k, hyperG0);
                % Sample from H_i
                R = sample_pdf(U_SS(c(k)));
                U_R(c(k)) = R;
            end

            if doPlot==1
                some_plot(y, hyperG0, U_R, m, c, k, i, cmap)
            end
        end

        % Update cluster locations U, we just sample using the sufficient 
        % statistics U_SS
        ind = find(m);
        for j=1:length(ind)
            R = sample_pdf(U_SS(ind(j)));
            U_R(ind(j)) = R;
        end

        if i>niter/2
            c_st(:, i-niter/2) = c;
        end

        if doPlot==2
            some_plot(y, hyperG0, U_R, m, c, k, i, cmap)
        end

    end
end

% Sample a new value for a coordinate to a table given an observation z. This
% can be an existing or new table coordinate. It corresponds to Neal's Eq. 3.6.
%
% The posterior predictive "pred" is called H_i, the posterior based on G_0
% (hyperG0) and single observation y_i (z). Here we do not yet sample from it,
% but we need this conditional probability to decide if we later want to sample
% from it. Note, that this is actually called the prior predictive
% distribution.
% The likelihoods are denoted by F(y_i,\theta_c)
function K = sample_c(m, alpha, z, hyperG0, U_R)

    c = find(m~=0); % gives indices of non-empty clusters

    n = m(c).*likelihoods(hyperG0.prior, repmat(z, 1, length(c)), U_R(c) )';

    n0 = pred(z, hyperG0);
    const = sum(n) + alpha*n0;

    p0 = alpha*n0/const; % probability of sampling a new item
    u=rand(1);
    if u<p0
        K = find(m==0, 1); % sample new
    else
        u1 = (u-p0);
        ind = find(cumsum(n/const)>=u1, 1 );
        K = c(ind);
    end
end

function some_plot(z, hyperG0, U_R, m, c, k, i, cmap)

    % If it's 
    if strcmp(hyperG0.prior, 'NIG')
        z=z(2:end,:);
    end

    ind=find(m);
    hold off;
    for j=1:length(ind)
        plot(z(1,c==ind(j)),z(2,c==ind(j)),'.','color',cmap(mod(5*ind(j),63)+1,:), 'markersize', 15);
        hold on
        %plot(U_mu(1,ind(j)),U_mu(2,ind(j)),'.','color',cmap(mod(5*ind(j),63)+1,:), 'markersize', 30);
        %plot(U_mu(1,ind(j)),U_mu(2,ind(j)),'ok', 'linewidth', 2, 'markersize', 10);
        plot(U_R(ind(j)).mu(1),U_R(ind(j)).mu(2),'.','color',cmap(mod(5*ind(j),63)+1,:), 'markersize', 30);
        plot(U_R(ind(j)).mu(1),U_R(ind(j)).mu(2),'ok', 'linewidth', 2, 'markersize', 10);
    end
    plot(z(1,k),z(2,k),'or', 'linewidth', 3)
    title(['i=' num2str(i) ',  k=' num2str(k) ', Nb of clusters: ' num2str(length(ind))]);
    xlabel('X');
    ylabel('Y');
    xlim([-1 1]*20);
    ylim([-1 1]*20);
    %pause(.01)

    save_plot=true;
    if (save_plot)
        name=sprintf('output/image%04i.jpg', i);
        name
        print(name, '-djpg');
    end
end
