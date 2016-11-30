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

    % n is the number of observations
    [d n] = size(y);

    if (n < d)
        warning("Are you sure that y is ordered with each observation as a column?");
    end

    % c_st is statistics to be returned
    c_st = zeros(n, niter/2);

    % U_SS is a structure array where U_SS(k) contains the sufficient
    % statistics associated to cluster k
    if (hyperG0.prior == 'NIW')
        U_SS = struct('prior', 'NIW', ...
            'mu', cell(n, 1), ...
            'kappa', cell(n, 1), ...
            'nu', cell(n, 1), ...
            'lambda', cell(n, 1));
    else
        U_SS = struct('prior', 'NIG', 'mu', cell(n, 1), 'a', cell(n, 1), 'b', cell(n, 1), 'Lambda', cell(n, 1));
    end

    % Just reserve enough tables/clusters
    % m is the number of observations at table k: #obs = m(k) 
    m = zeros(1,200);
    % c is assignment between observation i and table k: c(i) = k
    c = zeros(n, 1);

    % Initialisation
    for i=1:n
        % Sample new allocation uniform, using rand
        k = ceil(30*rand); 
        c(i) = k;
        m(k) = m(k) + 1;
        % Update sufficient statistics of table k with information from data i
        if m(k)>1
            U_SS(k) = update_SS(y, i, U_SS(k));
        else
            U_SS(k) = update_SS(y, i, hyperG0);
        end
    end

    % ind is indices to tables (by running find on tables m)
    ind = find(m);
    for j=1:length(ind)
        k = ind(j);
        R = sample_pdf(U_SS(k));
        U_R(k) = R;
    end

    % t is used to indicate iterations (not time)
    for t=2:niter
        for i=1:n
            % Update cluster assignment c(i) 
            k = c(i);

            % Remove data i from the partition
            m(k) = m(k) - 1;
	
            U_SS(k) = downdate_SS(y, i, U_SS(k));

            % New cluster assignment for observation i
            c(i) = sample_c(m, alpha, y(:,i), hyperG0, U_R);
            k = c(i);

            % Add data i to new found table
            m(k) = m(k) + 1;

            % Update or generate sufficient statistics
            if m(k)>1
                % There were already data items at this table
                % Update the sufficient statistics with the new data y(:,i)
                U_SS(k) = update_SS(y, i, U_SS(k));
            else
                % It is the first item at this table, calculate the sufficient
                % statistics using the base distribution G0 and the data item
                U_SS(k) = update_SS(y, i, hyperG0);
                % Sample from H_i
                R = sample_pdf(U_SS(k));
                U_R(k) = R;
            end

            if doPlot==1
                some_plot(y, hyperG0, U_R, m, c, i, t, cmap)
            end
        end

        % Update cluster locations U, we just sample using the sufficient 
        % statistics U_SS
        ind = find(m);
        for j=1:length(ind)
            k = ind(j);
            R = sample_pdf(U_SS(k));
            U_R(k) = R;
        end

        if t>niter/2
            c_st(:, t-niter/2) = c;
        end

        if doPlot==2
            some_plot(y, hyperG0, U_R, m, c, k, t, cmap)
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

    Z = repmat(z, 1, length(c));
    l = likelihoods(hyperG0.prior, Z, U_R(c) );

    if size(l, 2) != length(c) || size(l,1) != 1
        error("Likelihood vector should be of size 1 x C (row-vector)");
    end
    n = m(c) .* l;

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
        disp(name);
        print(name, '-djpg');
    end
end
