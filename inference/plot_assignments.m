% Plot mean values
function plot_assignments(z, hyperG0, U_R, m, c, k, i, cmap)
    switch (hyperG0.prior)
    case { 'NIG', 'DPM_Seg' }
        z=z(2:end,:);
    end
    
    % TODO: make into argument
    draw_lines = false;

    ind = find(m);
    hold off;

    for j=1:length(ind)
        color = cmap(mod(5*ind(j),63)+1,:);
        mu = U_R(ind(j)).mu;

        switch (hyperG0.prior)
        case { 'DPM_Seg' }
            x_a(j) = U_R(ind(j)).a;
            x_b(j) = U_R(ind(j)).b;
        otherwise
            x_a(j)=-25;
            x_b(j)=+25;
        end

        if (draw_lines)
            y_a(j) = [1 x_a(j)]*mu;
            y_b(j) = [1 x_b(j)]*mu;
            plot([x_a(j) x_b(j)], [y_a(j) y_b(j)], '-', 'color', color, 'linewidth', 5);
            hold on
        end

        switch (hyperG0.prior)
        case { 'DPM_Seg' }
            plot(x_a(j), y_a(j), '.', 'color',color, 'markersize', 30);
            plot(x_a(j), y_a(j), 'ok', 'linewidth', 2, 'markersize', 10);
            plot(x_b(j), y_b(j), '.', 'color',color, 'markersize', 30);
            plot(x_b(j), y_b(j), 'ok', 'linewidth', 2, 'markersize', 10);
        case 'NIW'
            plot(mu(1), mu(2), '.', 'color', color, 'markersize', 30);
            plot(mu(1), mu(2), 'ok', 'linewidth', 2, 'markersize', 10);
        end
        
        hold on

        plot(z(1,c==ind(j)), z(2,c==ind(j)), '.', 'color', color, 'markersize', 15);
        
        if (draw_lines)
            cust=find(c == ind(j));
            for f = 1:length(cust)
              plot([x_a(j) z(1,cust(f))],[y_a(j) z(2,cust(f))], '-', 'color', color);
              plot([x_b(j) z(1,cust(f))],[y_b(j) z(2,cust(f))], '-', 'color', color);
            end
        end
    end
    plot(z(1,k), z(2,k), 'or', 'linewidth', 3)
    title(['i=' num2str(i) ',  k=' num2str(k) ', Nb of clusters: ' num2str(length(ind))]);
    xlabel('X');
    ylabel('Y');
    %y_max=max([z(2,:),y_a,y_b, 25]);
    %y_min=min([z(2,:),y_a,y_b, -25]);
    %x_max=max([z(1,:),x_a,x_b, 25]);
    %x_min=min([z(1,:),x_a,x_b, -25]);
    y_max = 25;
    y_min = -25;
    x_max = 25;
    x_min = -25;
    xlim([x_min x_max]);
    ylim([y_min y_max]);

    pause(.01)

    save_plot=true;
    if (save_plot)
        name=sprintf('output/image%04i.jpg', i);
        disp(name);
        print(name, '-djpg');
    end
end

