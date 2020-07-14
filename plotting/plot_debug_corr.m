function [f1] = plot_debug_corr(eps, match, center, varargin)
    
    main_title = "Matched (r) v Center (b)";
    
    if ~isempty(varargin)
        main_title = varargin{1};
    end
    
    W = size(match,2);

    f1 = figure;
    
    % compare psycho-curves between match vs central
    ax1_comp = subplot(2,2,[1,3]);
    hold on
    
    title(main_title);
    xlabel("stimulus location (eps)");
    ylabel("Prob answer correctly");
    ylim([0 1]);
    
    plot(eps, mean(match, 2), 'r', 'LineWidth', 4);
    plot(eps, mean(center, 2), 'b', 'LineWidth', 4);
    
    legend("matched","center", 'Location','best');
    
    
    % look at individual Ws for Match condition
    ax2_match = subplot(2,2,2);
    hold on
    
    title("Match by W input (dashed \rightarrow R=1)");
    xlabel("stimulus location (eps)");
    ylabel("P(Correct Answer)");
    ylim([0 1]);
    
    for i = 2:(W-1)
        plot(eps,match(:,i), 'LineWidth',4);
    end

    for i = [1,W]
        plot(eps,match(:,i), '--', 'LineWidth',4);
    end
    
    % look at individual Ws for Center condition
    ax3_center = subplot(2,2,4);
    hold on
    
    title("Center by W input (dashed \rightarrow R=1)");
    xlabel("Stimulus Location (eps)");
    ylabel("P(Correct Answer)");
    ylim([0 1]);
    
    for i = 2:(W-1)
        plot(eps,match(:,i), 'LineWidth',4);
    end

    for i = [1,W]
        plot(eps,match(:,i), '--', 'LineWidth',4);
    end

end