function plotresults(t_train, t_test, y_train, y_test, Omega_train, Omega_test, Omega_train_init, Omega_test_init, optkappa, kappa_init)
    figure(1)
    tlo = tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    % === PLOT 1: Non-Optimal Reservoir (Train) ===
    ax1 = nexttile;
    hold on;
    plot(t_train, (Omega_train_init * kappa_init'), '-', 'DisplayName', 'Non-Optimal Reservoir', 'MarkerSize', 10, 'LineWidth', 4);
    plot(t_train, y_train, '--r', 'DisplayName', 'Training Signal', 'LineWidth', 4);
    xlabel('Time (t)');
    ylabel('Amplitude');
    box on;
    set(gca, 'fontsize', 20);
    legend('show', 'NumColumns', 1);
    xlim([20 30])
    
    % === PLOT 2: Non-Optimal Reservoir (Test) ===
    ax2 = nexttile;
    hold on;
    plot(t_test, (Omega_test_init * kappa_init'), '-', 'DisplayName', 'Non-Optimal Reservoir', 'MarkerSize', 10, 'LineWidth', 4);
    plot(t_test, y_test, '--r', 'DisplayName', 'Test Signal', 'LineWidth', 4);
    xlabel('Time (t)');
    ylabel('Amplitude');
    box on;
    set(gca, 'fontsize', 20);
    legend('show', 'NumColumns', 1);
    
    % === PLOT 3: Optimal Reservoir (Train) ===
    ax3 = nexttile;
    hold on;
    plot(t_train, Omega_train * optkappa, '-', 'DisplayName', 'Optimal Reservoir', 'MarkerSize', 10, 'LineWidth', 4);
    plot(t_train, y_train, '--r', 'DisplayName', 'Training Signal', 'LineWidth', 4);
    xlabel('Time (t)');
    ylabel('Amplitude');
    box on;
    set(gca, 'fontsize', 20);
    legend('show', 'NumColumns', 1);
    xlim([20 30])
    
    % === PLOT 4: Optimal Reservoir (Test) ===
    ax4 = nexttile;
    hold on;
    plot(t_test, Omega_test * optkappa, '-', 'DisplayName', 'Optimal Reservoir', 'MarkerSize', 10, 'LineWidth', 4);
    plot(t_test, y_test, '--r', 'DisplayName', 'Test Signal', 'LineWidth', 4);
    xlabel('Time (t)');
    ylabel('Amplitude');
    box on;
    set(gca, 'fontsize', 20);
    legend('show', 'NumColumns', 1);
    
    % === (a), (b), (c), (d) Labels ===
    labels = {'(a)', '(b)', '(c)', '(d)'};
    axes_list = {ax1, ax2, ax3, ax4};
    
    for i = 1:4
        ax = axes_list{i};
        xLimits = xlim(ax);
        yLimits = ylim(ax);
        xText = xLimits(1) + 0.02 * range(xLimits);
        yText = yLimits(2) - 0.05 * range(yLimits);
        text(ax, xText, yText, labels{i}, 'FontWeight', 'bold', 'FontSize', 20);
    end
end