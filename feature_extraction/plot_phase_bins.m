% For 8 phase bins (trough, early rising, rising, late rising, peak,
% early falling, falling, late falling)
plot_phase_bins_func(8);

% For 4 phase bins (trough, rising, peak, falling)
plot_phase_bins_func(4);

function plot_phase_bins_func(nBins)
    % nBins: 4 or 8 recommended

    figure; hold on; axis equal; box on;
    title(sprintf('Phase bins on circle (%d bins)', nBins));

    % Draw unit circle
    th = linspace(-pi, pi, 400);
    plot(cos(th), sin(th), 'k'); % circle outline

    % Compute centered bins (true trough/peak alignment)
    binWidth = 2*pi / nBins;

    % centers at -pi, -pi + binWidth, ..., +pi - binWidth
    centers = -pi : binWidth : pi - binWidth;

    % edges halfway between centers
    edges = centers - binWidth/2;
    edges(end+1) = edges(end) + binWidth;
    edges(1) = [];

    % Draw wedges
    for b = 1:nBins
        th_bin = linspace(edges(b), edges(b), 50);
        x = [0, cos(th_bin), 0];
        y = [0, sin(th_bin), 0];

        % filled wedge
        h = fill(x, y, 'k');
        h.FaceAlpha = 0.12;
        h.EdgeColor = 'none';
    end

      % Define semantic labels
    if nBins == 4
        phase_labels = ["Trough (1)", "Rising (2)", "Peak (3)", "Falling (4)"];
        desired_index = [1 2 3 4];  % matches order
    elseif nBins == 8
        phase_labels = ["Trough (1)", ...
                        "Early Rising (2)", ...
                        "Rising (3)", ...
                        "Late Rising (4)", ...
                        "Peak (5)", ...
                        "Early Falling (6)", ...
                        "Falling (7)", ...
                        "Late Falling (8)"];
        desired_index = 1:8;
    else
        error('Unsupported nBins. Use 4 or 8.');
    end
    
    % Plot center rays (THICK BLACK — bin meaning)
    for b = 1:nBins
        c = centers(b);
        plot([0 cos(c)], [0 sin(c)], 'k', 'LineWidth', 2)

        % label center
        r_label = 1.2;
        text(r_label*cos(c), r_label*sin(c), phase_labels(b), ...
            'HorizontalAlignment','center','FontWeight','bold');
    end

    % Plot edges (DOTTED — show boundaries)
    for b = 1:length(edges)
        e = edges(b);
        plot([0 cos(e)], [0 sin(e)], '--', 'Color', [0.2 0.2 0.2]);
    end

    % Formatting
    xlim([-1.4 1.4]);
    ylim([-1.4 1.4]);
    xlabel('cos(\theta)'); ylabel('sin(\theta)');
    set(gca,'XAxisLocation','origin','YAxisLocation','origin');
end
