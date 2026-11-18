function return_data = vis_sub_on_s1()
%VIS_SUB_ON_S1 provides a visualization of subtraction on S^1. My
%theoretical form is x-y=p_s(argmin_\bar{x}|\bar{x}-\bar{y}|-\bar{y}) where
%p_s is a covering map R\toS^1 and \bar notation selects representatives in
%R. This form simplifies in the computer by abusing the equivalence of
%true/false and 0/1.

    % Discretize space
    x = linspace(-1/2,1/2,1000);
    y = linspace(-1/2,1/2,1000);

    % Perform calculation
    d_r = x-y.';
    d_s = ((d_r < -1/2) - (d_r > 1/2)) + d_r;

    % Instantiate figure
    circ_diff_frame = figure(1);
    clf(circ_diff_frame);
    ax = axes(circ_diff_frame);

    % Plot results
    surf(ax, x,y,d_s);
    xlabel(ax,'x');
    ylabel(ax,'y');
    title(ax,'Graph of v(x,y) using S^1 ~ [-0.5,0.5)')
    shading interp;
    colorbar
    view(2)

    return_data = 1;

end