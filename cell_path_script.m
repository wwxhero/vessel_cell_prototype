s_spec 		 = [0,	1,		2,		3,		4,		5,		6,		7];
alpha_x_spec = [0,	0.5,	1.0,	0.8,	1.2,	1.6,	1.4,	0.9];
alpha_y_spec = [0,	1.0,	0.5,	0.8,	1.4,	1.2,	0.9,	1.6];
alpha_z_spec = [0,	0.9,	1.0,	0.5,	0.8,	0.9,	1.4,	1.2];
bata_spec	 = [0.2,0.19,	0.18,	0.16,	0.21,	0.21,	0.25,	0.24];
pp_a_x = spline(s_spec, alpha_x_spec);
pp_a_y = spline(s_spec, alpha_y_spec);
pp_a_z = spline(s_spec, alpha_z_spec);
pp_b   = spline(s_spec, bata_spec);
s 			 = 0:0.025:7;

alpha_x = ppval(pp_a_x, s);
alpha_y = ppval(pp_a_y, s);
alpha_z = ppval(pp_a_z, s);
bata 	= ppval(pp_b, s);

rk		= bata .* abs(cos(s/2));
rk = rk - 0.01;
thetak	= 2 * s;

[~, N_alpha, B_alpha] = ff_spline(pp_a_x, pp_a_y, pp_a_z, s_spec, s);

p_x 	= alpha_x + rk .* cos(thetak) .* N_alpha(1) + rk .* sin(thetak) .* B_alpha(1);
p_y 	= alpha_y + rk .* cos(thetak) .* N_alpha(2) + rk .* sin(thetak) .* B_alpha(2);
p_z 	= alpha_z + rk .* cos(thetak) .* N_alpha(3) + rk .* sin(thetak) .* B_alpha(3);

[~, n_s] = size(s);
alpha_s = zeros(3, n_s);
alpha_s(1, :) = alpha_x;
alpha_s(2, :) = alpha_y;
alpha_s(3, :) = alpha_z;
p_s	= zeros(3, n_s);
p_s(1, :) = p_x;
p_s(2, :) = p_y;
p_s(3, :) = p_z;

d_alpha_p = zeros(1, n_s);
err = zeros(1, n_s);
for i_s = 1 : n_s
	d_alpha_p(1, i_s) = norm(alpha_s(:, i_s) - p_s(:, i_s));
	err(1, i_s) = (bata(1, i_s) - d_alpha_p(1, i_s)); %make sure the cell is in vessel
end

figure('Name','cell_path');
plot3(alpha_x_spec, alpha_y_spec, alpha_z_spec, 'o'...
	, alpha_x, alpha_y, alpha_z, 'r' ...
	, p_x, p_y, p_z, 'b');



figure('Name','path_width');
bata_u = bata;
bata_l = -bata_u;
plot(s, bata_u, 'r' ...
    , s, bata_l, 'r' ...
    , s, d_alpha_p, 'b');


