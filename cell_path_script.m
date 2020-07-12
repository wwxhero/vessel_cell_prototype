s_spec 		 = [0,	1,		2,		3,		4,		5,		6,		7];
alpha_x_spec = [0,	0.5,	1.0,	0.8,	1.2,	1.6,	1.4,	0.9];
alpha_y_spec = [0,	1.0,	0.5,	0.8,	1.4,	1.2,	0.9,	1.6];
alpha_z_spec = [0,	0.9,	1.0,	0.5,	0.8,	0.9,	1.4,	1.2];
bata_spec	 = [0.2,0.19,	0.18,	0.16,	0.21,	0.21,	0.25,	0.24];
pp_a_x = spline(s_spec, alpha_x_spec);
pp_a_y = spline(s_spec, alpha_y_spec);
pp_a_z = spline(s_spec, alpha_z_spec);
pp_b   = spline(s_spec, bata_spec);
s 			 = 0:0.025:6.999;

alpha_x = ppval(pp_a_x, s);
alpha_y = ppval(pp_a_y, s);
alpha_z = ppval(pp_a_z, s);
bata 	= ppval(pp_b, s);



[T_alpha, N_alpha, B_alpha] = ff_spline(pp_a_x, pp_a_y, pp_a_z, s_spec, s);

[~, n_s] = size(s);
alpha_s = zeros(3, n_s);
alpha_s(1, :) = alpha_x;
alpha_s(2, :) = alpha_y;
alpha_s(3, :) = alpha_z;

phi = 0:0.1:2*pi+0.1;
[S,PHI] = meshgrid(s, phi);
[m, n] = size(S);
T = zeros(m, n, 3);
for i_phi = 1 : m
	for i_s = 1 : n
		T(i_phi, i_s, :) = alpha_s(:, i_s) ...
						+ bata(1, i_s) * cos(phi(1, i_phi)) * N_alpha(:, i_s) ...%
						+ bata(1, i_s) * sin(phi(1, i_phi)) * B_alpha(:, i_s);%bata(1, i_s) *
	end
end

figure('Name','cell_path');
%surf(T(:, :, 1), T(:, :, 2), T(:, :, 3));
mesh(T(:, :, 1), T(:, :, 2), T(:, :, 3), 'FaceAlpha',0);


for k = 1 : 10
	k_1_1 = random('Uniform', -1, +1);
	k_1_2 = random('Uniform', -pi, +pi);
	k_2_1 = random('Uniform', -1, +1);
	k_2_2 = random('Uniform', -pi, +pi);
	rk		= bata .* abs(cos(k_2_1 * s + k_2_2));
	rk = rk - 0.01;
	thetak	= k_1_1 * s + k_1_2;

	p_x 	= alpha_x + rk .* cos(thetak) .* N_alpha(1) + rk .* sin(thetak) .* B_alpha(1);
	p_y 	= alpha_y + rk .* cos(thetak) .* N_alpha(2) + rk .* sin(thetak) .* B_alpha(2);
	p_z 	= alpha_z + rk .* cos(thetak) .* N_alpha(3) + rk .* sin(thetak) .* B_alpha(3);

	hold on
	plot3(alpha_x_spec, alpha_y_spec, alpha_z_spec, 'o'...
		, alpha_x, alpha_y, alpha_z, 'r' ...
		, p_x, p_y, p_z, 'b');
end

len = 0.1;

for idx_s = 1 : 10: n_s
	ori = alpha_s(:, idx_s);
	t = T_alpha(:, idx_s);
	n = N_alpha(:, idx_s);
	b = B_alpha(:, idx_s);
	p_t = ori + t;
	p_n = ori + n;
	p_b = ori + b;
	hold on
	vectarrow(ori, p_t, len, 'r');
	hold on
	vectarrow(ori, p_n, len, 'g');
	hold on
	vectarrow(ori, p_b, len, 'b');
end
