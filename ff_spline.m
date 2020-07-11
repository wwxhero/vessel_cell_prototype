function [T_a, N_a, B_a] = ff_spline(pp_a_x, pp_a_y, pp_a_z, s_spec, s)
	[~, n_s] = size(s);
	T_a = zeros(3, n_s);
	N_a = zeros(3, n_s);
	B_a = zeros(3, n_s);
    epsilon = 0.0001;
	for i_s = 1 : n_s
		s_i = s(1, i_s);
		[~, n_spec] = size(s_spec);
		i = 0;
		for i_spec = 1: n_spec - 1
			if s_spec(1, i_spec) <= s_i ...
				&& 			 s_i < s_spec(1, i_spec + 1)
				i = i_spec;
				break;
			end
		end


		coefs_x = pp_a_x.coefs(i, :);
		coefs_y = pp_a_y.coefs(i, :);
		coefs_z = pp_a_z.coefs(i, :);
		s_l = s_i - s_spec(i);
		a_prime_x = 3 * coefs_x(1) * s_l^2 ...
					+ 2 * coefs_x(2) * s_l^1 ...
					+ 1 * coefs_x(3);
		a_prime_y = 3 * coefs_y(1) * s_l^2 ...
					+ 2 * coefs_y(2) * s_l^1 ...
					+ 1 * coefs_y(3);
		a_prime_z = 3 * coefs_z(1) * s_l^2 ...
					+ 2 * coefs_z(2) * s_l^1 ...
					+ 1 * coefs_z(3);
		a_prime = [a_prime_x, a_prime_y, a_prime_z];
		T_a_i = (a_prime/norm(a_prime))';
		a_prime_prime_x = 3 * 2 * coefs_x(1) * s_l ...
						+ 2 * 1 * coefs_x(2);
		a_prime_prime_y = 3 * 2 * coefs_y(1) * s_l ...
						+ 2 * 1 * coefs_y(2);
		a_prime_prime_z = 3 * 2 * coefs_z(1) * s_l ...
						+ 2 * 1 * coefs_z(2);
		a_prime_prime = [a_prime_prime_x, a_prime_prime_y, a_prime_prime_z];

		B_a_i_scaled = cross(a_prime, a_prime_prime);
		B_a_i = (B_a_i_scaled/norm(B_a_i_scaled))';
		dot_product = dot(T_a_i, B_a_i);
		assert(-epsilon < dot_product && dot_product < epsilon)
		N_a_i = cross(B_a_i, T_a_i);
		T_a(:, i_s) = T_a_i;
		N_a(:, i_s) = N_a_i;
		B_a(:, i_s) = B_a_i;
	end
end