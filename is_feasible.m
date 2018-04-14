%% check feasibility%% min c'x s.t. ax = b, x>= 0%% max b'v s.t. a'v + s = c, s>=0function result = is_feasible (x, v, s, a, b, c)  r1 = a * x - b;r2 = a' * v + s - c;[r1; r2]if (max([r1; r2]) > 1e-30)  result = 0;else   result = 1;endif;
endfunction
