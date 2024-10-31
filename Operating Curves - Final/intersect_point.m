%% Calculates point of intersection between Kp_max curve and reference isotherm

function [qm_star, Kl_star] = intersect_point(frnt_pts, Ke, qm, input)
            ref = [qm(frnt_pts); Ke(frnt_pts)];
            Kp_max = input.Kp_thresh;
            fun = @(x) pchip(ref(1,:),ref(2,:), x) - (Kp_max/x);
            x0 = 30;
            qm_star = fzero(fun, x0);
            Kl_star = Kp_max/qm_star;

end