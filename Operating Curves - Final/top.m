function [b, Load_top] = top(input, qm_star, Ke, qm, frnt_pts, ref_isotherm, h)
    wax = findobj('type','axes');
    tax = get(wax,'title');
    set(tax,'fontsize',15)
    frac = (length(ref_isotherm.SMA.MW)*1.43)/(1.43*length(ref_isotherm.SMA.MW)+7.4);
    totaltime = (1.43*length(ref_isotherm.SMA.MW)+7.4);
    msg = 'Calculating & Plotting Operating Curves';
    percent = round(100*frac);
    waitbar(frac, h, {msg; strcat(num2str(percent), '%')});
    %% Input Variables
    R = input.sep_fact_thresh;
    psi = input.psi;
    Ep = input.intra_ep;
    qm_base = linspace(0, qm_star, 10^3);
    b = logspace(log10(input.min_Vr_range), log10(input.max_Vr_range), 10^3);
    L = linspace(0, 20, 10^3);                                           % Load Range
    Load_top = zeros([1, length(L)]);
    
    %% Reference Isotherm Curve
    Kl_ref = pchip([0, qm(frnt_pts)], [0, Ke(frnt_pts)], qm_base);
    
    for i = 1:length(b)
        if ismember(i, linspace(1,length(b), 100)) == 1
        percent = round(100*(frac+(i/length(b))*(7.4/totaltime)));
        waitbar(frac+(i/length(b))*(7.4/totaltime), h,  {msg; strcat(num2str(percent), '%')});
        end
        for j = 1:length(L)
            %% Separation Factor Curve
            a = (b(i)*(1-R))*(1+Ep*(1/b(i))+psi)./(R*(L(j)-qm_base*(1-R)));
            a(a<0) = max(a);
            Kl_sep(j,:) = a;
        end
        if isempty(max(find(all(Kl_ref <= Kl_sep, 2)))) == 1 
            Load_ind(i) = 0; 
        else
        Load_ind(i)= max(find(all(Kl_ref <= Kl_sep, 2)));
        end
    end
    Load_top(find(Load_ind>0)) = L(Load_ind(Load_ind>0));
    delete(h)
end