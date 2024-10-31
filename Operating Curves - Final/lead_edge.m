%% Calculates reference isotherms points at the leading edge (worst case) conditions

function y = lead_edge(qm, Ke, Rsq, ref_isotherm)
    %% Find left-most point on qm vs Ke plot
    m=[qm;Ke;Rsq]';
    [~,start_ind] = min(qm);
    y(1) = start_ind;
    %% For each nearest neighbor
    for i = 1:length(m)
        %a = pdist2(m, m(start_ind,:));
        a = abs(m(:,1)-m(start_ind,1));
        [sor, sor_ind] = sort(a);
        
        %% Find closest point that increases Ke      
        for j = 2:length(sor)
            if m(sor_ind(j),2) >= m(start_ind,2) && m(sor_ind(j),1) >= m(start_ind,1) && m(start_ind,3) >= ref_isotherm.Rsqlim
                start_ind = sor_ind(j);
                break
            end
        end
    
        if start_ind ~= y(end)
            y(i+1) = start_ind;
        else
            break
        end

    end
end