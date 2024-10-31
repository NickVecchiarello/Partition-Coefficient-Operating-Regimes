%% Fits and compiles the reference isotherm data and compresses reference isotherms to within Ke_max and qm_max ranges

function [Ke_out, qm_out, Rsq_out, h] = iso(ref_isotherm, input)
    h = waitbar(0, {'Initializing'; '0%'});
    input.Ke_max = 10;                                                      %Max Ke tested
    input.qm_max = 300;                                                     %Max qm tested                                    
    opt = optimoptions('lsqcurvefit','Display','off'); 
    for i = 1:length(ref_isotherm.SMA.MW)
        if isfield(ref_isotherm, 'SMA') ==1
        %% Isotherm inputs
            Ke = ref_isotherm.SMA.Ke(i);
            qo = ref_isotherm.SMA.q0(i);
            z = ref_isotherm.SMA.z(i);
            sig = ref_isotherm.SMA.sig(i);
            MW = ref_isotherm.SMA.MW(i);
            
        %% Salt inputs
            start_salt = ref_isotherm.start_salt;
            end_salt = ref_isotherm.end_salt;
            n_salt = ref_isotherm.n_salt;
    
        %% Concentration inputs
            start_c = ref_isotherm.start_c;
            end_c = ref_isotherm.end_c;
            n_c = ref_isotherm.n_c;                                             
            Na = linspace(start_salt,end_salt,n_salt);                       %Salt concentrations for SMA
            c_temp=exp(linspace(log(start_c),log(end_c),n_c));               %Log spacing for protein concentration points to accurately capture linear part of isotherm 
            c_temp = 1000*c_temp/MW;                                         %Concentration of protein in mM
    
        %% Calculated SMA Isotherms
            Q_temp = zeros(length(Na), length(c_temp));
            for j=1:length(Na)
                Na_temp = Na(j);
                for k = 1:length(c_temp) 
                c = c_temp(k);
                opts = optimset('Diagnostics','off', 'Display','off');
                fun = @(q)MassAction(q, Ke, qo, z, Na_temp, sig, c);
                Q_temp(j,k) = fsolve(fun, 0, opts);
                end
            end

        %% Fit SMA Isotherms to Langmuir Isotherms
            Q = Q_temp*(MW/1000);
            C = c_temp*(MW/1000);
            for m=1:size(Q_temp,1)
                fun = @(x, C)(x(1).*x(2)*C./(1+x(2).*C));
                x0 = [100, 1];
                x = lsqcurvefit(fun, x0, C, Q(m,:), [], [], opt);
                qm_out(i,m) = x(1);
                Ke_out(i,m) = x(2);
                yfit = fun(x, C);                                               % Estimated  Regression Line
                SStot = sum((Q(m,:) - mean(Q(m,:))).^2);                                  % Total Sum-Of-Squares
                SSres = sum((Q(m,:) - yfit(:)').^2);                               % Residual Sum-Of-Squares
                Rsq(i,m) = 1-SSres/SStot;                                         % R^2
            end
        else 
        end
    C_out(i,:) = C;
    Q_out(i,:,:) = Q;
    frac = (i*1.43)/(1.43*length(ref_isotherm.SMA.MW)+7.4);
    msg = 'Calculating Reference Isotherms';
    percent = round(100*frac);
    close all
    waitbar(frac, h, {msg; strcat(num2str(percent), '%')});
    end

    %% Removes Ke and qm values that exceed the evaluated range of Ke and qm 
    Ke_out = reshape(Ke_out,[1,numel(Ke_out)]);
    qm_out = reshape(qm_out, [1,numel(qm_out)]);
    Rsq_out = reshape(Rsq,[1,numel(Rsq)]);
    qm_indicies = find(qm_out>input.qm_max);
    Rsq_ind = find(Rsq_out<=ref_isotherm.Rsqlim);
    comb_indicies = unique([qm_indicies, Rsq_ind]);
    Ke_out(comb_indicies)=[];
    qm_out(comb_indicies)=[];
    Rsq_out(comb_indicies)=[]; 

if isfield(ref_isotherm, 'Langmuir') ==1
    Ke_out =[Ke_out, ref_isotherm.Langmuir.kl'];
    qm_out =[qm_out, ref_isotherm.Langmuir.qm'];
    Rsq_out = [Rsq_out, ones([numel(ref_isotherm.Langmuir.kl), 1])'];
else
end
    
end
