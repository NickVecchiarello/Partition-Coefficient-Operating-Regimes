function [tr_app, C_equil, Sep_Fact] = equilibrium(input)
%% Calculates equilibrium concentrations (batch), separation factors (batch), and column retention times (column) as a function of Ke and qm for a given column loading and phase ratio
%% Import input variables
    charge = input.charge;                                                                              %Resin loading per well (mg protein/mL resin)       
    Vl = input.Vl;                                                                                      %Liquid volume per well (uL)
    Vr = input.Vr;                                                                                      %Resin volume per well (uL)
    Vh = input.Vh;                                                                                      %Holdup volume per well (uL)
    intra_ep = input.intra_ep;                                                                          %Intraparticle porosity
    inter_ep = input.inter_ep;                                                                          %Interparticle porosity
    V = input.V;                                                                                        %Column volume (mL)
    d = input.d;                                                                                        %Column diameter (cm)                                                                                
    Q = input.Q;                                                                                        %Column volumetric flow rate (mL/min)
    Ke_span = input.Ke_span;                                                                            
    qm_span = input.qm_span;
    Ke_max = input.Ke_max;
    qm_max = input.qm_max;
    
    %% Calculated parameters
    phi_col = (1-inter_ep)/inter_ep;                                                                    %Phase ratio of the column
    L = V/(3.14*(d/2)^2);                                                                               %Column length (cm)
    u = Q/(3.14*((d/2)^2));                                                                             %Linear flow rate (cm/min)
    v = u/inter_ep;                                                                                     %Interstitial velocity (cm/min)
    phi = Vr/Vl;                                                                                        %Phase ratio for batch
    Co = charge*phi;                                                                                    %Initial concentration needed to charge resin (mg/mL)
    psi = Vh/Vl;                                                                                        %Holdup ratio
    
    %% Output calculations
    Ke_temp = linspace(0, Ke_max, Ke_span);                                                             %Creates gridpoints for Ke values
    qm_temp = linspace(0,qm_max/(1-intra_ep),qm_span);                                                  %Creates gridpoints for Qm values
    slope_app = zeros(Ke_span,qm_span);                                                                 %Initialize initial isotherm slope
    C_equil = zeros(Ke_span,qm_span);                                                                   %Initialize protien equilibrium concentration
    Sep_Fact = zeros(Ke_span,qm_span);                                                                  %Initialize Separation Factor
    
    for i = 1:Ke_span
        Ke = Ke_temp(i);
        for j = 1:qm_span
            qm = qm_temp(j);
            a = Ke*(1+(intra_ep*phi)+(psi*Ke));
            b = 1+(intra_ep*phi)+psi+Ke*((qm*phi*(1-intra_ep))-Co);
            c = -Co;
            root = roots([a,b,c]);
            C_app = root(root>=0);                                                                      %Equilibrium concentration by mass balance from Langmuir isotherm
            q_app = (Ke*qm*C_app/(1+Ke*C_app));                                                         %Bound concentration by mass balance from Langmuir isotherm
            slope_app(i,j) = q_app/C_app;                                                               %Apparent slope of isotherm (secant slope)
            C_equil(i,j) = C_app;
            Sep_Fact(i,j) = 1/(1+Ke*C_app);                                                             %Separation Factor
        end
    end
    tr_app = (L/v)*(1+phi_col*slope_app+(phi_col*intra_ep*C_equil));                                    %Retention time on column
end

