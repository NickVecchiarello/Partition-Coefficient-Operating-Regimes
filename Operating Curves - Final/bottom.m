function [Load_bottom] = bottom(input, b)
    
    %% Input Variables
    psi = input.psi;
    Ep = input.intra_ep;
    Kp_max = input.Kp_thresh;
    Clod = input.conc_thresh;                                       

    %% Calculate Loading
    Load_bottom = b.*(Clod*(1+psi)) + Clod*(Kp_max+Ep);
end