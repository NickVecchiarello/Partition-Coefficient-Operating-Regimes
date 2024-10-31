function [ref_iso] = import_iso(x)
%Inputs all reference isotherm parameters and data from 'Reference
%Isotherms.xlsx'
import = importdata('Reference Isotherms.xlsx');
data = import.data;
textdata = import.textdata;

%Classify isotherms as SMA or Langmuir
isotype = zeros(1,size(data,1));
for i = 1:size(data,1)
        if strmatch(textdata(i+2,3), 'SMA') == 1
            isotype(i) = 1;                                                %Isotype: 1=SMA, 2=Langmuir
        elseif strmatch(textdata(i+2,3), 'Langmuir') == 1
            isotype(i) = 2;
        end
end

%% Compile SMA data values from Excel
if sum(isotype == 1) > 0
    ref_iso.SMA.MW = data(find(isotype == 1),1);                               %Molecular Weights
    ref_iso.SMA.Ke = data(find(isotype == 1),3);                               %Equilibrium Constants    
    ref_iso.SMA.q0 = data(find(isotype == 1),4);                               %Ionic Capacities   
    ref_iso.SMA.z = data(find(isotype == 1),5);                                %Characteristic Charges   
    ref_iso.SMA.sig = data(find(isotype == 1),6);                              %Steric Factor  
else

end
%% Compile Langmuir data values from Excel
if sum(isotype == 2) > 0
    ref_iso.Langmuir.kl = data(find(isotype == 2),7);                          %Langmuir Parameters 
    ref_iso.Langmuir.qm = data(find(isotype == 2),8);                          %Equilbirium Max Capacities
else
end