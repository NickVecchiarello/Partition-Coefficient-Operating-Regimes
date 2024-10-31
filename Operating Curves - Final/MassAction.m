function F = MassAction(q,Ke, qo, z, Na_temp, sig, c)
F = (Ke*(qo-(z+sig).*q)^z.*c/(Na_temp^z))-q;
end

