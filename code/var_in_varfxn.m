function v_of_v= var_in_varfxn( p, N0, V0, t, N )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

v_of_v= (1/N).*(V4_fxn_ODE(p,N0, V0, t)-(((N-3)./(N-1)).*(V_fxn(p, t, N0, V0).^2)));
for j = 1:length(v_of_v)
    if v_of_v(j)<0
        v_of_v(j)=1e-5;
end
end

