function H = trim_hessian(alpha,delta_e,params)

    dfdalphadalpha = 2 * (params.CMalpha^2 + (params.q^2 * params.S^2 * (params.CLalpha - (480.0 * alpha) / exp(400 * alpha^2))^2) / params.W^2 + (params.q * params.S * (-480.0 + 384000.0 * alpha^2) * (-params.W + params.q * params.S * (params.CL0 + params.CLalpha * alpha + params.CLdelta * delta_e + 0.6 / exp(400 * alpha^2)))) / (exp(400 * alpha^2) * params.W^2));
    dfdalphaddelta_e = 2 * params.CMalpha * params.CMdelta + (2 * params.CLdelta * params.q^2 * params.S^2 * (params.CLalpha - (480.0 * alpha) / exp(400 * alpha^2))) / params.W^2;
    dfddelta_edalpha = 2 * params.CMalpha * params.CMdelta + (2 * params.CLdelta * params.q^2 * params.S^2 * (params.CLalpha - (480.0 * alpha) / exp(400 * alpha^2))) / params.W^2;
    dfddelta_eddelta_e = 2 * (params.CMdelta^2 + (params.CLdelta^2 * params.q^2 * params.S^2) / params.W^2);
    
    H = [dfdalphadalpha dfdalphaddelta_e; dfddelta_edalpha dfddelta_eddelta_e];

end