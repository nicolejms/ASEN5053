function [alpha,beta,n] = FindVars(mu,a,s,c,lesserTOF,path)
    alpha = 2*asin(sqrt(s/(2*a)));
    beta =  2*asin(sqrt((s-c)/(2*a)));
    n = sqrt(mu/(a^3));

    if (path=="long")
        beta = -beta;
    end
    if not (lesserTOF)
       alpha = 2*pi - alpha; 
    end
end