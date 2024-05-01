function A = lambert(mu,a,s,c,angle,TOF,lesserTOF,path)

    deltaa = 100;
    options = optimoptions('fsolve','Display','off', 'FunctionTolerance',1e-7,'MaxIterations',100);
    a0 = a + deltaa;
    % error tolerance
    tol = 1e-7;
    diff = 1000;
    i = 1;
    while (diff > tol)  
        % Lambert's equation
        anew = fsolve(@(a)Lamberts(mu,a,s,c,TOF,lesserTOF,path),a0,options);
    
        % recalculate
        [alpha,beta,n] = FindVars(mu,anew,s,c,lesserTOF,path);
        % recheck TOF
        TOFnew =  1/n*((alpha-beta)-(sin(alpha)-sin(beta)));
        diff = abs(TOFnew-TOF) ;

        if i==100
            disp("could not find value for Lambert")
            break;
        end
        a0 = anew + 10;
        i = i+1;
    end
    disp("Number of iterations: " + i);
    A = anew;
end


function F = Lamberts(mu,a,s,c,TOF,lesserTOF,path)
    [alpha,beta,n] = FindVars(mu,a,s,c,lesserTOF,path);

    F = 1/n*((alpha - beta)-(sin(alpha)-sin(beta))) - TOF;
end

