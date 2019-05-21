function [profiles] = profile_likelihood(params, tsamp, N0, N,mudatavec, vardatavec, modelcode, factor, numpoints)

switch modelcode
    case 1 % birth-death model
        num_params = 2;
        profile = [];


for k = 1:num_params
profindex = k; % PROFILE b PARAMETER
    % Profile
    if k ==2 % since d is so small
    factoruse = 8*factor;%8e15*factor;%8*factor;%3e15%4*factor;
    else 
        factoruse = factor;
    end% increase factor for d
    profrangeDown = linspace((params(profindex)), (params(profindex)*(1-factoruse)),numpoints)'; 
    profrangeUp = linspace((params(profindex)), (params(profindex)*(1+factoruse)),numpoints)';
    % split into up and down so we can use last fitted value as starting value for next run
    profrange = [profrangeDown;profrangeUp];
    profrange = sort(profrange);
    currfvals = [];
    currparams = [];
    currflags = [];
    paramstemp = [];
    profile = [];
    for m = 1:length(profrange)
        [m] %track progress
        currp = profrange(m);
        % need to write a fitting function that only allows d to be
        % searched on
        [fvaltemp, paramstemp] = ML_fitnegLL(params, tsamp, N0, N, mudatavec, vardatavec, modelcode, profindex, currp);
        % fminsearch will out put the values of dguess that give the lowest
        % objfun_donly

        currfvals = [currfvals; fvaltemp];
        currparams = [currparams; [profrange(m),paramstemp]]; %storing the profiled value too, so the output parameter values are easy to run the model with
    end
    
    profile = horzcat(profrange, real(currfvals));
 


    profiles(:,:,profindex) = profile;

% 1st column is the parameter values that are "profiled"
% 2nd column is the negLL corresponding to each "profiled" parameter
% each dimension is for a parameter




 
end





case 2 % strong Allee on birth
        num_params = 3;
        profile = [];


for k = 1:num_params
profindex = k; % PROFILE b PARAMETER
    % Profile
    if k ==2 % since d is so small
    factor = 5*k*factor;
    end% increase factor for d
    if k ==3
        factor = factor;
    end
    profrangeDown = linspace((params(profindex)), (params(profindex)*(1-factor)),numpoints)'; 
    profrangeUp = linspace((params(profindex)), (params(profindex)*(1+factor)),numpoints)';
    % split into up and down so we can use last fitted value as starting value for next run
    profrange = [profrangeDown;profrangeUp];
    profrange = sort(profrange);
    currfvals = [];
    currparams = [];
    currflags = [];
    paramstemp = [];
    profile = [];
    for m = 1:length(profrange)
        [m] %track progress
        currp = profrange(m);
        % need to write a fitting function that only allows d to be
        % searched on
        [fvaltemp, paramstemp] = ML_fitnegLL(params, tsamp, N0, N, mudatavec, vardatavec, modelcode, profindex, currp);
        % fminsearch will out put the values of dguess that give the lowest
        % objfun_donly

        currfvals = [currfvals; fvaltemp];
        currparams = [currparams; [profrange(m),paramstemp]]; %storing the profiled value too, so the output parameter values are easy to run the model with
    end
    
    profile = horzcat(profrange, real(currfvals));
 


    profiles(:,:,profindex) = profile;

% 1st column is the parameter values that are "profiled"
% 2nd column is the negLL corresponding to each "profiled" parameter
% each dimension is for a parameter

end



case 3 % strong Allee on death
        num_params = 3;
        profile = [];


        for k = 1:num_params
        profindex = k; % PROFILE b PARAMETER
            % Profile
            if k ==1
                factor = 0.1*factor;
            end
            if k ==2 % since d is so small
            factor = 1;
            end% increase factor for d
            if k ==3
                factor = 0.1*factor;
            end
            profrangeDown = linspace((params(profindex)), (params(profindex)*(1-factor)),numpoints)'; 
            profrangeUp = linspace((params(profindex)), (params(profindex)*(1+factor)),numpoints)';
            % split into up and down so we can use last fitted value as starting value for next run
            profrange = [profrangeDown;profrangeUp];
            profrange = sort(profrange);
            currfvals = [];
            currparams = [];
            currflags = [];
            paramstemp = [];
            profile = [];
            for m = 1:length(profrange)
                [m] %track progress
                currp = profrange(m);
                % need to write a fitting function that only allows d to be
                % searched on
                [fvaltemp, paramstemp] = ML_fitnegLL(params, tsamp, N0, N, mudatavec, vardatavec, modelcode, profindex, currp);
                % fminsearch will out put the values of dguess that give the lowest
                % objfun_donly

                currfvals = [currfvals; fvaltemp];
                currparams = [currparams; [profrange(m),paramstemp]]; %storing the profiled value too, so the output parameter values are easy to run the model with
            end

            profile = horzcat(profrange, real(currfvals));



            profiles(:,:,profindex) = profile;

        % 1st column is the parameter values that are "profiled"
        % 2nd column is the negLL corresponding to each "profiled" parameter
        % each dimension is for a parameter
        end



case 4 % strong Allee on birth & death
        num_params = 3;
        profile = [];


        for k = 1:num_params
        profindex = k; % PROFILE b PARAMETER
            % Profile
            if k ==1
                factor = 0.08*factor;
            end
            if k ==2 % since d is so small
            factor = 1;
            end% increase factor for d
            if k ==3
                factor = 0.4*factor;
            end
            profrangeDown = linspace((params(profindex)), (params(profindex)*(1-factor)),numpoints)'; 
            profrangeUp = linspace((params(profindex)), (params(profindex)*(1+factor)),numpoints)';
            % split into up and down so we can use last fitted value as starting value for next run
            profrange = [profrangeDown;profrangeUp];
            profrange = sort(profrange);
            currfvals = [];
            currparams = [];
            currflags = [];
            paramstemp = [];
            profile = [];
            for m = 1:length(profrange)
                [m] %track progress
                currp = profrange(m);
                % need to write a fitting function that only allows d to be
                % searched on
                [fvaltemp, paramstemp] = ML_fitnegLL(params, tsamp, N0, N, mudatavec, vardatavec, modelcode, profindex, currp);
                % fminsearch will out put the values of dguess that give the lowest
                % objfun_donly

                currfvals = [currfvals; fvaltemp];
                currparams = [currparams; [profrange(m),paramstemp]]; %storing the profiled value too, so the output parameter values are easy to run the model with
            end

            profile = horzcat(profrange, real(currfvals));



            profiles(:,:,profindex) = profile;

        % 1st column is the parameter values that are "profiled"
        % 2nd column is the negLL corresponding to each "profiled" parameter
        % each dimension is for a parameter
        end

case 5 % weak/strong Allee on birth
        num_params = 4;
        profile = [];


        for k = 1:num_params
        profindex = k; % PROFILE b PARAMETER
            % Profile
            if k ==2 % since d is so small
            factor = 8;
            end% increase factor for d
           if k ==3 
                factor = 0.1;
            end
            if k==4
                factor = 0.1;
            end
            profrangeDown = linspace((params(profindex)), (params(profindex)*(1-factor)),numpoints)'; 
            profrangeUp = linspace((params(profindex)), (params(profindex)*(1+factor)),numpoints)';
            % split into up and down so we can use last fitted value as starting value for next run
            profrange = [profrangeDown;profrangeUp];
            profrange = sort(profrange);
            currfvals = [];
            currparams = [];
            currflags = [];
            paramstemp = [];
            profile = [];
            for m = 1:length(profrange)
                [m] %track progress
                currp = profrange(m);
                % need to write a fitting function that only allows d to be
                % searched on
                [fvaltemp, paramstemp] = ML_fitnegLL(params, tsamp, N0, N, mudatavec, vardatavec, modelcode, profindex, currp);
                % fminsearch will out put the values of dguess that give the lowest
                % objfun_donly

                currfvals = [currfvals; fvaltemp];
                currparams = [currparams; [profrange(m),paramstemp]]; %storing the profiled value too, so the output parameter values are easy to run the model with
            end

            profile = horzcat(profrange, real(currfvals));



            profiles(:,:,profindex) = profile;

        % 1st column is the parameter values that are "profiled"
        % 2nd column is the negLL corresponding to each "profiled" parameter
        % each dimension is for a parameter
        end
        
        
case 6 % weak/strong Allee on death
        num_params = 4;
        profile = [];


        for k = 1:num_params
        profindex = k; % PROFILE b PARAMETER
            % Profile
            if k ==2 % since d is so small
            factor = 1;
            end% increase factor for d
            if k ==3 
                factor = 1.5;
            end
            if k==4
                factor = 2.5;
            end
            
            if k==2 % special case bc d=0)
                profrangeUp= linspace((params(profindex)), ((params(profindex)+1e-4)*(1+factor)),numpoints)';
            else
                profrangeUp = linspace((params(profindex)), (params(profindex)*(1+factor)),numpoints)';
            end
            profrangeDown = linspace((params(profindex)), (params(profindex)*(1-factor)),numpoints)'; 
           
            % split into up and down so we can use last fitted value as starting value for next run
            profrange = [profrangeDown;profrangeUp];
            profrange = sort(profrange);
            currfvals = [];
            currparams = [];
            currflags = [];
            paramstemp = [];
            profile = [];
            for m = 1:length(profrange)
                [m] %track progress
                currp = profrange(m);
                % need to write a fitting function that only allows d to be
                % searched on
                [fvaltemp, paramstemp] = ML_fitnegLL(params, tsamp, N0, N, mudatavec, vardatavec, modelcode, profindex, currp);
                % fminsearch will out put the values of dguess that give the lowest
                % objfun_donly

                currfvals = [currfvals; fvaltemp];
                currparams = [currparams; [profrange(m),paramstemp]]; %storing the profiled value too, so the output parameter values are easy to run the model with
            end

            profile = horzcat(profrange, real(currfvals));



            profiles(:,:,profindex) = profile;

        % 1st column is the parameter values that are "profiled"
        % 2nd column is the negLL corresponding to each "profiled" parameter
        % each dimension is for a parameter
        end
        

        case 7 % weak/strong Allee on birth & death
        num_params = 4;
        profile = [];


        for k = 1:num_params
        profindex = k; % PROFILE b PARAMETER
            % Profile
            if k ==2 % since d is so small
            factor = 1;
            end% increase factor for d
            if k ==3 
                factor = 1.5;
            end
            if k==4
                factor = 2.5;
            end
            
            if k==2 % special case bc d=0)
                profrangeUp= linspace((params(profindex)), ((params(profindex)+1e-4)*(1+factor)),numpoints)';
            else
                profrangeUp = linspace((params(profindex)), (params(profindex)*(1+factor)),numpoints)';
            end
            profrangeDown = linspace((params(profindex)), (params(profindex)*(1-factor)),numpoints)'; 
           
            % split into up and down so we can use last fitted value as starting value for next run
            profrange = [profrangeDown;profrangeUp];
            profrange = sort(profrange);
            currfvals = [];
            currparams = [];
            currflags = [];
            paramstemp = [];
            profile = [];
            for m = 1:length(profrange)
                [m] %track progress
                currp = profrange(m);
                % need to write a fitting function that only allows d to be
                % searched on
                [fvaltemp, paramstemp] = ML_fitnegLL(params, tsamp, N0, N, mudatavec, vardatavec, modelcode, profindex, currp);
                % fminsearch will out put the values of dguess that give the lowest
                % objfun_donly

                currfvals = [currfvals; fvaltemp];
                currparams = [currparams; [profrange(m),paramstemp]]; %storing the profiled value too, so the output parameter values are easy to run the model with
            end

            profile = horzcat(profrange, real(currfvals));



            profiles(:,:,profindex) = profile;

        % 1st column is the parameter values that are "profiled"
        % 2nd column is the negLL corresponding to each "profiled" parameter
        % each dimension is for a parameter
        end
        
        
        

end
        
        
end
