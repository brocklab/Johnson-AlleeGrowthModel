function [ v4] = gen_model_v4(p,tsamp, Ninit, modelcode)
switch modelcode
    case 1 % b-d model
        b= p(1);
        d=p(2);
        timesimdatavec = [];
        % originally were running simulations-- we don't want to do this.
        % instead, take the solved system of ODEs and just use that
   
               for i = 1:length(Ninit)
                N0 = Ninit(i);
                V0 = 0;
                C_init(1)=N0;
                C_init(2)=N0.^2;
                C_init(3)= N0.^3;
                C_init(4)= N0.^4;
                C_init(5) = V0;
                C_init(6)= V0;

                f = @(t,C) [(b-d)*C(1); % dn/dt
                            2*C(2)*(b-d) + C(1)*(b+d);   % dn2/dt
                            3*C(3)*(b-d) + 3*C(2)*(b+d) + C(1)*(b-d); %dn3/dt
                            4*C(4)*(b-d)+ 6*C(3)*(b+d) + 4*C(2)*(b-d) + C(1)*(b+d); % dn4/dt
                            2*C(5)*(b-d) + (b+d)*C(1); %dV2dt
                            4*C(4)*(b-d)+ 6*C(3)*(b+d) + 4*C(2)*(b-d) + C(1)*(b+d)-4*(C(1).^3)*(b-d)*C(1)];%dV4dt
                            %4*C(4)*(b-d)+6*C(3)*(b+d) - 12*C(3)*(b-d)+ 4*C(2)*(b-d) - 12*C(2)*(b+d)-4*C(1)*(b-d)+C(1)*(b+d)]; %dV4dt
                options1 = odeset('Refine',1);  
                options = odeset(options1,'NonNegative',1:6);
                [t,C]=ode45(f, tsamp,C_init, options);
                mu_C= C(:,1);
                n2_C= C(:,2);
                n3_C = C(:,3);
                n4_C = C(:,4);
                v2_C = C(:,5);
                v4_C=C(:,6);

                Cstat(:,1,i) = mu_C;
                Cstat(:,2,i) = n2_C;
                Cstat(:,3,i) = v2_C;
                Cstat(:,4,i) = n3_C;
                Cstat(:,5,i) = n4_C;
                Cstat(:,6,i) = v4_C;
                timesimdatavec = vertcat(timesimdatavec, tsamp');
                end
    
    v4_exp = Cstat(:,6,:);
    v4_exp_long = reshape(v4_exp,size(v4_exp,1)*size(v4_exp,2)*size(v4_exp,3),1);
    v4 = v4_exp_long;
    
    i0 = find(timesimdatavec==0);
    v4(i0)= [];

    
    case 2 % strong Allee model with Allee on birth probability
        b= p(1);
        d=p(2);
        A = p(3);
        timesimdatavec = [];
        % originally were running simulations-- we don't want to do this.
        % instead, take the solved system of ODEs and just use that
   
               for i = 1:length(Ninit)
                N0 = Ninit(i);
                V0 = 0;
                C_init(1)=N0;
                C_init(2)=N0.^2;
                C_init(3)=V0;
                C_init(4)= N0.^3;
                C_init(5)= N0.^4;
                C_init(6) = V0;
            

                f = @(t,C) [((b-d)*C(1)-(b-d)*A); % dn/dt
                 2*C(2).*(b-d) - (2*C(1).*(b-d)*A) + (C(1).*(b+d))-((b-d)*A);   % dn2/dt
                (2*C(2).*(b-d) - (2*C(1).*(b-d)*A) + (C(1).*(b+d))-((b-d)*A)-(2*C(1).*((b-d)*C(1)-(b-d)*A))); % dV/dt
                (3*C(4).*(b-d)) + (3*C(2).*(b+d)-3*C(2)*(b-d)*A) - (3*C(1).*(b-d)*A)-((b-d)*A);%dn3dt
                (4*C(5).*(b-d)) + (6*C(4).*(b+d))-(4*C(4)*(b-d)*A)+(4*C(2).*(b-d))-(6*C(2).*(b-d)*A)+ (C(1).*(b+d))+...
                (4*C(1).*(b-d)*A)-((b-d)*A); %dn4dt
                (4*C(5).*(b-d)) + (6*C(4).*(b+d))-(4*C(4).*(b-d)*A)+(4*C(2).*(b-d))-(6*C(2).*(b-d)*A)+ (C(1).*(b+d))+...
                (4*C(1).*(b-d)*A)-((b-d)*A)-(4.*((C(1)).^3)*(((b-d).*C(1)-(b-d)*A)))];

                    options1 = odeset('Refine',1);  
                    options = odeset(options1,'NonNegative',1:6);
                    [t,C]=ode45(f, tsamp,C_init, options);
                    mu_C= C(:,1);
                    n2_C= C(:,2);
                    v2_C = C(:,3);
                    n3_C = C(:,4);
                    n4_C = C(:,5);
                    v4_C=C(:,6);

                Cstat(:,1,i) = mu_C;
                Cstat(:,2,i) = n2_C;
                Cstat(:,3,i) = v2_C;
                Cstat(:,4,i) = n3_C;
                Cstat(:,5,i) = n4_C;
                Cstat(:,6,i) = v4_C;
                timesimdatavec = vertcat(timesimdatavec, tsamp');
                
                end
    
    v4_exp = Cstat(:,6,:);
    v4_exp_long = reshape(v4_exp,size(v4_exp,1)*size(v4_exp,2)*size(v4_exp,3),1);
    v4 = v4_exp_long;
    
    i0 = find(timesimdatavec==0);
    v4(i0)= [];
    


     case 3 % strong Allee model with Allee on death
        b= p(1);
        d=p(2);
        A = p(3);
        timesimdatavec = [];
        % originally were running simulations-- we don't want to do this.
        % instead, take the solved system of ODEs and just use that
   
               for i = 1:length(Ninit)
                N0 = Ninit(i);
                V0 = 0;
                C_init(1)=N0;
                C_init(2)=N0.^2;
                C_init(3)= V0;
                C_init(4)= N0.^3;
                C_init(5)= N0.^4;
                C_init(6) = V0;
               
                
                f = @(t,C) [((b-d)*C(1)-(b-d)*A); % dn/dt
             2*C(2).*(b-d) - (2*C(1).*(b-d)*A) + (C(1).*(b+d))+((b-d)*A);   % dn2/dt
            (2*C(2).*(b-d) - (2*C(1).*(b-d)*A) + (C(1).*(b+d))+((b-d)*A)-(2.*C(1).*((b-d)*C(1)-(b-d)*A))); % dV/dt
            (3*C(4).*(b-d)) + (3*C(2).*(b+d))-(3*C(2).*(b-d)*A) + (3*C(1).*(b-d)*A) + ((b-d).*C(1))-((b-d)*A);%dn3dt
            (4*C(5).*(b-d)) + (6*C(4).*(b+d))-(4*C(4)*(b-d)*A)+(4*C(2).*(b-d))+(6*C(2).*(b-d)*A)+ (C(1).*(b+d))+...
            (4*C(1).*(b-d)*A)+((b-d)*A); %dn4dt
            (4*C(5).*(b-d)) + (6*C(4).*(b+d))-(4*C(4).*(b-d)*A)+(4*C(2).*(b-d))+(6*C(2).*(b-d)*A)+ (C(1).*(b+d))+...
            (4*C(1).*(b-d)*A)+((b-d)*A)-(4.*((C(1)).^3)*(((b-d).*C(1)-(b-d)*A)))];


                options1 = odeset('Refine',1);  
                options = odeset(options1,'NonNegative',1:6);
                [t,C]=ode45(f, tsamp,C_init, options);
                mu_C= C(:,1);
                n2_C= C(:,2);
                v2_C = C(:,3);
                n3_C = C(:,4);
                n4_C = C(:,5);
                v4_C=C(:,6);

                Cstat(:,1,i) = mu_C;
                Cstat(:,2,i) = n2_C;
                Cstat(:,3,i) = v2_C;
                Cstat(:,4,i) = n3_C;
                Cstat(:,5,i) = n4_C;
                Cstat(:,6,i) = v4_C;
                
                timesimdatavec = vertcat(timesimdatavec, tsamp');
                
                end
    
    v4_exp = Cstat(:,6,:);
    v4_exp_long = reshape(v4_exp,size(v4_exp,1)*size(v4_exp,2)*size(v4_exp,3),1);
    v4 = v4_exp_long;
    
    i0 = find(timesimdatavec==0);
    v4(i0)= [];   
    
    
    
    
    
    
    case 4 % strong Allee model with b-d 
        b= p(1);
        d=p(2);
        A = p(3);
        timesimdatavec = [];
        % originally were running simulations-- we don't want to do this.
        % instead, take the solved system of ODEs and just use that
   
               for i = 1:length(Ninit)
                N0 = Ninit(i);
                V0 = 0;
                C_init(1)=N0;
                C_init(2)=N0.^2;
                C_init(3)= V0;
                C_init(4)= N0.^3;
                C_init(5)= N0.^4;
                C_init(6) = V0;
               
                
                f = @(t,C) [((b-d)*C(1)-(b-d)*A); % dn/dt
                     2*C(2).*(b-d)- 2.*C(1)*(b-d)*A + C(1).*(b+d);   % dn2/dt
                      2*C(2).*(b-d)- 2.*C(1)*(b-d)*A + C(1).*(b+d)- 2.*C(1).*((b-d)*C(1)-(b-d)*A); % dV/dt
                     (3.*C(4).*(b-d)) + (3.*C(2).*(b+d))-(3.*C(2).*(b-d).*A) + (C(1).*(b-d)) - ((b-d)*A);%dn3dt
                    (4.*C(5).*(b-d)) + (6.*C(4).*(b+d))-(4.*C(4).*(b-d).*A)+(4.*C(2).*(b-d))+ (C(1).*(b+d))+...
                     (4.*C(1).*(b-d).*A); %dn4dt
                     (4.*C(5).*(b-d)) + (6.*C(4).*(b+d))-(4.*C(4).*(b-d).*A)+(4.*C(2).*(b-d))+ (C(1).*(b+d))+...
                     (4.*C(1).*(b-d).*A)- (4.*((C(1)).^3).*(((b-d).*C(1)-(b-d).*A)))];

                options1 = odeset('Refine',1);  
                options = odeset(options1,'NonNegative',1:6);
                [t,C]=ode45(f, tsamp,C_init, options);
                mu_C= C(:,1);
                n2_C= C(:,2);
                v2_C = C(:,3);
                n3_C = C(:,4);
                n4_C = C(:,5);
                v4_C=C(:,6);

                Cstat(:,1,i) = mu_C;
                Cstat(:,2,i) = n2_C;
                Cstat(:,3,i) = v2_C;
                Cstat(:,4,i) = n3_C;
                Cstat(:,5,i) = n4_C;
                Cstat(:,6,i) = v4_C;
                
                timesimdatavec = vertcat(timesimdatavec, tsamp');
                
                end
    
    v4_exp = Cstat(:,6,:);
    v4_exp_long = reshape(v4_exp,size(v4_exp,1)*size(v4_exp,2)*size(v4_exp,3),1);
    v4 = v4_exp_long;
    
    i0 = find(timesimdatavec==0);
    v4(i0)= [];
    
    
     case 5 % weak Allee with Allee on birth
        b= p(1);
        d=p(2);
        A = p(3);
        tau = p(4);
        timesimdatavec = [];
        % originally were running simulations-- we don't want to do this.
        % instead, take the solved system of ODEs and just use that
   
               for i = 1:length(Ninit)
                N0 = Ninit(i);
                V0 = 0;
                C_init(1)=N0;
                C_init(2)=N0.^2;
                C_init(3)= V0;
                C_init(4)= N0.^3;
                C_init(5)= N0.^4;
                C_init(6) = V0;
               
               f = @(t,C) [((b-d)*C(1)*(1-((A+tau)./(C(1)+tau))));  % dN/dt
               2.*C(2).*(b-d) + (b+d).*C(1) - 2.*C(2).*(b-d).*((A+tau)./(C(1)+tau))-(b-d).*C(1).*((A+tau)./(C(1)+tau));
               2.*C(2).*(b-d) + (b+d).*C(1) - 2.*C(2).*(b-d).*((A+tau)./(C(1)+tau))-(b-d).*C(1).*((A+tau)./(C(1)+tau))+...
               -2.*C(1).*((b-d)*C(1)*(1-((A+tau)./(C(1)+tau))));
               3.*C(4).*(b-d) - 3.*C(4).*(b-d).*((A+tau)./(C(1)+tau)) + 3.*C(2).*(b+d) - 3.*C(2).*(b-d).*((A+tau))./(C(1)+tau)+...
               + C(1).*(b-d) - C(1).*(b-d).*((A+tau)./(C(1)+tau)); %dn3/dt
               4.*C(5).*(b-d) - 4.*C(5).*(b-d).*((A+tau)./(C(1)+tau)) + 6.*C(4).*(b+d) - 6.*C(4).*(b-d).*((A+tau)./(C(1)+tau)) + 4.*C(2).*(b-d)+...
               -4.*C(2).*(b-d).*((A+tau)./(C(1)+tau)) + C(1).*(b+d) - C(1).*(b-d).*((A+tau)./(C(1)+tau)); %dn4/dt
                4.*C(5).*(b-d) - 4.*C(5).*(b-d).*((A+tau)./(C(1)+tau)) + 6.*C(4).*(b+d) - 6.*C(4).*(b-d).*((A+tau)./(C(1)+tau)) + 4.*C(2).*(b-d)+...
               -4.*C(2).*(b-d).*((A+tau)./(C(1)+tau)) + C(1).*(b+d) - C(1).*(b-d).*((A+tau)./(C(1)+tau))+...
               -4.*((C(1)).^3).*((b-d)*C(1)*(1-((A+tau)./(C(1)+tau))))]; %dV4dt
                
                options1 = odeset('Refine',1);  
                options = odeset(options1,'NonNegative',1:6);
                [t,C]=ode45(f, tsamp,C_init, options);
                mu_C= C(:,1);
                n2_C= C(:,2);
                v2_C = C(:,3);
                n3_C = C(:,4);
                n4_C = C(:,5);
                v4_C=C(:,6);

                Cstat(:,1,i) = mu_C;
                Cstat(:,2,i) = n2_C;
                Cstat(:,3,i) = v2_C;
                Cstat(:,4,i) = n3_C;
                Cstat(:,5,i) = n4_C;
                Cstat(:,6,i) = v4_C;
                
                timesimdatavec = vertcat(timesimdatavec, tsamp');
                
                end
    
    v4_exp = Cstat(:,6,:);
    v4_exp_long = reshape(v4_exp,size(v4_exp,1)*size(v4_exp,2)*size(v4_exp,3),1);
    v4 = v4_exp_long;
    
    i0 = find(timesimdatavec==0);
    v4(i0)= [];
    
    
    
    
    
    case 6 % weak Allee with Allee on death
        b= p(1);
        d=p(2);
        A = p(3);
        tau = p(4);
        timesimdatavec = [];
        % originally were running simulations-- we don't want to do this.
        % instead, take the solved system of ODEs and just use that
   
               for i = 1:length(Ninit)
                N0 = Ninit(i);
                V0 = 0;
                C_init(1)=N0;
                C_init(2)=N0.^2;
                C_init(3)= V0;
                C_init(4)= N0.^3;
                C_init(5)= N0.^4;
                C_init(6) = V0;
               
               f = @(t,C) [((b-d)*C(1)*(1-((A+tau)./(C(1)+tau))));  % dN/dt
               2.*C(2).*(b-d) + (b+d).*C(1) - 2.*(C(1).^2).*(b-d).*((A+tau)./(C(1)+tau))+(b-d).*C(1).*((A+tau)./(C(1)+tau));
               2.*C(2).*(b-d) + (b+d).*C(1) - 2.*C((2)).*(b-d).*((A+tau)./(C(1)+tau))+(b-d).*C(1).*((A+tau)./(C(1)+tau))+...
               -2.*C(1).*((b-d)*C(1)*(1-((A+tau)./(C(1)+tau))));
               3.*C(4).*(b-d) - 3.*(C(2).*C(1)).*(b-d).*((A+tau)./(C(1)+tau)) + 3.*C(2).*(b+d) + 3.*(C(1).^2).*(b-d).*((A+tau))./(C(1)+tau)+...
               + C(1).*(b-d) - C(1).*(b-d).*((A+tau)./(C(1)+tau)); %dn3/dt
               % fix n4-V4!!!
                4.*C(5).*(b-d) - 4.*(C(4).*C(1)).*(b-d).*((A+tau)./(C(1)+tau)) + 6.*C(4).*(b+d) + 6.*(C(2).*C(1))*(b-d).*((A+tau)./(C(1)+tau)) + 4.*C(2).*(b-d)+...
                -4.*(C(1).^2).*(b-d).*((A+tau)./(C(1)+tau)) + C(1).*(b+d) + C(1).*(b-d).*((A+tau)./(C(1)+tau)); %dn4/dt
                 4.*C(5).*(b-d) - 4.*C(5).*(b-d).*((A+tau)./(C(1)+tau)) + 6.*C(4).*(b+d) + 6.*C(4).*(b-d).*((A+tau)./(C(1)+tau)) + 4.*C(2).*(b-d)+...
                -4.*C(2).*(b-d).*((A+tau)./(C(1)+tau)) + C(1).*(b+d) + C(1).*(b-d).*((A+tau)./(C(1)+tau))+...
                -4.*((C(1)).^3).*((b-d)*C(1)*(1-((A+tau)./(C(1)+tau))))]; %dV4dt
                options1 = odeset('Refine',1);  
                options = odeset(options1,'NonNegative',1:6);
                [t,C]=ode45(f, tsamp,C_init, options);
                mu_C= C(:,1);
                n2_C= C(:,2);
                v2_C = C(:,3);
                n3_C = C(:,4);
                n4_C = C(:,5);
                v4_C=C(:,6);

                Cstat(:,1,i) = mu_C;
                Cstat(:,2,i) = n2_C;
                Cstat(:,3,i) = v2_C;
                Cstat(:,4,i) = n3_C;
                Cstat(:,5,i) = n4_C;
                Cstat(:,6,i) = v4_C;
                
                timesimdatavec = vertcat(timesimdatavec, tsamp');
                
                end
    
    v4_exp = Cstat(:,6,:);
    v4_exp_long = reshape(v4_exp,size(v4_exp,1)*size(v4_exp,2)*size(v4_exp,3),1);
    v4 = v4_exp_long;
    
    i0 = find(timesimdatavec==0);
    v4(i0)= [];
    
    
    
    case 7 % weak Allee with Allee on birth & death
        b= p(1);
        d=p(2);
        A = p(3);
        tau = p(4);
        timesimdatavec = [];
        % originally were running simulations-- we don't want to do this.
        % instead, take the solved system of ODEs and just use that
   
               for i = 1:length(Ninit)
                N0 = Ninit(i);
                V0 = 0;
                C_init(1)=N0;
                C_init(2)=N0.^2;
                C_init(3)= V0;
                C_init(4)= N0.^3;
                C_init(5)= N0.^4;
                C_init(6) = V0;
               
               f = @(t,C) [((b-d)*C(1)*(1-((A+tau)./(C(1)+tau))));  % dN/dt
           2.*C(2).*(b-d) + (b+d).*C(1) - 2.*(C(1).^2).*(b-d).*((A+tau)./(C(1)+tau)); %dN2/dt
           2.*C(2).*(b-d) + (b+d).*C(1) - 2.*C(2).*(b-d).*((A+tau)./(C(1)+tau))-2.*C(1).*((b-d).*C(1).*(1-((A+tau)./(C(1)+tau)))); %dV/dt
           3.*C(4).*(b-d) - 3.*(C(2)*C(1)).*(b-d).*((A+tau)./(C(1)+tau)) + 3.*C(2).*(b+d) + C(1).*(b-d)+...
           -C(1).*(b-d).*((A+tau)./(C(1)+tau)); %dn3/dt
           4.*C(5).*(b-d) - 4.*(C(4)*C(1)).*(b-d).*((A+tau)./(C(1)+tau)) + 6.*C(4).*(b+d) + 4.*C(2).*(b-d)+...
           -4.*(C(1).^2).*(b-d).*((A+tau)./(C(1)+tau)) + C(1).*(b+d); %dn4/dt
            4.*C(5).*(b-d) - 4.*C(5).*(b-d).*((A+tau)./(C(1)+tau)) + 6.*C(4).*(b+d) + 4.*C(2).*(b-d)+...
           -4.*C(2).*(b-d).*((A+tau)./(C(1)+tau)) + C(1).*(b+d)+...
           -4.*((C(1)).^3).*((b-d)*C(1)*(1-((A+tau)./(C(1)+tau))))]; %dV4/dt
                
            
                options1 = odeset('Refine',1);  
                options = odeset(options1,'NonNegative',1:6);
                [t,C]=ode45(f, tsamp,C_init, options);
                mu_C= C(:,1);
                n2_C= C(:,2);
                v2_C = C(:,3);
                n3_C = C(:,4);
                n4_C = C(:,5);
                v4_C=C(:,6);

                Cstat(:,1,i) = mu_C;
                Cstat(:,2,i) = n2_C;
                Cstat(:,3,i) = v2_C;
                Cstat(:,4,i) = n3_C;
                Cstat(:,5,i) = n4_C;
                Cstat(:,6,i) = v4_C;
                
                timesimdatavec = vertcat(timesimdatavec, tsamp');
                
                end
    
    v4_exp = Cstat(:,6,:);
    v4_exp_long = reshape(v4_exp,size(v4_exp,1)*size(v4_exp,2)*size(v4_exp,3),1);
    v4 = v4_exp_long;
    
    i0 = find(timesimdatavec==0);
    v4(i0)= [];
end
end