function[V_m, I_Na, I_K, g_Na, g_K, n, m, h, I_ion, I_L, I_c, I_m, I_i, ...
    t, x] = Parte2HH(delta_t, t_final, delta_x, comprimento, I_p, t_0, t_corrente , alinea_e)

% Método de Euler para resolver as equações de Hodgkin-Huxley
% Propagação dum Potencial de Ação ao longo de um Axónio Não-Mielinizado
%
% Entradas:
%   Passo temporal e espacial (delta_t, delta_x) 
%   Duração da simulação (t_final)
%   Vetor de estímulo (I_p) 
%   Duração do estimulo (t_corrente)
%   Instante inicial do estímulo (t_0)
%   Comprimento do Axónio (comprimento)
%   Opção para a Alínea E (alinea_e)
%
% Saídas: 
%   Gating probabilities (n,m,h) 
%   Condutâncias (g_Na,g_K)
%   Potencial de membrana (V_m)
%   Correntes (I_Na,I_K,I_Ion,I_L,I_c,I_m,I_i) 
%   Vetor temporal e espacial (t, x)

    parametros;
    
    
    % Equações para alfas e betas
   
    alfa_n = @(v_m) 0.01*(10-v_m)./(exp((10-v_m)./10)-1);
    alfa_m = @(v_m) 0.1*(25-v_m)./(exp((25-v_m)./10)-1);
    beta_h = @(v_m) 1./(exp((30-v_m)/10)+1);
    beta_n = @(v_m) 0.125*exp(-v_m/80);
    beta_m = @(v_m) 4*exp(-v_m/18);
    alfa_h = @(v_m) 0.07*exp(-v_m/20);

    

    for i = 1:t_pontos-1                                                   % Ciclo for temporal 
        
        for j = 1:x_pontos                                                 % Ciclo for espacial
            
            %Definir corrente transmembranar em função da posição espacial
            %e tendo em conta as condições de fronteira
            
            if j == 1                                                      % Extremidade inicial
                
                I_m(i,j) = 100/(2*pi*a*(r_i+r_e))*((V_m(i,j+1)- ...
                    V_m(i,j))/delta_x^2 - r_e*I_p(i,j)); 
                   
            elseif j == x_pontos && i < (t_0 + t_corrente)/delta_t         % Extremidade final durante simulação da corrente de estimulação
          
                I_m(i,x_pontos) = 100/(2*pi*a*(r_i+r_e))*((V_m(i,x_pontos)- ...
                    V_m(i,x_pontos-1))/delta_x^2 - r_e*I_p(i,x_pontos));
                
            end
            
            
            if j ~=1 && j < x_pontos 
                
                I_m(i,j) = 100/(2*pi*a*(r_i+r_e))*((V_m(i,j+1)- ...
                    2*V_m(i,j)+V_m(i,j-1))/delta_x^2 - r_e*I_p(i,j));
                
            end
            
            % Condutâncias dos canais iónicos
            g_K(i,j) = gi_K*n(i,j)^4;
            g_Na(i,j) = gi_Na*m(i,j)^3*h(i,j);
            
            % Correntes iónicas e de leakage
            I_K(i,j) = g_K(i,j)*(V_m(i,j) -  E_K);
            I_Na(i,j) = g_Na(i,j)*(V_m(i,j) -  E_Na);
            I_L(i,j) = g_L*(V_m(i,j) - E_L);
            I_ion(i,j) = I_Na(i,j) + I_K(i,j) + I_L(i,j); 
            
            % Potencial de Membrana
            I_c(i,j) = (I_m(i,j) - I_ion(i,j));
            dV = (I_c(i,j))*delta_t/C_m;
            v_m(i+1,j) = v_m(i,j) + dV;
            V_m(i+1,j) = v_m(i+1,j) + V_r;
            I_i(i,j+1) = I_i(i,j) - 2*pi*a*I_m(i,j)*delta_x;
          
           
            % Visualizar a resposta passiva, no caso da alinea e) ser
            % selecionada
            if alinea_e == 1
            
                a_n = 0;
                b_n = 0;
                a_m = 0;
                b_m = 0;
                a_h = 0;
                b_h = 0;
            
            else
                % Definir singularidades de alfa m e alfa n
                if v_m(i+1,j) == 10 || v_m(i+1,j) == 25
                
                    a_m = 1;
                    a_n = 0.1;
                
                else
                    %Avançar os valores de alfas e betas
                    a_n = alfa_n(v_m(i+1,j));    
                    a_m = alfa_m(v_m(i+1,j));    
                    b_h = beta_h(v_m(i+1,j));    
                    b_n = beta_n(v_m(i+1,j));    
                    b_m = beta_m(v_m(i+1,j));    
                    a_h = alfa_h(v_m(i+1,j));    
                
                end
            end
         
            % Avançar os valores para n, m e h
            n(i+1,j) = n(i,j) + (a_n*(1-n(i,j))-b_n*(n(i,j)))*delta_t;
            m(i+1,j) = m(i,j) + (a_m*(1-m(i,j))-b_m*(m(i,j)))*delta_t;
            h(i+1,j) = h(i,j) + (a_h*(1-h(i,j))-b_h*(h(i,j)))*delta_t;
            
            if i == t_pontos-1 && j < x_pontos
               
                x(j+1) = x(j) + delta_x;                                   %Avançar o vetor espacial caso o vetor temporal ja tenha chegado ao ultimo elemento mas o vetor espacial não
            
            end

        end     
        
        t(i+1) = t(i) + delta_t;                                           %Avançar vetor temporal   

    end
end
