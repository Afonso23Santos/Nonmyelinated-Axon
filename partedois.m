%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               BIOELETRICIDADE - PARTE 2                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%Engenharia Biom�dica e Biof�sica, 2019
%Afonso Pourfarzaneh, 48249
%�ngelo Nunes, 48254
%Maria Lopes, 49941

clear all
clear vars
clc
%% Al�nea a) Definir mesh ratio.
delta_t = input('Passo temporal (ms):  ');                                 %Passo temporal
    
if isempty(delta_t)
   disp('Valor atribu�do: 0.02 ms.')
   delta_t = 0.02;
end      
    
delta_x = input('Passo espacial (cm):  ');                                 %Passo espacial
    
if isempty(delta_x)
   disp('Valor atribu�do: 0.1 cm.')
   delta_x = 0.1;
end

parametros;

mesh_ratio = delta_t/(c_m*r_i*delta_x^2);                                  %Retirado do Plonsey

fprintf('Valor da mesh ratio �: %d \n', mesh_ratio)
%% Al�nea c) e d) Gr�ficos para os nodos 0 e 600, e 300.

%Par�metros pedidos ao utilizador 

t_corrente = input('Dura��o da corrente de estimula��o (ms): ');
if isempty(t_corrente)
    disp('Valor atribu�do: 0.1 ms.')
    t_corrente = 0.1;
end

t_0 = input('Instante a que a corrente � aplicado (ms): ');
if isempty(t_0)
    disp('Valor atribu�do: 5 ms.')
    t_0 = 5;
end


I_s = 3.95;                                                                %uA/cm - Amplitude da corrente de est�mulo
                                                                 
I_p = zeros(t_pontos,x_pontos);                                            %Vetor da corrente de estimulo

loc_impulso = ...
    input("Impulso nas extremidades (op��o 0) ou no centro(op��o 1): ");                                           

if loc_impulso == 0                                                        % Ciclo if para preencher o vetor da corrente de estimulo com a amplitude definida
    
    %Para o est�mulo nas extremidades - nodos 0 e 600
    I_p(t_0/delta_t:(t_0 + t_corrente)/delta_t,1) = -I_s; 
    I_p(t_0/delta_t:(t_0 + t_corrente)/delta_t,x_pontos) = -I_s;

elseif loc_impulso == 1
    
    %Para o estimulo no centro - nodo 300
    I_p(t_0/delta_t:(t_0 + t_corrente)/delta_t,x_pontos/2) = -I_s;
end

%% Al�nea e) - Resposta passiva

alinea_e = ...
    input('Observar a resposta de uma fibra passiva?(sim = 1,n�o = 0): ');

[V_m, I_Na, I_K, g_Na, g_K, n, m, h, I_ion, I_L, I_c, I_m, I_i, t, x] = ...
    Parte2HH(delta_t, t_final, delta_x, comprimento, I_p, t_0,...
    t_corrente, alinea_e);
%% Al�nea b) - C�lculo do limiar da corrente de estimula��o

alinea_b = ...
    input('Calcular o limiar da corrente de estimula��o?(sim = 1, n�o = 0) ');

if alinea_b == 1
    
    % Definir o vetor da corrente de estimula��o
    
    vetor_corrente = I_s-0.05:0.01:I_s+0.05;
    i_pontos = size(vetor_corrente);
    V_maximo = zeros(1,i_pontos(2));
    
    for i = 1:i_pontos(2)
        
        corrente_estimulo = zeros(t_pontos,x_pontos);                      %Vetor de est�mulo a ser preenchido pelos valores da corrente de est�mulo

        corrente_estimulo(t_0/delta_t:(t_0 + ...
            t_corrente)/delta_t,x_pontos/2) = - vetor_corrente(i);         %1 ms de impulso

        V_m = Parte2HH(delta_t, t_final, delta_x, comprimento,...
            corrente_estimulo, t_0, t_corrente, alinea_e);                 %Calcular potencial da membrana gerado
        V_maximo(i) = max(max(V_m));
    
    end
    
    figure(1)
    vetor_corrente = vetor_corrente/(comprimento*delta_x);                 %Corrente por unidade de comprimento
    plot(vetor_corrente,V_maximo,'>-')
    xlabel('Corrente de Estimulo (mA/cm)');
    ylabel('Potencial membranar (mV)');
    title(['Limiar para o Potencial de A��o', ...
        ' (para x = 15cm, dt = 0.02 ms, t_{dur} = 0.1 ms)'])
    
    disp(vetor_corrente(find(V_maximo >= 0, 1)))                           %Indica o limiar calculado
    
end
%% Gr�ficos

% Caso o impulso seja nas extremidades e a alinea b) n�o seja escolhida ou
% caso a alinea e) seja escolhida mas a alinea b) n�o

if loc_impulso == 0 && alinea_b  == 0 || alinea_e == 1 && alinea_b == 0    

    figure(3)
    subplot(121)
    s4 = surf(x,t,V_m,'FaceAlpha',0.5);
    s4.EdgeColor = 'none'; 
    xlabel('x (cm)'), ylabel('t (ms)'), zlabel('V_m (mV)')
    title('Potencial Transmembranar ao longo da posi��o e do tempo')
    colormap hot
    
    subplot(122)
    imagesc([x(1) x(end)],[t(1) t(end)],V_m)
    xlabel('x (cm)')
    ylabel('t (ms)')
    title(['Potencial Transmembranar ao longo da posi��o e do tempo', ...
        ' (Visto de cima)'])
    colormap hot
    
    figure (2)
    title('Fluxo de Correntes I�nicas e Transmembranar (t = 15 ms)')
    yyaxis('left')
    plot(x(1:x_pontos/5), I_Na(15/delta_t,1:x_pontos/5)/1000)
    xlabel('x (cm)')
    ylabel('I_{i�nicas} (mA)')
    hold on
    plot(x(1:x_pontos/5), I_K(15/delta_t,1:x_pontos/5)/1000)
    yyaxis('right')
    plot(x(1:x_pontos/5),I_m(15/delta_t,1:x_pontos/5)/1000)
    legend('I_{Na}', 'I_{K}', 'I_m')
    ylabel('I_{m} (mA)')
    
    figure(1)
    yyaxis('left')
    plot(x(1:x_pontos/5),I_i(15/delta_t,1:x_pontos/5))
    ylabel('I (mA)')
    hold on
    yyaxis('right')
    plot(x(1:x_pontos/5),V_m(15/delta_t,1:x_pontos/5))
    axis([x(1) x(x_pontos/5) -80 60])
    legend('I_i','V_m')
    xlabel('x (cm)')
    title('Potencial Transmembranar e Corrente Axial Intracelular (t = 15 ms)')
    ylabel('V (mV)')
    
end

% Caso o impulso seja no centro

if loc_impulso == 1 && alinea_b == 0 || alinea_e == 1 && alinea_b == 0
    
    figure(5)
    title('Fluxo de Correntes I�nicas e Transmembranar (x = 16 cm)')
    yyaxis('left')
    plot(t, I_Na(:,x_pontos/2+1/delta_x)/1000)
    hold on
    plot(t, I_K(:,x_pontos/2+1/delta_x)/1000)
    ylabel('I_{i�nicas} (mA)')
    yyaxis('right')
    plot(t, I_m(:,x_pontos/2+1/delta_x)/1000)
    legend('I_{Na}', 'I_{K}', 'I_m')
    xlabel('t')
    ylabel('I_{m} (mA)') 
    
    figure(4)
    yyaxis('right')
    plot(t, V_m(:,x_pontos/2+1/delta_x))
    axis([t(1) t(end) -80 50])
    ylabel('V (mV)')
    hold on
    yyaxis('left')
    plot(t, I_i(:,x_pontos/2+1/delta_x)/1000)
    ylabel('I (mA)')
    legend('I_i','V_m')
    xlabel('t')
    title(['Potencial Transmembranar e Corrente Axial Intracelular', ...
        ' (x = 16 cm)'])
    
end
