%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               BIOELETRICIDADE - PARTE 1                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%Engenharia Biomédica e Biofísica, 2019
%Afonso Pourfarzaneh, 48249
%Ângelo Nunes, 48254
%Maria Lopes, 49941
clearvars
close all
clc

%% Parâmetros

% Definir as constantes
V_r = -60;                                                                 %(mV) - Potencial de repouso
C_m = 1.0;                                                                 %(uF/cm^2) - Capacidade da membrana
 
% Definir potenciais
E_K = -72.10;                                                              %(mV) - Potencial de Nernst para o K
E_Na = 52.4;                                                               %(mV) - Potencial de Nernst para o Na
E_L = -49.187;                                                             %(mV) - Leakage equilibrium potential                            

%Definir condutâncias
gi_Na = 120.0;                                                             %(mS/cm^2) - Condutância máxima de Na
gi_K = 36.0;                                                               %(mS/cm^2) - Condutância máxima de K
g_L = 0.3;                                                                 %(mS/cm^2) - Condutância de leakage

%% Recolha de inputs do utilizador para definir: Is, t0_Is, duration_Is, t_total, delta_t

Is = input('Indicar intensidade da corrente aplicada (uA/cm^2): ');
if isempty(Is)  
    disp(['Vazio: 10 ' char(956) 'A/cm^2.'])
    Is = 10;
end

t0_Is = input('Indicar tempo inicial (ms): ');
if isempty(t0_Is)
    disp('Vazio: 5 ms.')
    t0_Is = 5;
end

duration_Is = input('Indicar duração da corrente aplicada (ms): ');
if isempty(duration_Is)
    disp('Vazio: 1 ms.')
    duration_Is = 1;
end

delta_t = input('Indicar passo temporal (ms): '); 
if isempty(delta_t)
    disp('Vazio: 0.01 ms.')
    delta_t = 0.01;
end

t_total = input('Indicar duração total da simulação(ms)- mínimo de 30 ms: ');
if isempty(t_total)
    disp('O valor utilizado será de 30 ms.')
    t_total = 30;
end

while t_total < 30
    disp('Minímo de 30 ms!');
    t_total = input('Indicar duração total da simulação(ms)- mínimo de 30 ms: ');
end

%% Definir vetor de tempo

Npoints = t_total/delta_t;                                                 %número de pontos
time_vector = zeros(1,Npoints);                                            %ms
length_time = length (time_vector); 
n_0 = t0_Is/delta_t;                                                       %número de pontos antes da aplicação da corrente

% Inicialização de vector Im - corrente da membrana - atenção que esta varia
% consoante Is
Im = zeros(1,Npoints);
%%  Inicialiazação de valores inicias de vetores Vm, vm, n, m, h, g_Na, g_K,correntes

vm = time_vector;                                                          %Potencial de membrana reduzido
Vm = vm + V_r;                                                             %Potencial de membrana
n = time_vector; n(1) = 0.31768;                                           %Gating probability n 
m = time_vector; m(1) = 0.05293;                                           %Gating probability m
h = time_vector; h(1) = 0.59612;                                           %Gating probability h
g_Na = time_vector;                                                        %Condutância dos canais de Na
g_K = time_vector;                                                         %Condutância dos canais de K
I_K = time_vector;                                                         %Corrente dos canais de K
I_Na = time_vector;                                                        %Corrente dos canais de Na                                  
I_L = time_vector;                                                         %Corrente de leakage
I_ionica = time_vector;                                                    %Corrente Iónica
Ic = time_vector;                                                          %Corrente Capacitiva 

%% Alínea E
% Período Refratário

num_impulsos = input('Quantos impulsos pretende: ');
if isempty(num_impulsos)
    num_impulsos = 1;
end

intervalo = 0;

if num_impulsos > 1
    intervalo = input('Qual o intervalo temporal entre impulsos (ms): ');
end

maxV = Is;

for i = 1:num_impulsos
    
    Im((intervalo+duration_Is)*(i-1)/delta_t+n_0:((intervalo+duration_Is)*...
        (i-1))/delta_t+duration_Is/delta_t+n_0) = Is*ones(1,duration_Is/delta_t+1);

end
%% Alínea C
% Gráficos

% Método de Euler
    for i = 1:Npoints-1
        
        g_K(i) = gi_K*n(i)^4;
        g_Na(i) = gi_Na*m(i)^3*h(i);
        I_K(i) = g_K(i)*(Vm(i) -  E_K);
        I_Na(i) = g_Na(i)*(Vm(i) -  E_Na);
        I_L(i) = g_L*(Vm(i) - E_L);
        I_ionica(i) = I_Na(i) + I_K(i) + I_L(i); 
        Ic(i) = C_m*(Im(i) - I_ionica(i));
        vm(i+1) = vm(i) + (Ic(i))*delta_t/C_m;
        Vm(i+1) = vm(i+1) + V_r;
        
        if vm(i+1) == 10 || vm(i+1) == 25
            alfa_m = 1;
            alfa_n = 0.1;
        else
            [alfa_n,alfa_m,alfa_h,beta_n,beta_m,beta_h] = Parte1HH(vm(i+1));
        end
        
        n(i+1) = n(i) + (alfa_n*(1-n(i))-beta_n*(n(i)))*delta_t;
        m(i+1) = m(i) + (alfa_m*(1-m(i))-beta_m*(m(i)))*delta_t;
        h(i+1) = h(i) + (alfa_h*(1-h(i))-beta_h*(h(i)))*delta_t;
        
        time_vector(i+1) = time_vector(i) + delta_t;

    end

figure(1)
subplot(324)
plot(time_vector, I_Na/1000, time_vector, I_K/1000), 
legend('I_{Na^+}', 'I_{K^+}'), xlabel('t (ms)'), ylabel('I (mA)')
title('Correntes Iónicas')
grid on

subplot(326)
plot(time_vector, g_Na, time_vector, g_K) 
legend('g_{Na^+}', 'g_{K^+}'), xlabel('t (ms)'), ylabel('g (mS)')
title('Condutâncias Iónicas')
grid on

subplot(221)
plot(time_vector,Vm), xlabel('t (ms)'), ylabel('Vm (mV)')
title('Potencial de Membrana')
grid on

subplot(223)
plot(time_vector,n, time_vector, m, time_vector, h),legend('n', 'm', 'h'),xlabel('t (ms)')
title('Gating Probabilities')
grid on

subplot(322)
plot(time_vector,Im(1,1:Npoints)/1000),xlabel('t (ms)'), ylabel('I (mA)')
axis([0 t_total 0 Is/1000+0.001])
title('Corrente de Estímulo')
grid on

figure(2)
subplot(121), plot(time_vector,I_ionica/1000),xlabel('t (ms)'), ylabel('I (mA)')
title('Corrente Iónica')
subplot(222),plot(time_vector, I_L/1000),xlabel('t (ms)'), ylabel('I (mA)')
title('Corrente de leakage')
subplot(224), plot(time_vector,Ic/1000),xlabel('t (ms)'), ylabel('I (mA)') 
title('Corrente Capacitiva')
%% Alínea D
% Limiar de geração de potencias de ação

% Só devolve o limiar de geração de potencial se a Corrente de Estímulo for
% superior ao limiar
    
Int = (0:0.01:Is);
 
for j = 1:length(Int)
 
    Im = zeros(1,Npoints);
    Im(n_0:duration_Is/delta_t+n_0) = Int(j)*ones(1,duration_Is/delta_t+1);
        
    for i = 1:Npoints-1
        
        g_K(i) = gi_K*n(i)^4;
        g_Na(i) = gi_Na*m(i)^3*h(i);
        I_K(i) = g_K(i)*(Vm(i) -  E_K);
        I_Na(i) = g_Na(i)*(Vm(i) -  E_Na);
        I_L(i) = g_L*(Vm(i) - E_L);
        I_ionica(i) = I_Na(i) + I_K(i) + I_L(i); 
        Ic(i) = C_m*(Im(i) - I_ionica(i));
        vm(i+1) = vm(i) + (Ic(i))*delta_t/C_m;
        Vm(i+1) = vm(i+1) + V_r;
        
        if vm(i+1) == 10 || vm(i+1) == 25
            alfa_m = 1;
            alfa_n = 0.1;
        else
            [alfa_n,alfa_m,alfa_h,beta_n,beta_m,beta_h] = Parte1HH(vm(i+1));
        end
        
        n(i+1) = n(i) + (alfa_n*(1-n(i))-beta_n*(n(i)))*delta_t;
        m(i+1) = m(i) + (alfa_m*(1-m(i))-beta_m*(m(i)))*delta_t;
        h(i+1) = h(i) + (alfa_h*(1-h(i))-beta_h*(h(i)))*delta_t;
        
        time_vector(i+1) = time_vector(i) + delta_t;

    end
    
    maxV(j) = max(Vm);
    
end

k = find(maxV>0,1);

if k ~= 0
    
    figure(3)
    plot(Int,maxV), xlabel('I_m (mA)'), ylabel('max(V_m) (mV)')
    hold on
    plot(Int(k), maxV(k),'r*')
    title('Limiar de Geração de Potencias de Ação')
    disp(Int(k))
else
    
    disp(['A Corrente de Estímulo não é suficiente para induzir um' ...
        , ' Potencial de Ação.'])
    
end

