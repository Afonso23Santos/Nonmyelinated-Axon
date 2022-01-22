%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               BIOELETRICIDADE - PARTE 2                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Definir as constantes e par�metros
    V_r = -60;                                                             %(mV) - Potencial de repouso
    E_K = -72.10;                                                          %(mV) - Potencial de Nernst para o K
    E_Na = 52.4;                                                           %(mV) - Potencial de Nernst para o Na
    E_L = -49.187;                                                         %(mV) - Leakage equilibrium potential
    C_m = 1.0;                                                             %(uF/cm^2) - Capacidade da membrana
    gi_Na = 120.0;                                                         %(mS/cm^2) - Condut�ncia m�xima de Na
    gi_K = 36.0;                                                           %(mS/cm^2) - Condut�ncia m�xima de K
    g_L = 0.3;                                                             %(mS/cm^2) - Condut�ncia de leakage
    r_i = 0.1061e5;                                                        %(ohm/cm) - Resist�ncia interna da membrana
    r_e = 0.2357e4;                                                        %(ohm/cm) - Resist�ncia externa da membrana
    a = 0.03;                                                              %(cm) - Raio  
    c_m = 0.1885;                                                          %Capacitancia membranar

    
    comprimento = 30;                                                      %Comprimento da fibra em centimetros
    t_final = 30;                                                          %Dura��o total da simula��o em milisegundos
    t_pontos = t_final/delta_t;                                            %N�mero de pontos temporais
    x_pontos = comprimento/delta_x;                                        %N�mero de pontos espaciais
    
%% Criar as matrizes    
    t = zeros(1,t_pontos);                                                 %Vetor temporal
    x = zeros(1,x_pontos);                                                 %Vetor espacial
    I_m = zeros(t_pontos,x_pontos);                                        %Corrente transmembranar
    g_Na = zeros(t_pontos,x_pontos);                                       %Condut�ncia dos canais de Na
    g_K = zeros(t_pontos,x_pontos);                                        %Condut�ncia dos canais de K
    I_K = zeros(t_pontos,x_pontos);                                        %Corrente dos canais de K
    I_Na = zeros(t_pontos,x_pontos);                                       %Corrente dos canais de Na
    n = zeros(t_pontos,x_pontos);                                          %Gating probability n
    m = zeros(t_pontos,x_pontos);                                          %Gating probability m
    h = zeros(t_pontos,x_pontos);                                          %Gating probability h
    I_L = zeros(t_pontos,x_pontos);                                        %Corrente de leakage
    v_m = zeros(t_pontos,x_pontos);                                        %Potencial de menbrana reduzido
    V_m = v_m + V_r;                                                       %Potencial de membrana
    I_ion = zeros(t_pontos,x_pontos);                                      %Corrente I�nica
    I_c = zeros(t_pontos,x_pontos);                                        %Corrente Capacitiva
    I_i = zeros(t_pontos,x_pontos);                                        %Corrente Intracelular

    n(1,:) = 0.31768;                                                      %Condi��es iniciais para gating probability n
    m(1,:) = 0.05293;                                                      %Condi��es iniciais para gating probability m
    h(1,:) = 0.59612;                                                      %Condi��es iniciais para gating probability h  