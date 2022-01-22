function[alfa_n,alfa_m,alfa_h,beta_n,beta_m,beta_h] = Parte1HH(vm)
   
   %Alinea a
    alfa_n = 0.01*(10-vm)./(exp((10-vm)./10)-1);
    alfa_m = 0.1*(25-vm)./(exp((25-vm)./10)-1);
    alfa_h = 0.07*exp(-vm/20);
    beta_n = 0.125*exp(-vm/80);
    beta_m = 4*exp(-vm/18);
    beta_h = 1./(exp((30-vm)/10)+1);
    
end

%%