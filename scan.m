k_init=1e-18;

while (k_init<1e-16)
    k_fin=k_init;
    while (k_fin<6e-16)
        if(exist(strcat('.\dir_', num2str(k_init),'_',num2str(k_fin))))
            k_fin=1.5*k_fin;
        else
        mkdir(strcat('.\dir_', num2str(k_init),'_',num2str(k_fin)));
        P=setdefaultparams_steak(k_init,k_fin);
        steak_sim_dirichlet(P);
        k_fin=1.5*k_fin;
        end
    end
    k_init=1.5*k_init;
end