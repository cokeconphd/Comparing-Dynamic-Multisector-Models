1) Main code: main_n26_logs.m 

Simulates pairwise correlations (CORRELATION=1)
Simulates spectral densities (SPECTRALST=1)
Simulates Impulse Responses (VAR=1)

2) Other codes:

i) n_1k_eq_ss.m: equations to solve the nonlinear system of equations when RHO<-1
ii) simu_1st_original.m: simulates series from the models
iii) solab_original: Paul Kleins code to solve the system of first difference equations
iv) VARmakexy: construct VAR matrices to just run OLS later

3) Matlab data

i)   alfa26: capital shares
ii)  av_sell_theta26: average capital flow matrix (THETA)
iii) gamma26: Input output table

4) Excel spreadsheets

i) m_pii_26_"model".xls: these contains the policy rules and law of motion of states
 from carvalho and pierre's models with random walk TPFs

ii) m_pii_26_"model"_099.xls: these contains the policy rules and law of motion of states
 from carvalho and pierre's models with stationary TFP processes (varho=0.99)

iii) ipl2: Series of IP index for different sectors 1972-2007

iv) 1997DetailedItemOutput: Contains the observed sectors shares



