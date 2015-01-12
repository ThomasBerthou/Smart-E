 #-*- coding: utf-8 -*-
#modele de bâtiment de type R6C2 developpé et validé dans la thèse de Thomas Berthou 2013
"""
Created on Wed Apr 09 15:01:10 2014
Modèle de bâtiement développé dans le thèse de Berthou 2013
@author: Thomas B
"""
def R6C2_c(data, x, debut, fin, init, delta, p_max_chauf, p_max_clim):
    import numpy as np
    f=data[debut:fin,0]
    fm=data[debut:fin,1]
    tc_chauf=data[debut:fin,2]
    te=data[debut:fin,3]
    occ=data[debut:fin,4]
    ve=data[debut:fin,5]
    tc_clim=data[debut:fin,6] #consigne de climatisation
  
    ph=np.zeros(fin-debut)
    tih=np.zeros(fin-debut)
    th=np.zeros(fin-debut)
    ts=np.zeros(fin-debut)
    tm=np.zeros(fin-debut)
    p1=np.zeros(fin-debut)
    p_min=0
    inter=3600*5 #pèriode entre deux interventions humaines sur les volets ou les fenètres
    #en seconde (3 heures)
    #initialisation
    v=x[10]
    a=x[9]
    #ph[0]=1000
    ph[0]=init[0]
    tih[0]=init[1]
    th[0]=init[2]
    ts[0]=init[3]
    tm[0]=init[4]
    k=0
    for i in range(1,fin-debut):
        
        k=k-1*delta
        #modèle simplifié de fermeture des volets
        if tih[i-1] > tc_clim[i] and occ[i] > 0 and k <=0 :
            f[i-1:i+inter/delta-1] = 0.1 * f[i-1:i+inter/delta-1]
            k = inter
            
        #modèle simplifié d'ouverture des fenètres    
        if tih[i-1] > tc_clim[i] and te[i-1] < tih[i-1] and occ[i]<0:
            ve[i] = 1
        
        th[i-1]=(tm[i-1]/x[4]+ te[i-1]/x[5]+fm[i-1])/(1/x[4]+1/x[5])
        ts[i-1]=(tih[i-1]/x[2]+tm[i-1]/x[3]+f[i-1]+a*x[8]*occ[i])/(1/x[2]+1/x[3])
    
        if tih[i-1] <= tc_chauf[i]:
                
            p1[i] = x[0]*(tc_chauf[i]-tih[i-1])/delta+(tih[i-1]-ts[i-1])/x[2]+(tih[i-1]-te[i-1])/x[6]-(1-a)*x[8]*occ[i]+v*ve[i]*(tih[i-1]-te[i-1])
    
            ph[i] = max(min(p1[i],p_max_chauf),p_min)
            
        elif tih[i-1] >= tc_clim[i]:
                
            p1[i] = x[0]*(tc_clim[i]-tih[i-1])/delta+(tih[i-1]-ts[i-1])/x[2]+(tih[i-1]-te[i-1])/x[6]-(1-a)*x[8]*occ[i]+v*ve[i]*(tih[i-1]-te[i-1])
    
            ph[i] = min(max(p1[i],-p_max_clim),p_min)     
                
        else:
            ph[i]=0
    
        tm[i] = ((th[i-1]-tm[i-1])/x[4]+(ts[i-1]-tm[i-1])/x[3])*delta/x[1]+tm[i-1]
    
        tih[i] = ((ts[i-1]-tih[i-1])/x[2]+(te[i-1]-tih[i-1])/x[6]+(ph[i]+(1-a)*x[8]*occ[i]-v*ve[i]*(tih[i-1]-te[i-1])))*delta/x[0]+tih[i-1]  
    
    

    return ph, tih , tm, f
