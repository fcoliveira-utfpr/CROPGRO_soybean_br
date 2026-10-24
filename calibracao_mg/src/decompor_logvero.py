"""Decompõe o log-post do melhor estado de cada cadeia por fonte (experimentos com o intercepto aleatório marginalizado,
catálogos, agregados e prior da literatura). Uso: DELTA8_FIXO=11 PRIOR_LIT=1 .venv/bin/python src/decompor_logvero.py a.h5 b.h5"""
import sys; sys.path.insert(0, __import__('os').path.dirname(__file__))
import calibrar_fenologia as cf, emcee, numpy as np, pandas as pd
c=cf.Calibracao(); o=c.obs; ix={n:i for i,n in enumerate(cf.NOMES)}
def comp(th):
    pr=c.prever(th); out={}
    s={k:np.exp(th[ix[k]]) for k in ('ln_s_r1','ln_s_r5','ln_s_r7','ln_s_cat','ln_s_agr')}
    sm={'r1':s['ln_s_r1'],'r5':s['ln_s_r5'],'r7':s['ln_s_r7'],'r8':s['ln_s_r7']}
    t2=np.exp(2*th[ix['ln_tau']])
    R=c.y[c.i_exp]-np.array([pr[e][i] for i,e in zip(c.i_exp,c.est_exp)])
    S=np.array([sm[e] for e in c.est_exp]); W=1/S**2
    fonte=o.fonte.values[c.i_exp]
    for f in np.unique(fonte):
        k=fonte==f; G=c.G[k]; G=G[:,G.sum(0)>0]
        A,B=W[k]@G,(W[k]*R[k])@G
        out['exp '+f[:22]]=-0.5*((W[k]*R[k]**2).sum()-(t2*B**2/(1+t2*A)).sum())-0.5*(np.log(S[k]**2).sum()+np.log1p(t2*A).sum())
    for agr,nome in [(False,'catalogos'),(True,'agregados')]:
        ll=0
        for e,g,y,a in zip(c.cat_est,c.cat_grupos,c.cat_y,c.cat_agr):
            if a!=agr: continue
            sg=s['ln_s_agr'] if a else s['ln_s_cat']; r=y-pr[e][g].mean()
            ll+=-0.5*(r/sg)**2-np.log(sg)
        out[nome]=ll
    pri=0
    for k,(mu,dp) in cf.PRIOR_LIT.items():
        for gm in (5,9):
            v=np.clip(th[ix['a_'+k]]+th[ix['b_'+k]]*(gm-7),*cf.FAIXAS[k]); pri+=-0.5*((v-mu)/dp)**2
    out['prior_lit']=pri
    return pd.Series(out)
res={}
for s in sys.argv[1:]:
    b=emcee.backends.HDFBackend(s,read_only=True); s=s.split('/')[-1]
    lp=b.get_log_prob(flat=True); th=b.get_chain(flat=True)[np.argmax(lp)]
    res[s]=comp(th); print(s,lp.max().round(1),res[s].sum().round(1))
d=pd.DataFrame(res); d['dif']=d.iloc[:,-1]-d.iloc[:,0]; print(d.round(1).sort_values('dif').to_string())
