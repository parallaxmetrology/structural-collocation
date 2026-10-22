import numpy as np, pandas as pd, json, os
OUT = os.path.dirname(os.path.abspath(__file__)) + '/'
def nf(C, x, y, z):
    return 1-C[x, y]*C[x, z]/(C[x, x]*C[y, z])
def ev(C, x, y, z):
    return C[x, x]-C[x, y]*C[x, z]/C[y, z]
res = {}
w = pd.read_csv(OUT+'a2_windows.csv')
res['a2_n'] = len(w)
for v in ['wet', 'cy', 'cx', 'disp', 'isl', 'logmean']:
    C = np.cov(w[[o+'_'+v for o in ('R', 'IR', 'MW')]].values.T)
    r = {'nf_R': nf(C, 0, 1, 2), 'nf_IR': nf(C, 1, 0, 2), 'nf_MW': nf(C, 2, 0, 1)}
    r['diff'] = r['nf_R']-r['nf_MW']; r['ev_MW'] = ev(C, 2, 0, 1)
    res['a2_'+v] = r
p = pd.read_csv(OUT+'t7_points.csv')
med = p.rpeak.median(); res['median_peak'] = med
for st, sub in (('weak', p[p.rpeak <= med]), ('strong', p[p.rpeak > med])):
    for v in ['cy', 'cx']:
        X = sub[[o+'_'+v for o in ('R', 'IR', 'MW', 'G')]].values.T
        C = np.cov(X); c = np.corrcoef(X)
        R, I, M, G = 0, 1, 2, 3
        r = {'n': len(sub)}
        r['RGM_nf_R'] = nf(C, R, I, M); r['RGM_nf_IR'] = nf(C, I, R, M); r['RGM_nf_MW'] = nf(C, M, R, I)
        r['RGM_diff_R_MW'] = r['RGM_nf_R']-r['RGM_nf_MW']
        r['RMG_nf_R'] = nf(C, R, M, G); r['RMG_nf_MW'] = nf(C, M, R, G); r['RMG_nf_G'] = nf(C, G, R, M)
        r['RMG_diff_R_MW'] = r['RMG_nf_R']-r['RMG_nf_MW']
        r['RIRG_nf_R'] = nf(C, R, I, G); r['RIRG_nf_IR'] = nf(C, I, R, G); r['RIRG_nf_G'] = nf(C, G, R, I)
        r['IRMWG_nf_IR'] = nf(C, I, M, G); r['IRMWG_nf_MW'] = nf(C, M, I, G); r['IRMWG_nf_G'] = nf(C, G, I, M)
        # error cov of I and M with R and G as instruments
        cim = C[I, M]-C[I, R]*C[M, G]/C[R, G]
        r['cov_IR_MW'] = cim
        eI = ev(C, I, R, G); eM = ev(C, M, R, G)
        r['c_IR_MW_norm_RG'] = cim/np.sqrt(eI*eM)
        r['c_IR_MW_alt'] = (C[I, M]-C[I, G]*C[M, R]/C[G, R])  # symmetric partner
        eI2 = ev(C, I, R, M); eM2 = ev(C, M, R, I)
        r['c_IR_MW_norm_triplet_RGM'] = cim/np.sqrt(ev(C, I, R, M)*ev(C, M, R, I)) if False else None
        r['t1'] = c[R, G]*c[I, M]-c[R, I]*c[G, M]
        r['t2'] = c[R, G]*c[I, M]-c[R, M]*c[G, I]
        r['t3'] = c[R, I]*c[G, M]-c[R, M]*c[G, I]
        r['r'] = {a+'-'+b: c[i, j] for (a, i), (b, j) in [(('R', 0), ('IR', 1)), (('R', 0), ('MW', 2)), (('IR', 1), ('MW', 2)), (('R', 0), ('G', 3)), (('IR', 1), ('G', 3)), (('MW', 2), ('G', 3))]}
        res[st+'_'+v] = r
json.dump(res, open(OUT+'reimpl_results.json', 'w'), indent=1, default=float)
