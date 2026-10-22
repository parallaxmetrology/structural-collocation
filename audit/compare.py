import json, os
d = '/Users/russellparrish/Documents/russell/precip-disagreement/'
m = json.load(open(d+'audit/reimpl_results.json')); a = json.load(open(d+'data/probe_a2_results.json'))
t = json.load(open(d+'data/t7_results.json')); r = json.load(open(d+'data/t7r_results.json'))
rows = []
for v in ['wet', 'cy', 'cx', 'disp', 'isl', 'logmean']:
    for k in ['nf_R', 'nf_IR', 'nf_MW']:
        rows.append(('A2 '+v, k, m['a2_'+v][k], a[v]['point'][k]))
    rows.append(('A2 '+v, 'nf_R-nf_MW', m['a2_'+v]['diff'], a[v]['point']['nf_R']-a[v]['point']['nf_MW']))
rows.append(('median_peak', '', m['median_peak'], t['median_peak']))
for st in ['weak', 'strong']:
    for v in ['cy', 'cx']:
        x = m[st+'_'+v]; o = t[st][v]['point']; q = r[st][v]['point']; n = f'T7 {st} {v}'
        for k in ['RGM_nf_R', 'RGM_nf_IR', 'RGM_nf_MW', 'RGM_diff_R_MW', 'RMG_nf_R', 'RMG_nf_MW', 'RMG_nf_G', 'RMG_diff_R_MW']:
            rows.append((n, k, x[k], o[k]))
        rows.append((n, 'c_IR_MW', x['c_IR_MW_norm_RG'], o['c_IR_MW']))
        for k, k2 in [('RIRG_nf_R', 'nf_R|RIRG'), ('RIRG_nf_IR', 'nf_IR|RIRG'), ('RIRG_nf_G', 'nf_G|RIRG'), ('IRMWG_nf_IR', 'nf_IR|IRMWG'), ('IRMWG_nf_MW', 'nf_MW|IRMWG'), ('IRMWG_nf_G', 'nf_G|IRMWG'), ('t1', 't1'), ('t2', 't2'), ('t3', 't3')]:
            rows.append((n.replace('T7', 'T7R'), k, x[k], q[k2]))
        rows.append((n.replace('T7', 'T7R'), 'R2 nfR(RIRG)-nfR(RMWG)', x['RIRG_nf_R']-x['RMG_nf_R'], q['R2_nfR_RIRG_minus_RMWG']))
out = ['| block | quantity | re-impl | original | diff |', '|---|---|---|---|---|']
mx = 0
for b, k, x, o in rows:
    dd = x-o; mx = max(mx, abs(dd))
    out.append(f'| {b} | {k} | {x:.4f} | {o:.4f} | {dd:+.5f} |')
open(d+'audit/match_table.md', 'w').write('\n'.join(out)+f'\n\nmax abs diff {mx:.2e} over {len(rows)} quantities\n')
print(mx, len(rows), sorted(rows, key=lambda r: -abs(r[2]-r[3]))[:5])
