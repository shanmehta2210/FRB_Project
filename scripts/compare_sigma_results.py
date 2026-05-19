import csv, math, numpy as np

def load_csv(path):
    with open(path) as f:
        return {row['FRB']: row for row in csv.DictReader(f)}

old = load_csv('galfit_metrics_summary.csv')
sig = load_csv('galfit_sigma_metrics_summary.csv')

def safe(v):
    try:
        v = v.replace('*','')
        return float(v)
    except:
        return None

def inc(q):
    if q is None: return None
    q0 = 0.2
    if q <= q0: return 90.0
    arg = (q**2 - q0**2)/(1 - q0**2)
    if arg < 0 or arg > 1: return None
    return math.degrees(math.acos(math.sqrt(arg)))

frbs = sorted(old.keys())

hdr = "{:<12} {:>10} {:>14} {:>6} {:>6} {:>6} {:>7} {:>7} {:>7} {:>7} {:>7} {:>7} {:>8} {:>8} {:>8}"

# ====== NO-PSF ======
print('='*130)
print('NO-PSF COMPARISON: without sigma (old) vs with sigma (new)')
print('='*130)
print(hdr.format('FRB','chi2_old','chi2_new','n_old','n_new','dn','Re_old','Re_new','dRe','ba_old','ba_new','dba','inc_old','inc_new','dinc'))
print('-'*130)

dba_np, dinc_np, dn_np, dre_np = [], [], [], []

for frb in frbs:
    o, s = old[frb], sig[frb]
    chi_o = safe(o['chi2nu_nopsf']); chi_s = safe(s['chi2nu_nopsf'])
    n_o = safe(o['n_nopsf']); n_s = safe(s['n_nopsf'])
    re_o = safe(o['re_nopsf']); re_s = safe(s['re_nopsf'])
    ba_o = safe(o['b_a_nopsf']); ba_s = safe(s['b_a_nopsf'])
    inc_o = inc(ba_o); inc_s = inc(ba_s)
    dn = (n_s - n_o) if n_o and n_s else None
    dre = (re_s - re_o) if re_o and re_s else None
    dba = (ba_s - ba_o) if ba_o and ba_s else None
    di = (inc_s - inc_o) if inc_o is not None and inc_s is not None else None
    if dba is not None: dba_np.append(dba)
    if di is not None: dinc_np.append(di)
    if dn is not None: dn_np.append(dn)
    if dre is not None: dre_np.append(dre)
    co = f'{chi_o:.3f}' if chi_o and chi_o < 100 else (f'{chi_o:.0f}' if chi_o else '?')
    cs = f'{chi_s:.3f}' if chi_s and chi_s < 100 else (f'{chi_s:.0f}' if chi_s else '?')
    print(f'{frb:<12} {co:>10} {cs:>14} {n_o:>6.2f} {n_s:>6.2f} {dn:>+6.2f} {re_o:>7.2f} {re_s:>7.2f} {dre:>+7.2f} {ba_o:>7.2f} {ba_s:>7.2f} {dba:>+7.3f} {inc_o:>8.1f} {inc_s:>8.1f} {di:>+8.1f}')

print()
print(f'NO-PSF STATISTICS (n={len(dba_np)}):')
print(f'  delta b/a:  mean={np.mean(dba_np):+.4f}  std={np.std(dba_np):.4f}  median={np.median(dba_np):+.4f}  max|d|={max(abs(x) for x in dba_np):.4f}')
print(f'  delta inc:  mean={np.mean(dinc_np):+.2f}  std={np.std(dinc_np):.2f}  median={np.median(dinc_np):+.2f}  max|d|={max(abs(x) for x in dinc_np):.2f} deg')
print(f'  delta n:    mean={np.mean(dn_np):+.3f}  std={np.std(dn_np):.3f}  median={np.median(dn_np):+.3f}')
print(f'  delta Re:   mean={np.mean(dre_np):+.3f}  std={np.std(dre_np):.3f}  median={np.median(dre_np):+.3f}')

# ====== WITH-PSF ======
print()
print('='*130)
print('WITH-PSF COMPARISON: without sigma (old) vs with sigma (new)')
print('='*130)
print(hdr.format('FRB','chi2_old','chi2_new','n_old','n_new','dn','Re_old','Re_new','dRe','ba_old','ba_new','dba','inc_old','inc_new','dinc'))
print('-'*130)

dba_p, dinc_p, dn_p, dre_p = [], [], [], []

for frb in frbs:
    o, s = old[frb], sig[frb]
    chi_o = safe(o['chi2nu_psf']); chi_s = safe(s['chi2nu_psf'])
    n_o = safe(o['n_psf']); n_s = safe(s['n_psf'])
    re_o = safe(o['re_psf']); re_s = safe(s['re_psf'])
    ba_o = safe(o['b_a_psf']); ba_s = safe(s['b_a_psf'])
    inc_o = inc(ba_o); inc_s = inc(ba_s)
    dn = (n_s - n_o) if n_o and n_s else None
    dre = (re_s - re_o) if re_o and re_s else None
    dba = (ba_s - ba_o) if ba_o and ba_s else None
    di = (inc_s - inc_o) if inc_o is not None and inc_s is not None else None
    if dba is not None: dba_p.append(dba)
    if di is not None: dinc_p.append(di)
    if dn is not None: dn_p.append(dn)
    if dre is not None: dre_p.append(dre)
    co = f'{chi_o:.3f}' if chi_o and chi_o < 100 else (f'{chi_o:.0f}' if chi_o else '?')
    cs = f'{chi_s:.3f}' if chi_s and chi_s < 100 else (f'{chi_s:.0f}' if chi_s else '?')
    ba_o_s = f'{ba_o:>7.2f}' if ba_o else '      ?'
    ba_s_s = f'{ba_s:>7.2f}' if ba_s else '      ?'
    inc_o_s = f'{inc_o:>8.1f}' if inc_o is not None else '       ?'
    inc_s_s = f'{inc_s:>8.1f}' if inc_s is not None else '       ?'
    di_s = f'{di:>+8.1f}' if di is not None else '       ?'
    dba_s = f'{dba:>+7.3f}' if dba is not None else '      ?'
    dn_s = f'{dn:>+6.2f}' if dn is not None else '     ?'
    dre_s = f'{dre:>+7.2f}' if dre is not None else '      ?'
    print(f'{frb:<12} {co:>10} {cs:>14} {n_o:>6.2f} {n_s:>6.2f} {dn_s} {re_o:>7.2f} {re_s:>7.2f} {dre_s} {ba_o_s} {ba_s_s} {dba_s} {inc_o_s} {inc_s_s} {di_s}')

print()
print(f'WITH-PSF STATISTICS (n={len(dba_p)}):')
print(f'  delta b/a:  mean={np.mean(dba_p):+.4f}  std={np.std(dba_p):.4f}  median={np.median(dba_p):+.4f}  max|d|={max(abs(x) for x in dba_p):.4f}')
print(f'  delta inc:  mean={np.mean(dinc_p):+.2f}  std={np.std(dinc_p):.2f}  median={np.median(dinc_p):+.2f}  max|d|={max(abs(x) for x in dinc_p):.2f} deg')
print(f'  delta n:    mean={np.mean(dn_p):+.3f}  std={np.std(dn_p):.3f}  median={np.median(dn_p):+.3f}')
print(f'  delta Re:   mean={np.mean(dre_p):+.3f}  std={np.std(dre_p):.3f}  median={np.median(dre_p):+.3f}')

# ====== CHI2 DIAGNOSTIC ======
print()
print('='*130)
print('CHI2/NU DIAGNOSTIC: Sigma-run galaxies with chi2/nu > 100 (likely invvar scale issue)')
print('='*130)
for frb in frbs:
    s = sig[frb]
    c1 = safe(s['chi2nu_nopsf']); c2 = safe(s['chi2nu_psf'])
    if (c1 and c1 > 100) or (c2 and c2 > 100):
        o = old[frb]
        co1 = safe(o['chi2nu_nopsf']); co2 = safe(o['chi2nu_psf'])
        co1s = f'{co1:.3f}' if co1 else '?'
        co2s = f'{co2:.3f}' if co2 else '?'
        print(f'  {frb}: sigma chi2 = {c1:.1f} / {c2:.1f}  |  old chi2 = {co1s} / {co2s}')

# Count how many changed b/a by > 0.05
big_ba_nopsf = sum(1 for x in dba_np if abs(x) > 0.05)
big_ba_psf = sum(1 for x in dba_p if abs(x) > 0.05)
big_inc_nopsf = sum(1 for x in dinc_np if abs(x) > 5)
big_inc_psf = sum(1 for x in dinc_p if abs(x) > 5)

print()
print('='*130)
print('IMPACT SUMMARY')
print('='*130)
print(f'  No-PSF: {big_ba_nopsf}/23 galaxies shifted b/a by >0.05, {big_inc_nopsf}/23 shifted inclination by >5 deg')
print(f'  W/ PSF: {big_ba_psf}/23 galaxies shifted b/a by >0.05, {big_inc_psf}/23 shifted inclination by >5 deg')
