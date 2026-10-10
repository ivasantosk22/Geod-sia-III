# -*- coding: utf-8 -*-
"""Controle de qualidade de observações GNSS no estilo do TEQC (+qc), em Python puro (GPS L1/L2).
Reimplementação a partir das definições publicadas (Estey & Meertens, 1999); não é o executável original."""
import numpy as np
from matplotlib.figure import Figure
import matplotlib.dates as mdates
import observaveis as O
from nucleo import _png, GPS0, ler_nav

C, F1, F2 = 299792458.0, 1575.42e6, 1227.60e6
ALFA = (F1 / F2) ** 2                       # (77/60)^2
MU, WE = 3.986005e14, 7.2921151467e-5


# ------------------------------------------------------------------ leitura (obs + navegação)
def carregar_todos(dados, nome, nav_extra=None):
    """Lê todos os dias/estações do arquivo enviado; usa a navegação do próprio zip ou a enviada à parte."""
    fs = O.carregar_todos(dados, nome, com_nav=True)
    for f in fs:
        txt = nav_extra or f.pop('nav_txt', None); f.pop('nav_txt', None)
        f['efs'] = [e for e in ler_nav(txt) if e.sat[0] == 'G'] if txt else []
    return fs


def carregar(dados, nome, nav_extra=None):
    return carregar_todos(dados, nome, nav_extra)[0]


# ------------------------------------------------------------------ geometria (elevação / azimute)
def _geod(x, y, z):
    a, fl = 6378137.0, 1 / 298.257223563; e2 = fl * (2 - fl)
    lon = np.arctan2(y, x); p = np.hypot(x, y); lat = np.arctan2(z, p * (1 - e2))
    for _ in range(6):
        N = a / np.sqrt(1 - e2 * np.sin(lat) ** 2); h = p / np.cos(lat) - N; lat = np.arctan2(z, p * (1 - e2 * N / (N + h)))
    return lat, lon


def _azel(sat, rx):
    lat, lon = _geod(*rx); d = sat - np.array(rx)
    e = -np.sin(lon) * d[..., 0] + np.cos(lon) * d[..., 1]
    n = -np.sin(lat) * np.cos(lon) * d[..., 0] - np.sin(lat) * np.sin(lon) * d[..., 1] + np.cos(lat) * d[..., 2]
    u = np.cos(lat) * np.cos(lon) * d[..., 0] + np.cos(lat) * np.sin(lon) * d[..., 1] + np.sin(lat) * d[..., 2]
    return np.degrees(np.arctan2(e, n)) % 360, np.degrees(np.arctan2(u, np.hypot(e, n)))


def _pos(e, sow):
    m = e.m; Crs, dn, M0 = m[0, 1:4]; Cuc, ec, Cus, ra = m[1]; Toe, Cic, O0, Cis = m[2]; i0, Crc, w, Od = m[3]; idot = m[4, 0]
    A = ra ** 2; tk = (sow - Toe + 302400) % 604800 - 302400
    M = M0 + (np.sqrt(MU / A ** 3) + dn) * tk; E = M.copy()
    for _ in range(10):
        E = M + ec * np.sin(E)
    ph = np.arctan2(np.sqrt(1 - ec ** 2) * np.sin(E), np.cos(E) - ec) + w; s2, c2 = np.sin(2 * ph), np.cos(2 * ph)
    u = ph + Cuc * c2 + Cus * s2; r = A * (1 - ec * np.cos(E)) + Crc * c2 + Crs * s2; i = i0 + idot * tk + Cic * c2 + Cis * s2
    Om = O0 + (Od - WE) * tk - WE * Toe; x, y = r * np.cos(u), r * np.sin(u)
    return np.stack([x * np.cos(Om) - y * np.sin(Om) * np.cos(i), x * np.sin(Om) + y * np.cos(Om) * np.cos(i), y * np.sin(i)], -1)


def elevacao(f):
    nep, ns = f['A'].shape[:2]; az = np.full((nep, ns), np.nan); el = az.copy()
    if not f['efs'] or not f['xyz']:
        return az, el
    sow = ((f['ep'][0] - GPS0).total_seconds() + f['t']) % 604800
    for j, s in enumerate(f['sats']):
        rs = [e for e in f['efs'] if e.sat == s]
        if not rs:
            continue
        toe = np.array([e.m[2, 0] for e in rs]); d = (sow[:, None] - toe[None, :] + 302400) % 604800 - 302400
        k = np.abs(d).argmin(1); ok = np.abs(d[np.arange(nep), k]) <= 4 * 3600; pos = np.full((nep, 3), np.nan)
        for q in np.unique(k[ok]):
            m = ok & (k == q); pos[m] = _pos(rs[q], sow[m])
        az[:, j], el[:, j] = _azel(pos, f['xyz'])
    return az, el


# ------------------------------------------------------------------ controle de qualidade
def analisar(f, mask=10.0, iod_lim=400.0, mp_lim=5.0, min_arco=20):
    A = f['A']; t = f['t']; nep, ns = A.shape[:2]
    l1, l2, p1, p2, s1, s2 = A[:, :, 0] * C / F1, A[:, :, 1] * C / F2, A[:, :, 2], A[:, :, 3], A[:, :, 4], A[:, :, 5]
    okph = np.isfinite(l1) & np.isfinite(l2)                                   # fase em L1 e L2
    ok1, ok2 = okph & np.isfinite(p1), okph & np.isfinite(p2)                  # + código C1 (MP1) / + código C2 (MP2)
    passo = float(np.median(np.diff(t))) if nep > 1 else 30.0
    az, el = elevacao(f); tem_el = bool(np.isfinite(el).any())
    esp = (el >= mask) if tem_el else ok1.copy()
    ac = (np.nan_to_num(el, nan=-90) >= mask) if tem_el else np.ones_like(ok1)
    mp1 = p1 - (1 + 2 / (ALFA - 1)) * l1 + (2 / (ALFA - 1)) * l2                               # MP1 (m), com viés
    mp2 = p2 - (2 * ALFA / (ALFA - 1)) * l1 + (2 * ALFA / (ALFA - 1) - 1) * l2                 # MP2 (m), com viés
    ion = (l1 - l2) / (ALFA - 1)                                                                 # atraso ionosférico em L1 (m), a menos de constante
    dt = np.diff(t)[:, None] if nep > 1 else np.ones((0, 1))
    cons = np.zeros_like(okph); cons[1:] = okph[1:] & okph[:-1] & (dt[:, 0] <= 2.5 * passo)[:, None]
    iod = np.full_like(ion, np.nan); iod[1:] = (ion[1:] - ion[:-1]) / dt; iod[~cons] = np.nan   # m/s
    dm = np.full_like(mp1, np.nan); dm[1:] = np.maximum(np.abs(mp1[1:] - mp1[:-1]), np.abs(mp2[1:] - mp2[:-1]))
    sl_iod = cons & (np.abs(iod) * 6000 > iod_lim); sl_mp = cons & (dm > mp_lim); slip = sl_iod | sl_mp
    # arcos contínuos: novo arco no início, após lacuna ou salto; viés removido pela média do arco
    arc = np.cumsum(okph & (~cons | slip), axis=0) * okph
    mp1c, mp2c = np.full_like(mp1, np.nan), np.full_like(mp2, np.nan)
    for j in range(ns):
        for mp, mc, okk in ((mp1, mp1c, ok1), (mp2, mp2c, ok2)):
            v = okk[:, j]
            if not v.any():
                continue
            ids = arc[v, j]; n = np.bincount(ids); med = np.bincount(ids, weights=mp[v, j]) / np.maximum(n, 1)
            mc[v, j] = np.where(n[ids] >= min_arco, mp[v, j] - med[ids], np.nan)
    sel1, sel2 = ok1 & ac, ok2 & ac
    rms = lambda x, m: float(np.sqrt(np.nanmean(x[m & np.isfinite(x)] ** 2))) if (m & np.isfinite(x)).any() else np.nan
    mean = lambda x, m: float(np.nanmean(x[m & np.isfinite(x)])) if (m & np.isfinite(x)).any() else np.nan
    n_esp, n_obs, n_sl = int(esp.sum()), int((sel1 & esp).sum()), int((slip & okph & ac).sum())
    res = dict(estacao=f['estacao'], data=f['ep'][0], fim=f['ep'][-1], intervalo=passo, nep=nep, nsat=int(okph.any(0).sum()), tem_el=tem_el,
               esperado=n_esp, obs=n_obs, compl=100 * n_obs / n_esp if n_esp else np.nan, mp1=rms(mp1c, sel1), mp2=rms(mp2c, sel2),
               sn1=mean(s1, okph & ac), sn2=mean(s2, okph & ac), sl_iod=int((sl_iod & ac).sum()), sl_mp=int((sl_mp & ac).sum()), slips=n_sl,
               o_slps=n_obs / n_sl if n_sl else np.inf, sl_1000=1000 * n_sl / max(n_obs, 1), el_media=mean(el, okph),
               rec=f.get('rec', ''), ant=f.get('ant', ''))
    sats = []
    for j, s in enumerate(f['sats']):
        e_j = int(esp[:, j].sum()); o_j = int((sel1[:, j] & esp[:, j]).sum())
        sats.append([s, mean(el[:, j], okph[:, j]), e_j, o_j, 100 * o_j / e_j if e_j else np.nan, rms(mp1c[:, j], sel1[:, j]), rms(mp2c[:, j], sel2[:, j]),
                     mean(s1[:, j], okph[:, j] & ac[:, j]), mean(s2[:, j], okph[:, j] & ac[:, j]), int((slip[:, j] & okph[:, j] & ac[:, j]).sum())])
    return dict(res=res, sats=sats, az=az, el=el, mp1c=mp1c, mp2c=mp2c, iod=iod * 6000, slip=slip, sl_iod=sl_iod, ok=ok1, esp=esp,
                s1=s1, s2=s2, mask=mask, iod_lim=iod_lim, tem_el=tem_el)


CAB_RES = ['Estação', 'Data (início)', 'Épocas', 'Intervalo (s)', 'Obs. esperadas', 'Obs. completas', 'Completude (%)', 'MP1 rms (m)', 'MP2 rms (m)',
           'SN1 médio', 'SN2 médio', 'Saltos (IOD+MP)', 'Obs. por salto', 'Saltos / 1000 obs']
CAB_SAT = ['PRN', 'Elev. média (°)', 'Esperadas', 'Completas', 'Completude (%)', 'MP1 rms (m)', 'MP2 rms (m)', 'SN1', 'SN2', 'Saltos']


def _f(x, d=2):
    return '—' if x is None or (isinstance(x, float) and not np.isfinite(x)) else f'{x:.{d}f}'


def linha_res(R):
    r = R['res']
    return [r['estacao'], f"{r['data']:%d/%m/%Y}", r['nep'], f"{r['intervalo']:g}", r['esperado'] if r['tem_el'] else '—', r['obs'], _f(r['compl'], 1),
            _f(r['mp1']), _f(r['mp2']), _f(r['sn1'], 1), _f(r['sn2'], 1), r['slips'], '∞' if not np.isfinite(r['o_slps']) else f"{r['o_slps']:.0f}", _f(r['sl_1000'], 2)]


def linhas_sat(R):
    return [[s[0], _f(s[1], 1), s[2], s[3], _f(s[4], 1), _f(s[5]), _f(s[6]), _f(s[7], 1), _f(s[8], 1), s[9]] for s in R['sats']]


# ------------------------------------------------------------------ gráficos
def graficos(f, R, sat, mask):
    out = []; ep = f['ep']; j = f['sats'].index(sat); tit = f"{f['estacao']} · {ep[0]:%d/%m/%Y}"
    el, az, m1, m2 = R['el'], R['az'], R['mp1c'], R['mp2c']
    if R['tem_el']:
        g = Figure(figsize=(11, 5.2)); sc = None
        for k, (nome, mp) in enumerate((('MP1', m1), ('MP2', m2))):
            ax = g.add_subplot(1, 2, k + 1, projection='polar'); ax.set_theta_zero_location('N'); ax.set_theta_direction(-1)
            m = np.isfinite(mp) & np.isfinite(el); idx = np.where(m)
            sc = ax.scatter(np.radians(az[m][::3]), 90 - el[m][::3], c=np.abs(mp[m][::3]), s=2, cmap='plasma', vmin=0, vmax=1)
            ax.set_ylim(0, 90); ax.set_yticks([0, 30, 60, 90]); ax.set_yticklabels(['90°', '60°', '30°', '0°']); ax.set_title(f'{nome} no céu — {tit}', fontsize=10)
        g.colorbar(sc, ax=g.axes, shrink=.8, label='|MP| (m)'); out.append(('Multipath no céu (azimute × elevação)', _png(g)))

        faixas = np.arange(0, 91, 5); cen = faixas[:-1] + 2.5
        def por_faixa(x, fun):
            return [fun(x[(el >= a) & (el < a + 5) & np.isfinite(x)]) if ((el >= a) & (el < a + 5) & np.isfinite(x)).sum() > 5 else np.nan for a in faixas[:-1]]
        g = Figure(figsize=(10, 7)); a1, a2 = g.subplots(2, 1, sharex=True)
        a1.plot(cen, por_faixa(m1, lambda v: np.sqrt(np.mean(v ** 2))), '-o', label='MP1'); a1.plot(cen, por_faixa(m2, lambda v: np.sqrt(np.mean(v ** 2))), '-o', label='MP2')
        a1.set(ylabel='RMS do multipath (m)', title=f'Multipath e relação sinal-ruído × elevação — {tit}'); a1.axvline(mask, color='gray', ls='--'); a1.legend(); a1.grid(alpha=.3)
        a2.plot(cen, por_faixa(R['s1'], np.mean), '-o', label='SN1'); a2.plot(cen, por_faixa(R['s2'], np.mean), '-o', label='SN2')
        a2.set(ylabel='SNR médio (unid. do receptor)', xlabel='Elevação (°)'); a2.axvline(mask, color='gray', ls='--'); a2.legend(); a2.grid(alpha=.3)
        g.tight_layout(); out.append(('Multipath e SNR por elevação', _png(g)))

        g = Figure(figsize=(11, 5)); a1, a2 = g.subplots(1, 2, gridspec_kw={'width_ratios': [1, 2]})
        pe = [int(((el >= a) & (el < a + 5)).sum()) for a in faixas[:-1]]; po = [int(((el >= a) & (el < a + 5) & R['ok']).sum()) for a in faixas[:-1]]
        a1.plot([100 * o / e if e else np.nan for o, e in zip(po, pe)], cen, '-o'); a1.set(xlabel='Completude (%)', ylabel='Elevação (°)', xlim=(0, 105), title='Completude por elevação'); a1.grid(alpha=.3)
        a2.bar([s[0] for s in R['sats']], [0 if not np.isfinite(s[4]) else s[4] for s in R['sats']], color='#2e8b57'); a2.set(ylabel='Completude (%)', ylim=(0, 105), title='Completude por satélite (acima da máscara)')
        a2.tick_params(axis='x', labelrotation=90, labelsize=7); a2.grid(axis='y', alpha=.3); g.tight_layout(); out.append(('Completude dos dados', _png(g)))
    g = Figure(figsize=(10, 8)); a1, a2, a3 = g.subplots(3, 1, sharex=True); E = np.array(ep)
    for ax, mp, nome in ((a1, m1, 'MP1 (m)'), (a2, m2, 'MP2 (m)')):
        v = np.isfinite(mp[:, j]); ax.plot(E[v], mp[v, j], '.', ms=2, color='#1f77b4'); ax.set_ylabel(nome); ax.grid(alpha=.3)
    a1.set_title(f'Séries do {sat} — {tit}'); v = np.isfinite(R['iod'][:, j]); a3.plot(E[v], R['iod'][v, j], '.', ms=2, color='#2ca02c')
    s = R['slip'][:, j] & v; a3.plot(E[s], R['iod'][s, j], 'rx', ms=7, label='salto'); a3.set(ylabel='IOD (cm/min)', xlabel='Hora (GPS)'); a3.grid(alpha=.3); a3.legend()
    fmt = '%H:%M' if (ep[-1] - ep[0]).total_seconds() <= 86400 else '%d/%m\n%H:%M'; a3.xaxis.set_major_formatter(mdates.DateFormatter(fmt))
    g.tight_layout(); out.append((f'Séries de MP1, MP2 e IOD do {sat}', _png(g)))
    return out