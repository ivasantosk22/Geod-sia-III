# -*- coding: utf-8 -*-
"""Movimento do polo (Aula 4): leitura de Boletins B/C04 do IERS, coordenadas polares (θ, d), distância linear,
média móvel, indicador de descontinuidade e gráficos (xp e yp, trajetória, polar, 3D)."""
import io, re, urllib.request
from datetime import date, timedelta
import numpy as np
import matplotlib.dates as mdates
from matplotlib.figure import Figure
from nucleo import _png

R_MEDIO = 6371000.0                        # raio médio da Terra (m)
MAS = np.pi / (180 * 3600 * 1000)          # 1 mas em radianos
MJD0 = date(1858, 11, 17)
EVENTOS = [('Chile (Maule), Mw 8,8', '2010-02-27'), ('Japão (Tohoku), Mw 9,0', '2011-03-11'), ('Sumatra, Mw 8,6', '2012-04-11'),
           ('Chile (Iquique), Mw 8,2', '2014-04-01'), ('Chile (Illapel), Mw 8,3', '2015-09-16'), ('Alasca, Mw 8,2', '2021-07-29'),
           ('Turquia, Mw 7,8', '2023-02-06'), ('Kamchatka, Mw 8,8 (29/07 23h24 UTC; 1º valor de 0h UTC: 30/07)', '2025-07-30')]


def mjd_data(m):
    return MJD0 + timedelta(days=int(m))


def data_mjd(d):
    return (d - MJD0).days


def soma_meses(d, n):
    m = d.year * 12 + d.month - 1 + n; y, mm = divmod(m, 12); mm += 1
    ult = [31, 29 if (y % 4 == 0 and (y % 100 or y % 400 == 0)) else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31][mm - 1]
    return date(y, mm, min(d.day, ult))


# ------------------------------------------------------------------ leitura dos boletins
def _num(t):
    try:
        return float(t)
    except ValueError:
        return None


def ler_boletim(texto):
    """Retorna dict(linhas=[(mjd, x_mas, y_mas, ut1utc_ms, xerr, yerr, final)], taiutc=[(m0, m1, s)], formato)."""
    L = texto.replace('\r', '').split('\n'); rows = []; fmt = None; prelim = False
    cab = next((l for l in L if l.strip()), '')
    if cab.startswith('MJD;Year;Month;Day;Type;x_pole'):                 # Boletim B em CSV (IERS, ';' como separador, mas)
        col = {n: i for i, n in enumerate(cab.strip().split(';'))}
        for l in L[1:]:
            t = l.strip().split(';')
            if len(t) < 9 or not t[0].isdigit():
                continue
            g = lambda n: _num(t[col[n]]) if n in col and col[n] < len(t) and t[col[n]] != '' else None
            x, y = g('x_pole'), g('y_pole')
            if x is None or y is None:
                continue
            u = g('UT1-UTC')
            rows.append((int(t[0]), x, y, u if u is not None else np.nan, g('sigma_x_pole') or np.nan, g('sigma_y_pole') or np.nan, t[col['Type']] == 'final'))
        return dict(linhas=rows, taiutc=[], formato='Boletim B (CSV, mas)')
    ini = next((i for i, l in enumerate(L) if re.search(r'1\s*-\s*DAILY FINAL VALUES OF x, y', l)), None)
    if ini is not None:                                               # Boletim B (2009+): seção 1, em mas e ms
        fim = next((i for i in range(ini + 1, len(L)) if re.match(r'^\s*2\s*-\s', L[i])), len(L)); fmt = 'Boletim B (mas)'
        for l in L[ini:fim]:
            prelim = prelim or 'Preliminary' in l
            t = l.split()
            if len(t) >= 7 and re.fullmatch(r'\d{4}', t[0]) and re.fullmatch(r'\d{5}', t[3]):
                v = [_num(x) for x in t[4:]]
                if None not in v[:3]:
                    rows.append((int(t[3]), v[0], v[1], v[2], v[5] if len(v) > 5 else np.nan, v[6] if len(v) > 6 else np.nan, not prelim))
    else:
        i2 = next((i for i, l in enumerate(L) if re.search(r'2\s*-\s*SMOOTHED VALUES OF', l)), None)
        if i2 is not None:                                            # Boletim B antigo (<2009): "MÊS dia MJD x y UT1..." em " e s
            i3 = next((i for i in range(i2 + 1, len(L)) if re.match(r'^\s*3\s*-\s', L[i])), len(L)); fmt = 'Boletim B antigo (")'
            for l in L[i2:i3]:
                t = l.split()
                if len(t) >= 6 and re.fullmatch(r'[A-Za-z]{3}', t[0]) and re.fullmatch(r'\d{5}', t[2]):
                    v = [_num(x) for x in t[3:6]]
                    if None not in v:
                        rows.append((int(t[2]), v[0] * 1000, v[1] * 1000, v[2] * 1000, np.nan, np.nan, True))
        else:                                                         # série C04 / outras: "ano mês dia MJD x y UT1-UTC ..."
            fmt = 'série tabular (C04)'
            for l in L:
                t = l.split()
                if len(t) >= 6 and re.fullmatch(r'\d{4}', t[0]) and re.fullmatch(r'\d{5}', t[3]) and re.fullmatch(r'\d{1,2}', t[1]):
                    v = [_num(x) for x in t[4:]]
                    if None not in v[:2]:
                        rows.append((int(t[3]), v[0], v[1], v[2] if len(v) > 2 and v[2] is not None else np.nan, np.nan, np.nan, True))
            if rows and np.nanmedian(np.abs([r[1] for r in rows])) < 2:   # arcsec / s -> mas / ms
                rows = [(m, x * 1000, y * 1000, u * 1000, xe, ye, f) for m, x, y, u, xe, ye, f in rows]
    tai = [(int(a), int(b), int(s)) for s, a, b in re.findall(r'TAI\s*-\s*UTC\s*=\s*(\d+)\s*s\.?\s*during period \[(\d+)\s*-\s*(\d+)\]', texto)]
    return dict(linhas=rows, taiutc=tai, formato=fmt)


def juntar(boletins):
    """Une vários boletins; em dia repetido, valor final prevalece sobre preliminar."""
    d = {}
    for b in boletins:
        for r in b['linhas']:
            if r[0] not in d or (r[6] and not d[r[0]][6]):
                d[r[0]] = r
    ordem = sorted(d)
    A = np.array([d[m][1:6] for m in ordem], float) if ordem else np.zeros((0, 5))
    return dict(mjd=np.array(ordem, int), x=A[:, 0], y=A[:, 1], ut1=A[:, 2], xe=A[:, 3], ye=A[:, 4], final=np.array([d[m][6] for m in ordem], bool),
                taiutc=[t for b in boletins for t in b['taiutc']])


def ut1_utc_ms(dados, mjd):
    """UT1-UTC (ms) interpolado linearmente em um MJD fracionário (None se fora do intervalo)."""
    m = dados['mjd']
    if len(m) < 2 or mjd < m[0] or mjd > m[-1]:
        return None
    return float(np.interp(mjd, m, dados['ut1']))


def tai_utc_boletim(dados, mjd):
    for a, b, s in dados['taiutc']:
        if a <= mjd <= b:
            return s
    return None


# ------------------------------------------------------------------ boletins necessários / download
def numero_boletim(ano, mes):
    """Nº do Boletim B (formato 2009+) publicado no mês dado: 290 = abril/2012, um por mês."""
    return 290 + (ano - 2012) * 12 + (mes - 4)


def boletins_necessarios(ini, fim):
    """Valores finais de um mês K estão no boletim publicado no mês K+2."""
    a, b = soma_meses(date(ini.year, ini.month, 1), 2), soma_meses(date(fim.year, fim.month, 1), 2)
    return list(range(numero_boletim(a.year, a.month), numero_boletim(b.year, b.month) + 1))


def baixar_boletins(numeros, timeout=12):
    """Tenta baixar bulletinb-N.txt do IERS Data Center. Retorna ({n: texto}, [erros])."""
    ok, erros = {}, []
    for n in numeros:
        try:
            req = urllib.request.Request(f'https://datacenter.iers.org/data/207/bulletinb-{n}.txt', headers={'User-Agent': 'efemerides-app/1.0'})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                ok[n] = r.read().decode('latin-1')
        except Exception as e:
            erros.append(f'bulletinb-{n}: {e}')
            if len(erros) >= 2:
                break                                                 # sem rede: não insistir em todos
    return ok, erros


# ------------------------------------------------------------------ análise
def analisar(dados, evento, meses=6, janela_ma=7):
    ini, fim = soma_meses(evento, -meses), soma_meses(evento, meses)
    m0, m1 = data_mjd(ini), data_mjd(fim)
    sel = (dados['mjd'] >= m0) & (dados['mjd'] <= m1)
    if sel.sum() < 10:
        raise ValueError(f'poucos dados no período {ini:%d/%m/%Y} a {fim:%d/%m/%Y} ({int(sel.sum())} dias). '
                         f'Boletins necessários: {boletins_necessarios(ini, fim)[0]} a {boletins_necessarios(ini, fim)[-1]}')
    mjd = dados['mjd'][sel]; x, y = dados['x'][sel], dados['y'][sel]
    d = np.hypot(x, y); th = np.arctan2(y, x) % (2 * np.pi)
    d_m = d * MAS * R_MEDIO
    dx, dy = np.diff(x), np.diff(y); seq = np.diff(mjd) == 1
    var = np.full(len(x), np.nan); var[1:] = np.where(seq, np.hypot(dx, dy), np.nan)               # variação diária (mas/dia)
    var_cm = var * MAS * R_MEDIO * 100
    w = max(3, int(janela_ma) | 1)                                                               # janela ímpar
    def ma(v):
        o = np.full(len(v), np.nan); k = np.ones(w) / w
        if len(v) >= w:
            o[w // 2: len(v) - w // 2] = np.convolve(v, k, mode='valid')
        return o
    mx, my, md = ma(x), ma(y), ma(d)
    ev = int(np.argmin(np.abs(mjd - data_mjd(evento))))
    antes, depois = np.arange(len(x)) < ev, np.arange(len(x)) >= ev
    def est(v, m):
        t = (mjd[m] - mjd[m][0]).astype(float)
        return v[m].mean(), v[m].std(ddof=1), (np.polyfit(t, v[m], 1)[0] if m.sum() > 2 else np.nan)
    # indicador de descontinuidade: 2ª diferença (aceleração) no evento ± 2 dias contra o restante da janela
    a = np.full(len(x), np.nan); a[1:-1] = np.hypot(x[2:] - 2 * x[1:-1] + x[:-2], y[2:] - 2 * y[1:-1] + y[:-2])
    perto = np.abs(np.arange(len(x)) - ev) <= 2; longe = ~(np.abs(np.arange(len(x)) - ev) <= 4) & np.isfinite(a)
    base = np.nanmedian(a[longe]); mad = 1.4826 * np.nanmedian(np.abs(a[longe] - base)) or np.nan
    a_ev = np.nanmax(a[perto]) if np.isfinite(a[perto]).any() else np.nan
    v_ev = np.nanmax(var[perto]) if np.isfinite(var[perto]).any() else np.nan
    est_tab = []
    for nome, v in (('xp (mas)', x), ('yp (mas)', y), ('d (mas)', d)):
        e1, e2 = est(v, antes), est(v, depois)
        est_tab.append([nome, f'{e1[0]:.2f}', f'{e1[1]:.2f}', f'{e1[2]:+.3f}', f'{e2[0]:.2f}', f'{e2[1]:.2f}', f'{e2[2]:+.3f}'])
    ind = dict(a_evento=a_ev, a_mediana=base, a_sigma=mad, razao=(a_ev - base) / mad if np.isfinite(mad) else np.nan, var_evento=v_ev,
               var_mediana=np.nanmedian(var), xerr=np.nanmean(dados['xe'][sel]), yerr=np.nanmean(dados['ye'][sel]), dia_evento=mjd_data(mjd[ev]))
    return dict(mjd=mjd, x=x, y=y, d=d, th=th, d_m=d_m, var=var, var_cm=var_cm, mx=mx, my=my, md=md, ev=ev, evento=evento, ini=ini, fim=fim,
                est=est_tab, ind=ind, janela=w, ut1=dados['ut1'][sel], final=dados['final'][sel], n=len(x), lacunas=int((~seq).sum()))


CAB_EST = ['', 'Média antes', 'Desvio antes', 'Tendência antes (mas/dia)', 'Média depois', 'Desvio depois', 'Tendência depois (mas/dia)']


def csv_polo(R):
    out = ['data,MJD,xp_mas,yp_mas,d_mas,theta_rad,theta_graus,d_m,variacao_mas_dia,variacao_cm_dia,media_movel_xp_mas,media_movel_yp_mas,media_movel_d_mas']
    f = lambda v: '' if not np.isfinite(v) else f'{v:.6g}'
    for i in range(R['n']):
        out.append(f"{mjd_data(R['mjd'][i]):%Y-%m-%d},{R['mjd'][i]},{R['x'][i]:.3f},{R['y'][i]:.3f},{R['d'][i]:.3f},{R['th'][i]:.6f},"
                   f"{np.degrees(R['th'][i]):.3f},{R['d_m'][i]:.4f},{f(R['var'][i])},{f(R['var_cm'][i])},{f(R['mx'][i])},{f(R['my'][i])},{f(R['md'][i])}")
    return '\n'.join(out)


# ------------------------------------------------------------------ gráficos
def graficos(R):
    out = []; dt = [mjd_data(m) for m in R['mjd']]; ev = R['ev']; dev = dt[ev]
    tit = f"Evento em {R['evento']:%d/%m/%Y}  (IERS, mas)"
    def datas(ax):
        ax.axvline(dev, color='crimson', ls='--', lw=1.2, label='evento'); ax.grid(alpha=.3)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%Y'))
    g = Figure(figsize=(10, 6.5)); a1, a2 = g.subplots(2, 1, sharex=True)
    for ax, v, mv, nome in ((a1, R['x'], R['mx'], 'xp (mas)'), (a2, R['y'], R['my'], 'yp (mas)')):
        ax.plot(dt, v, '.', ms=3, color='#1f77b4', label='diário'); ax.plot(dt, mv, '-', color='#ff7f0e', lw=1.8, label=f"média móvel ({R['janela']} d)")
        ax.set_ylabel(nome); datas(ax)
    a1.set_title('Componentes do polo: xp e yp — ' + tit); a2.set_xlabel('Data (UTC)'); a1.legend(ncol=3, fontsize=8); g.tight_layout()
    out.append(('Componentes xp e yp separadas', _png(g)))

    t = (R['mjd'] - R['mjd'][ev]).astype(float)
    g = Figure(figsize=(7.5, 7)); ax = g.subplots(); ax.plot(R['x'], R['y'], '-', color='gray', lw=.6)
    sc = ax.scatter(R['x'], R['y'], c=t, cmap='coolwarm', s=14, zorder=3); g.colorbar(sc, label='dias desde o evento')
    ax.plot(R['x'][0], R['y'][0], 's', color='green', ms=9, zorder=4, label='início'); ax.plot(R['x'][-1], R['y'][-1], '^', color='black', ms=9, zorder=4, label='fim')
    ax.plot(R['x'][ev], R['y'][ev], '*', color='gold', mec='k', ms=16, zorder=5, label='evento')
    ax.set(xlabel='xp (mas)', ylabel='yp (mas)', title='Trajetória do polo (xp × yp)\n' + tit); ax.set_aspect('equal', 'datalim'); ax.grid(alpha=.3); ax.legend(); g.tight_layout()
    out.append(('Trajetória xp × yp (dispersão)', _png(g)))

    g = Figure(figsize=(11, 5)); ap = g.add_subplot(1, 2, 1, projection='polar'); ac = g.add_subplot(1, 2, 2)
    s1 = ap.scatter(R['th'], R['d'], c=t, cmap='coolwarm', s=12); ap.plot(R['th'][ev], R['d'][ev], '*', color='gold', mec='k', ms=15)
    ap.set_title('Coordenadas polares (θ, d)', fontsize=10); ap.set_rlabel_position(135)
    ac.scatter(R['th'], R['d'], c=t, cmap='coolwarm', s=12); ac.plot(R['th'][ev], R['d'][ev], '*', color='gold', mec='k', ms=15)
    ac.set(xlabel='Ângulo θ (rad)', ylabel='Distância angular d (mas)', title='d × θ'); ac.grid(alpha=.3); g.colorbar(s1, ax=g.axes, shrink=.8, label='dias desde o evento')
    out.append(('Coordenadas polares (θ, d)', _png(g)))

    g = Figure(figsize=(8, 7)); ax = g.add_subplot(111, projection='3d')
    ax.plot(R['x'], R['y'], t, color='gray', lw=.7); ax.scatter(R['x'], R['y'], t, c=t, cmap='coolwarm', s=8)
    ax.plot([R['x'][ev]], [R['y'][ev]], [0], '*', color='gold', mec='k', ms=15, label='evento')
    ax.set(xlabel='xp (mas)', ylabel='yp (mas)', zlabel='tempo (dias desde o evento)', title='Polo no tempo (3D)\n' + tit); ax.legend(); g.tight_layout()
    out.append(('Representação 3D (z = tempo)', _png(g)))

    g = Figure(figsize=(10, 6.5)); a1, a2 = g.subplots(2, 1, sharex=True)
    a1.plot(dt, R['d'], '.', ms=3, color='#1f77b4'); a1.plot(dt, R['md'], '-', color='#ff7f0e', lw=1.8); a1.set_ylabel('d = √(xp²+yp²) (mas)'); datas(a1)
    a1.set_title('Distância do polo ao CIO e variação diária — ' + tit)
    a2.plot(dt, R['var_cm'], '.-', ms=3, lw=.6, color='#2ca02c'); a2.set_ylabel('variação diária (cm)'); a2.set_xlabel('Data (UTC)'); datas(a2)
    a2.text(.01, .95, '1 mas ≈ 3,09 cm na superfície', transform=a2.transAxes, va='top', fontsize=8); g.tight_layout()
    out.append(('Distância angular e variação diária', _png(g)))
    return out