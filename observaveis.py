# -*- coding: utf-8 -*-
"""Qualidade das observáveis GPS L1/L2 (RINEX de observação da RBMC): combinações Lc, Lg, Pc,
saltos de ciclo, disponibilidade. Aceita .zip do IBGE, .gz, CRINEX (Hatanaka) e RINEX 2/3."""
import gzip, io, zipfile
from datetime import datetime, timedelta
import numpy as np
import matplotlib.dates as mdates
from matplotlib.figure import Figure
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import crinex
from nucleo import _png

F1, F2, C = 1575.42e6, 1227.60e6, 299792458.0
PREF = (['L1', 'L1C', 'L1W', 'L1P', 'L1X'], ['L2', 'L2W', 'L2P', 'L2C', 'L2X', 'L2L'],
        ['P1', 'C1W', 'C1', 'C1C', 'C1P'], ['P2', 'C2W', 'C2P', 'C2', 'C2C', 'C2X', 'C2L'],
        ['S1', 'S1C', 'S1W', 'S1P', 'S1X'], ['S2', 'S2W', 'S2P', 'S2C', 'S2X', 'S2L'])


def _eh_obs(b):
    cab = b[:400].decode('latin-1', 'ignore')
    return 'CRINEX' in cab or 'OBSERVATION DATA' in cab


def _achar(dados):
    """Procura o arquivo de observação, entrando em .zip/.gz aninhados (ex.: zip do portal com o zip da estação)."""
    if dados[:2] == b'\x1f\x8b':
        dados = gzip.decompress(dados)
    if dados[:2] == b'PK':
        z = zipfile.ZipFile(io.BytesIO(dados))
        for n in z.namelist():
            r = _achar(z.read(n))
            if r is not None:
                return r
        return None
    return dados if _eh_obs(dados) else None


def extrair(dados):
    """Devolve o texto do RINEX de observação (descompacta gz, zip, zip aninhado e Hatanaka)."""
    dados = _achar(dados)
    if dados is None:
        raise ValueError('não encontrei arquivo de observação (RINEX/CRINEX) no arquivo enviado')
    if b'CRINEX' in dados[:100]:
        try:
            return crinex.crx2rnx(dados).replace('\r', '')          # Python puro (não depende de executáveis)
        except Exception as e:
            try:                                                      # alternativa, se existir a biblioteca hatanaka
                import hatanaka
                dados = hatanaka.decompress(dados)
            except Exception:
                raise ValueError(f'não consegui descompactar o CRINEX (Hatanaka): {e}')
    return dados.decode('latin-1').replace('\r', '')


def _num(s):
    s = s.strip()
    try:
        return float(s)
    except ValueError:
        return np.nan


def ler_obs(txt):
    L = txt.split('\n')
    v3 = float(L[0][:9]) >= 3
    tipos, marker, i, ult = {}, '', 0, None
    xyz, rec, ant = None, '', ''
    while 'END OF HEADER' not in L[i]:
        l, lab = L[i], L[i][60:].strip()
        if lab == 'MARKER NAME':
            marker = l[:60].strip()
        elif lab == 'APPROX POSITION XYZ':
            xyz = [float(l[k:k + 14]) for k in (0, 14, 28)]
        elif lab == 'REC # / TYPE / VERS':
            rec = l[20:40].strip()
        elif lab == 'ANT # / TYPE':
            ant = l[20:40].strip()
        elif v3 and lab == 'SYS / # / OBS TYPES':
            if l[0] != ' ':
                ult = l[0]; tipos[ult] = []
            tipos[ult] += l[7:60].split()
        elif not v3 and lab == '# / TYPES OF OBSERV':
            tipos.setdefault('G', []); tipos['G'] += l[6:60].split()
        i += 1
    i += 1
    tg = tipos['G']
    cands = [[tg.index(t) for t in p if t in tg] for p in PREF]
    cols = sorted({k for c_ in cands for k in c_})
    ep, regs = [], []

    def vals(texto):
        return [_num(texto[16 * k:16 * k + 14]) for k in cols]
    while i < len(L):
        l = L[i]
        if v3:
            if not l.startswith('>'):
                i += 1; continue
            t = datetime(int(l[2:6]), int(l[7:9]), int(l[10:12]), int(l[13:15]), int(l[16:18])) + timedelta(seconds=float(l[18:29]))
            flag, ns = int(l[31]), int(l[32:35]); i += 1
            if flag > 1:
                i += ns; continue
            ei = len(ep); ep.append(t)
            for _ in range(ns):
                s = L[i]; i += 1
                if s[0] == 'G':
                    regs.append((ei, s[:3].replace(' ', '0'), vals(s[3:])))
        else:
            if len(l) < 32:
                i += 1; continue
            yy = int(l[1:3])
            t = datetime(2000 + yy if yy < 80 else 1900 + yy, int(l[4:6]), int(l[7:9]), int(l[10:12]), int(l[13:15])) + timedelta(seconds=float(l[15:26]))
            flag, ns = int(l[28:29]), int(l[29:32]); i += 1
            sv = l[32:68]
            while len(sv) < 3 * ns:
                sv += L[i][32:68]; i += 1
            if flag > 1:
                i += ns; continue
            ei = len(ep); ep.append(t)
            nl = -(-len(tg) // 5)
            for k in range(ns):
                s = sv[3 * k:3 * k + 3]
                txt_s = ''.join(L[i + j].ljust(80)[:80] for j in range(nl)); i += nl
                if s[0] in ' G':
                    regs.append((ei, 'G' + s[1:3].replace(' ', '0'), vals(txt_s)))
    sats = sorted({r[1] for r in regs})
    A = np.full((len(ep), len(sats), len(cols)), np.nan)
    ix = {s: k for k, s in enumerate(sats)}
    for ei, s, v in regs:
        A[ei, ix[s]] = v
    # por grupo (L1, L2, código 1, código 2) usa o tipo que realmente tem mais dados (empate: ordem de preferência)
    cont = np.isfinite(A).sum((0, 1))
    esc = [max(c_, key=lambda k: cont[cols.index(k)]) if c_ else None for c_ in cands]
    A = np.stack([np.full(A.shape[:2], np.nan) if k is None else A[:, :, cols.index(k)] for k in esc], axis=2)
    t0 = ep[0]
    return dict(estacao=marker[:4].upper(), ep=ep, t=np.array([(x - t0).total_seconds() for x in ep]), sats=sats, A=A,
                obs=[None if k is None else tg[k] for k in esc], xyz=xyz, rec=rec, ant=ant)


def carregar(dados, nome):
    f = ler_obs(extrair(dados))
    f['nome'] = nome
    f['estacao'] = f['estacao'] or nome[:4].upper()
    return f


def combinar(f, limiar):
    L1, L2, P1, P2 = (f['A'][:, :, k] for k in range(4))
    a, b = L1 * C / F1, L2 * C / F2
    k = F1 ** 2 - F2 ** 2
    Lg, Lc, Pc = a - b, (F1 ** 2 * a - F2 ** 2 * b) / k, (F1 ** 2 * P1 - F2 ** 2 * P2) / k
    d = np.full_like(Lg, np.nan); d[1:] = Lg[1:] - Lg[:-1]
    dt = np.diff(f['t']); passo = np.median(dt) if len(dt) else 30
    d[1:][dt > 2.5 * passo] = np.nan                       # lacuna = novo arco, não é salto
    return dict(Lg=Lg, Lc=Lc, Pc=Pc, d=d, slip=np.abs(d) > limiar, passo=passo)


def resumo(f, limiar):
    c = combinar(f, limiar); A = f['A']
    n1, n12 = np.isfinite(A[:, :, 0]).sum(), np.isfinite(A[:, :, 1]).sum()
    nc, ns = np.isfinite(c['d']).sum(), c['slip'].sum()
    return [f['estacao'], f['ep'][0].strftime('%d/%m/%Y'), len(f['ep']), '%g' % c['passo'], len(f['sats']), int(n1),
            '%.1f' % (100 * n12 / max(n1, 1)), int(ns), '%.2f' % (100 * ns / max(nc, 1))]


RES_CAB = ['Estação', 'Data', 'Épocas', 'Intervalo (s)', 'Satélites GPS', 'Obs. L1', 'L2 válidas / L1 (%)', 'Saltos', 'Taxa de saltos (%)']


def csv_dados(f, limiar):
    c = combinar(f, limiar); A = f['A']
    out = ['época,satélite,L1_ciclos,L2_ciclos,cod1_m,cod2_m,Lc_m,Lg_m,Pc_m,dLg_m,salto']
    for e in range(len(f['ep'])):
        for s in range(len(f['sats'])):
            if np.isfinite(A[e, s, 0]) or np.isfinite(A[e, s, 1]):
                v = [A[e, s, 0], A[e, s, 1], A[e, s, 2], A[e, s, 3], c['Lc'][e, s], c['Lg'][e, s], c['Pc'][e, s], c['d'][e, s]]
                out.append('%s,%s,%s,%d' % (f['ep'][e].strftime('%Y-%m-%d %H:%M:%S'), f['sats'][s],
                                            ','.join('' if np.isnan(x) else '%.4f' % x for x in v), int(c['slip'][e, s])))
    return '\n'.join(out)


def _tempo(ax, ep):
    fmt = '%H:%M' if (ep[-1] - ep[0]).total_seconds() <= 86400 else '%d/%m\n%H:%M'
    ax.xaxis.set_major_formatter(mdates.DateFormatter(fmt))


def graficos(fs, idx, sat, limiar):
    """Retorna [(título, png)]."""
    f = fs[idx]; c = combinar(f, limiar); ep = f['ep']; sats = f['sats']
    si = sats.index(sat); tit = f"{f['estacao']} · {ep[0]:%d/%m/%Y}"
    out = []
    g = Figure(figsize=(10, 6.5)); a1, a2 = g.subplots(2, 1, sharex=True)
    lg = c['Lg'][:, si]; ok = np.isfinite(lg); base = lg[ok][0] if ok.any() else 0.0
    a1.plot(np.array(ep)[ok], lg[ok] - base, '.', ms=2, color='#1f77b4')
    sl = c['slip'][:, si]; a1.plot(np.array(ep)[sl], lg[sl] - base, 'rx', ms=7, label='possível salto de ciclo')
    a1.set(ylabel='Lg − Lg₀ (m)', title=f'Observáveis do {sat} — {tit}'); a1.legend(loc='upper right'); a1.grid(alpha=.3)
    d = c['Lc'][:, si] - c['Pc'][:, si]; ok2 = np.isfinite(d)
    a2.plot(np.array(ep)[ok2], d[ok2], '.', ms=2, color='#2ca02c'); a2.set(ylabel='Lc − Pc (m)', xlabel='Hora (GPS)'); a2.grid(alpha=.3)
    _tempo(a2, ep); g.tight_layout(); out.append((f'Séries de Lg e Lc−Pc do {sat}', _png(g)))

    g = Figure(figsize=(10, 4)); ax = g.subplots()
    ax.plot(ep, np.isfinite(f['A'][:, :, 0]).sum(1), color='#1f77b4', label='L1')
    ax.plot(ep, (np.isfinite(f['A'][:, :, 0]) & np.isfinite(f['A'][:, :, 1])).sum(1), color='#d62728', label='L1 e L2')
    ax.set(ylabel='Nº de satélites', xlabel='Hora (GPS)', title=f'Satélites rastreados por época — {tit}'); ax.grid(alpha=.3); ax.legend(); _tempo(ax, ep)
    g.tight_layout(); out.append(('Satélites rastreados por época', _png(g)))

    A = f['A']; n_ok = np.isfinite(A[:, :, 0]).astype(int) + np.isfinite(A[:, :, 1]).astype(int)
    g = Figure(figsize=(10, 6)); ax = g.subplots()
    ax.pcolormesh(mdates.date2num(ep), np.arange(len(sats)), n_ok.T, cmap=ListedColormap(['white', '#ffb347', '#2e8b57']), vmin=0, vmax=2, shading='nearest')
    ax.set_yticks(range(len(sats))); ax.set_yticklabels(sats, fontsize=7); ax.set(xlabel='Hora (GPS)', title=f'Disponibilidade de dados — {tit}'); _tempo(ax, ep)
    ax.legend(handles=[Patch(color='#2e8b57', label='L1 e L2'), Patch(color='#ffb347', label='só uma portadora'), Patch(fc='white', ec='gray', label='sem dado')], loc='upper right', fontsize=7)
    g.tight_layout(); out.append(('Mapa de disponibilidade', _png(g)))

    g = Figure(figsize=(10, 7)); a1, a2 = g.subplots(2, 1, gridspec_kw={'height_ratios': [2, 1]})
    E = np.repeat(np.array(ep)[:, None], len(sats), 1); d = c['d']; ok = np.isfinite(d)
    a1.plot(E[ok & ~c['slip']], d[ok & ~c['slip']], '.', ms=1.5, color='gray'); a1.plot(E[c['slip']], d[c['slip']], 'r.', ms=4, label='acima do limiar')
    for s in (-limiar, limiar): a1.axhline(s, color='r', lw=.8, ls='--')
    a1.set(ylim=(-6 * limiar, 6 * limiar), ylabel='ΔLg entre épocas (m)', title=f'Diagnóstico de saltos de ciclo (limiar {limiar:g} m) — {tit}'); a1.legend(); a1.grid(alpha=.3); _tempo(a1, ep)
    taxa = 100 * c['slip'].sum(0) / np.maximum(np.isfinite(d).sum(0), 1)
    a2.bar(sats, taxa, color='#d62728'); a2.set(ylabel='Saltos (%)'); a2.tick_params(axis='x', labelsize=6, rotation=90); a2.grid(axis='y', alpha=.3)
    g.tight_layout(); out.append(('Diagnóstico de saltos de ciclo', _png(g)))

    return out


# ---------------------------------------------------------------- comparação entre duas estações
CAB_COMP = ['Estação', 'Dias', 'Épocas', 'Satélites/época (média)', 'L2 válidas / L1 (%)', 'Saltos', 'Taxa de saltos (%)']
COR = ('#1f77b4', '#d62728')
pct = lambda a, b: 100 * a / b if b else np.nan


def metricas(f, limiar):
    c = combinar(f, limiar); A = f['A']
    f1, f2 = np.isfinite(A[:, :, 0]), np.isfinite(A[:, :, 1]); d = np.isfinite(c['d']); sl = c['slip']
    h = np.array([e.hour for e in f['ep']])
    return dict(data=f['ep'][0].date(), nep=len(f['ep']), n1=int(f1.sum()), n12=int(f2.sum()), nc=int(d.sum()), ns=int(sl.sum()),
                sat={s: np.array([f1[:, i].sum(), f2[:, i].sum(), d[:, i].sum(), sl[:, i].sum()]) for i, s in enumerate(f['sats'])},
                hora=np.array([[f1[h == k].sum(), f2[h == k].sum(), d[h == k].sum(), sl[h == k].sum()] for k in range(24)]),
                nep_h=np.array([(h == k).sum() for k in range(24)]))


def _soma(ms):
    t = dict(dias=len({m['data'] for m in ms}), nep=sum(m['nep'] for m in ms), n1=sum(m['n1'] for m in ms), n12=sum(m['n12'] for m in ms),
             nc=sum(m['nc'] for m in ms), ns=sum(m['ns'] for m in ms), hora=sum(m['hora'] for m in ms), nep_h=sum(m['nep_h'] for m in ms), sat={})
    for m in ms:
        for s, v in m['sat'].items():
            t['sat'][s] = t['sat'].get(s, 0) + v
    return t


def comparar(fs, ea, eb, limiar):
    """Compara duas estações (agregando todos os arquivos/dias de cada uma). Retorna (cabeçalho, linhas, [(título, png)])."""
    est = (ea, eb)
    ms = {e: [metricas(f, limiar) for f in fs if f['estacao'] == e] for e in est}
    T = {e: _soma(ms[e]) for e in est}
    val = lambda t: (t['n1'] / max(t['nep'], 1), pct(t['n12'], t['n1']), pct(t['ns'], t['nc']))
    rows = [[e, T[e]['dias'], T[e]['nep'], '%.2f' % val(T[e])[0], '%.1f' % val(T[e])[1], T[e]['ns'], '%.2f' % val(T[e])[2]] for e in est]
    va, vb = val(T[ea]), val(T[eb])
    rows.append([f'{ea} − {eb}', '', '', '%+.2f' % (va[0] - vb[0]), '%+.1f' % (va[1] - vb[1]), '', '%+.2f' % (va[2] - vb[2])])
    out = []
    g = Figure(figsize=(10, 4)); ax = g.subplots(1, 3)
    for k, (tt, yl) in enumerate((('Satélites por época', 'nº médio'), ('L2 válidas / L1', '%'), ('Taxa de possíveis saltos', '%'))):
        ax[k].bar(est, [val(T[e])[k] for e in est], color=COR); ax[k].set(title=tt, ylabel=yl); ax[k].grid(axis='y', alpha=.3)
    g.suptitle(f'{ea} × {eb} — visão geral'); g.tight_layout(); out.append(('Visão geral', _png(g)))

    g = Figure(figsize=(10, 8)); ax = g.subplots(3, 1, sharex=True); hrs = np.arange(24)
    for e, cor in zip(est, COR):
        H = T[e]['hora']
        ax[0].plot(hrs, H[:, 0] / np.maximum(T[e]['nep_h'], 1), '-o', ms=3, color=cor, label=e)
        ax[1].plot(hrs, [pct(a, b) for a, b in zip(H[:, 1], H[:, 0])], '-o', ms=3, color=cor, label=e)
        ax[2].plot(hrs, [pct(a, b) for a, b in zip(H[:, 3], H[:, 2])], '-o', ms=3, color=cor, label=e)
    for a, yl in zip(ax, ('Satélites por época', 'L2 válidas / L1 (%)', 'Taxa de saltos (%)')):
        a.set_ylabel(yl); a.grid(alpha=.3); a.legend()
    ax[2].set(xlabel='Hora do dia (tempo GPS)', xticks=range(0, 24, 2)); ax[0].set_title('Comportamento ao longo do dia (média dos dias enviados)')
    g.tight_layout(); out.append(('Variação ao longo do dia', _png(g)))

    sats = sorted(set(T[ea]['sat']) | set(T[eb]['sat']))
    g = Figure(figsize=(11, 6)); ax = g.subplots(2, 1, sharex=True); x = np.arange(len(sats))
    for k, (e, cor) in enumerate(zip(est, COR)):
        v = [T[e]['sat'].get(s, np.zeros(4)) for s in sats]
        ax[0].bar(x + (k - .5) * .4, [0 if np.isnan(pct(a[1], a[0])) else pct(a[1], a[0]) for a in v], .4, color=cor, label=e)
        ax[1].bar(x + (k - .5) * .4, [0 if np.isnan(pct(a[3], a[2])) else pct(a[3], a[2]) for a in v], .4, color=cor, label=e)
    ax[0].set_ylabel('L2 válidas / L1 (%)'); ax[1].set_ylabel('Taxa de saltos (%)'); ax[0].legend(); ax[0].set_title('Desempenho por satélite')
    ax[1].set_xticks(x); ax[1].set_xticklabels(sats, rotation=90, fontsize=7)
    for a in ax: a.grid(axis='y', alpha=.3)
    g.tight_layout(); out.append(('Desempenho por satélite', _png(g)))

    if all(len({m['data'] for m in ms[e]}) >= 2 for e in est):
        g = Figure(figsize=(10, 6)); ax = g.subplots(2, 1, sharex=True)
        for e, cor in zip(est, COR):
            M = sorted(ms[e], key=lambda m: m['data'])
            ax[0].plot([m['data'] for m in M], [pct(m['n12'], m['n1']) for m in M], '-o', color=cor, label=e)
            ax[1].plot([m['data'] for m in M], [pct(m['ns'], m['nc']) for m in M], '-o', color=cor, label=e)
        ax[0].set_ylabel('L2 válidas / L1 (%)'); ax[1].set_ylabel('Taxa de saltos (%)'); ax[0].legend(); ax[0].set_title('Evolução dia a dia')
        for a in ax: a.grid(alpha=.3)
        ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%d/%m')); g.tight_layout(); out.append(('Evolução dia a dia', _png(g)))

    comuns = sorted({f['ep'][0].date() for f in fs if f['estacao'] == ea} & {f['ep'][0].date() for f in fs if f['estacao'] == eb})
    if comuns:
        g = Figure(figsize=(10, 4)); ax = g.subplots()
        for e, cor in zip(est, COR):
            f = next(x for x in fs if x['estacao'] == e and x['ep'][0].date() == comuns[0])
            ax.plot(f['ep'], np.isfinite(f['A'][:, :, 0]).sum(1), color=cor, label=e, lw=1)
        ax.set(ylabel='Nº de satélites (L1)', xlabel='Hora (GPS)', title=f'Satélites rastreados em {comuns[0]:%d/%m/%Y}'); ax.grid(alpha=.3); ax.legend()
        _tempo(ax, f['ep']); g.tight_layout(); out.append(('Satélites por época no mesmo dia', _png(g)))
    return CAB_COMP, rows, out