# -*- coding: utf-8 -*-
"""Resultados do IBGE-PPP: lê os arquivos .sum e .pos do zip de resposta, calcula σ3D = √(σφ²+σλ²+σh²), tempo de convergência
e gera os gráficos comparativos do Trabalho Prático (tempo de rastreio × sigma, L1 × L1/L2, curvas época a época)."""
import io, re, unicodedata, zipfile
from datetime import datetime
import numpy as np
from matplotlib.figure import Figure
from nucleo import _png

CORES = {1: '#c62828', 2: '#f9a825', 3: '#2e7d32', None: '#555555'}
NOMES_CEN = {1: 'Cenário 1 (alta obstrução)', 2: 'Cenário 2 (obstrução parcial)', 3: 'Cenário 3 (campo aberto)', None: 'Sem cenário'}


def _sem_acento(t):
    return ''.join(c for c in unicodedata.normalize('NFKD', t) if not unicodedata.combining(c))


def _f(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


# ------------------------------------------------------------------ .sum
def ler_sum(texto):
    t = _sem_acento(texto.replace('\r', ''))
    def campo(rot):
        m = re.search(rot + r'\s*:\s*([^\n]+)', t, re.I)
        return m.group(1).strip() if m else ''
    d = dict(marco=campo(r'Nome do Marco'), modo_op=campo(r'Modo de Opera[a-z]+ do Usuario'), obs_proc=campo(r'Observa[a-z]+ processada'),
             freq=campo(r'Frequencia observada'), orbitas=campo(r'Orbitas dos satelites'), iono=campo(r'Modelo de Ionosfera'),
             elev=_f(campo(r'Angulo de Eleva[a-z]+ \(grau\)')), marco_cfa=_f(campo(r'Distancia Marco->CFA \(m\)')),
             intervalo=_f(campo(r'Intervalo de observa[a-z]+ \(seg\)')), epocas=_f(campo(r'Numero de epocas processadas')),
             satelites=_f(campo(r'Numero de satelites processados')), obs_ok=_f(campo(r'Numero de observa[a-z]+ processadas')),
             obs_rej=_f(campo(r'Numero de observa[a-z]+ rejeitadas')), res_cod=_f(campo(r'Residuos da pseudodistancia')),
             res_fase=_f(campo(r'Residuos da fase da portadora')), antena=(re.search(r'Modelo da Antena\s+(\S+)', t, re.I) or [None, ''])[1])
    for rot, k in ((r'Inicio', 'inicio'), (r'Fim', 'fim')):
        m = re.search(rot + r'\s*:\s*(\d{4})/(\d\d)/(\d\d)\s+(\d\d):(\d\d):([\d.]+)', t)
        d[k] = datetime(int(m[1]), int(m[2]), int(m[3]), int(m[4]), int(m[5]), int(float(m[6]))) if m else None
    i3 = re.search(r'3\.3\s+Coordenadas', t); i4 = re.search(r'3\.4\s+Diferen', t)
    blk = t[i3.start(): i4.start() if i4 and i4.start() > i3.start() else len(t)] if i3 else t
    dms = r'(-?\d+)\s+(\d+)\s+([\d.]+)\s+(-?\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)'
    la = re.search(r'Latitude\s*\(gms\)\s+' + dms, blk); lo = re.search(r'Longitude\s*\(gms\)\s+' + dms, blk)
    al = re.search(r'Alt\.?\s*Geo\.?\s*\(m\)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+([\d.]+)', blk)
    d['s_lat'] = float(la[7]) if la else None; d['s_lon'] = float(lo[7]) if lo else None; d['s_h'] = float(al[3]) if al else None
    d['lat'] = f'{la[1]} {la[2]} {la[3]}' if la else ''; d['lon'] = f'{lo[1]} {lo[2]} {lo[3]}' if lo else ''; d['h'] = _f(al[1]) if al else None
    f = d['freq'].upper()
    d['modo'] = 'L1/L2' if ('L3' in f or '&' in d['iono']) else ('L1' if f else None)
    if d['inicio'] and d['fim']:
        d['dur_min'] = round((d['fim'] - d['inicio']).total_seconds() / 60 + (d['intervalo'] or 0) / 60, 1)
    return d


# ------------------------------------------------------------------ .pos
def ler_pos(texto, intervalo=None):
    """Lê o .pos (época a época) usando os nomes das colunas do cabeçalho (DLAT DLON DHGT … SLAT SLON SHGT). Retorna dict ou None."""
    L = [l for l in texto.replace('\r', '').split('\n') if l.strip()]
    hi = next((i for i, l in enumerate(L) if re.search(r'\bSLAT\b', l, re.I) and re.search(r'\bSHGT\b', l, re.I)), None)
    if hi is None:
        return None
    cols = [c.upper() for c in re.sub(r'^[#%*\s]+', '', L[hi]).split()]
    need = ['DLAT', 'DLON', 'DHGT', 'SLAT', 'SLON', 'SHGT']
    if any(c not in cols for c in need):
        return None
    rows = [t for t in (l.split() for l in L[hi + 1:]) if len(t) >= len(cols) - 4 and _f(t[-1]) is not None]   # linhas de dados: terminam em número
    if not rows:
        return None
    delta = len(rows[0]) - len(cols)                                  # tokens extras antes das colunas numéricas (data/hora)
    g = lambda nome: np.array([_f(r[cols.index(nome) + delta]) if cols.index(nome) + delta < len(r) else np.nan for r in rows], float)
    tk = [next((x for x in r if re.fullmatch(r'\d{1,2}:\d{2}:\d{2}(?:\.\d+)?', x)), None) for r in rows]
    if all(tk):
        sec = np.array([sum(float(p) * m for p, m in zip(x.split(':'), (3600, 60, 1))) for x in tk]); sec = np.unwrap(sec, period=86400)
        t = (sec - sec[0]) / 60
    else:
        t = np.arange(len(rows)) * (intervalo or 1.0) / 60
    return dict(t=t, dlat=g('DLAT'), dlon=g('DLON'), dh=g('DHGT'), slat=g('SLAT'), slon=g('SLON'), sh=g('SHGT'), n=len(rows))


def _tempo_ate(t, ok):
    bad = np.where(~ok)[0]
    if len(bad) == 0:
        return float(t[0])
    return None if bad[-1] == len(ok) - 1 else float(t[bad[-1] + 1])


def convergencia(pos):
    s3 = np.sqrt(pos['slat'] ** 2 + pos['slon'] ** 2 + pos['sh'] ** 2)
    d3 = np.sqrt((pos['dlat'] - pos['dlat'][-1]) ** 2 + (pos['dlon'] - pos['dlon'][-1]) ** 2 + (pos['dh'] - pos['dh'][-1]) ** 2)
    out = {}
    for lim in (1.0, 0.5, 0.1, 0.05):
        out[('sig', lim)] = _tempo_ate(pos['t'], s3 <= lim); out[('dif', lim)] = _tempo_ate(pos['t'], d3 <= lim)
    return out, s3, d3


# ------------------------------------------------------------------ carga dos arquivos enviados
def _stem(n):
    return re.sub(r'\.[A-Za-z0-9]+$', '', n.replace('\\', '/').split('/')[-1])


def inferir(nome, d):
    n = _sem_acento(nome).upper()
    m = re.search(r'(?:CENARIO|CEN|C)[_\- ]?([123])(?!\d)', n)
    cen = int(m[1]) if m else None
    m = re.search(r'(?<![0-9A-Z])(\d+)[_\-]?(MIN|H)(?![A-Z])', n)
    dur = (float(m[1]) * (60 if m[2] == 'H' else 1)) if m else (d.get('dur_min') if d else None)
    modo = (d or {}).get('modo') or ('L1' if re.search(r'(?:^|[_\- ])L1(?:[_\- ]|$)', n) and 'L1L2' not in n else None)
    return cen, dur, modo


def carregar(arquivos):
    """arquivos: lista de (nome, bytes) — zips do IBGE-PPP (ou .sum/.pos soltos). Retorna lista de itens."""
    grupos = {}
    def add(n, b, ns=''):
        if b[:2] == b'PK':
            z = zipfile.ZipFile(io.BytesIO(b))
            for m in z.namelist():
                add(m, z.read(m), _stem(n))
            return
        ext = n.lower().rsplit('.', 1)[-1]
        if ext in ('sum', 'pos'):
            k = _stem(n)
            if k in grupos and ext in grupos[k]:
                k = f'{ns}_{k}' if ns else k
            grupos.setdefault(k, {})[ext] = b.decode('latin-1')
    for n, b in arquivos:
        add(n, b)
    itens = []
    for k, g in sorted(grupos.items()):
        item = dict(nome=k, sum=None, pos=None, erro='')
        try:
            if 'sum' in g:
                item['sum'] = ler_sum(g['sum'])
            if 'pos' in g:
                item['pos'] = ler_pos(g['pos'], (item['sum'] or {}).get('intervalo'))
            if item['sum'] is None or item['sum']['s_lat'] is None:
                item['erro'] = 'não encontrei os sigmas (seção 3.3 do .sum): digite-os na tabela manual'
        except Exception as e:
            item['erro'] = f'erro ao ler: {e}'
        cen, dur, modo = inferir(k, item['sum'])
        item['rot'] = dict(cen=cen, dur=dur, modo=modo)
        itens.append(item)
    return itens


def manual(texto):
    """Linhas 'cenário;duração_min;modo;σφ;σλ;σh[;conv_min]' (vírgula ou ponto decimal)."""
    out = []
    for i, l in enumerate(texto.replace('\r', '').split('\n')):
        if not l.strip() or l.strip().startswith('#'):
            continue
        p = [x.strip().replace(',', '.') for x in l.split(';')]
        if len(p) < 6:
            raise ValueError(f'linha {i + 1}: use cenário;duração(min);modo;σφ;σλ;σh[;convergência(min)]')
        cen = int(float(p[0])) if p[0] else None; modo = 'L1' if p[2].upper().replace(' ', '') == 'L1' else 'L1/L2'
        s = dict(s_lat=float(p[3]), s_lon=float(p[4]), s_h=float(p[5]), modo=modo)
        out.append(dict(nome=f'manual_{i + 1}', sum=s, pos=None, erro='', conv_manual=_f(p[6]) if len(p) > 6 else None,
                        rot=dict(cen=cen, dur=float(p[1]), modo=modo)))
    return out


# ------------------------------------------------------------------ tabela
CAB = ['Arquivo', 'Cenário', 'Duração (min)', 'Modo', 'σφ (m)', 'σλ (m)', 'σh (m)', 'σ3D (m)', 'σ plan. (m)', 'Épocas', 'Satélites', 'Obs. rejeitadas (%)',
       'Res. código (m)', 'Res. fase (cm)', 'Conv. σ3D<1 m (min)', 'Conv. σ3D<0,1 m (min)', 'Antena / altura (m)']


def sigma3d(a, b, c):
    return float(np.sqrt(a ** 2 + b ** 2 + c ** 2))


def linhas(itens):
    out = []
    for it in itens:
        s, r, p = it['sum'] or {}, it['rot'], it['pos']
        ok = s.get('s_lat') is not None and s.get('s_lon') is not None and s.get('s_h') is not None
        c1 = c01 = None
        if p:
            cv, _, _ = convergencia(p); c1, c01 = cv[('sig', 1.0)], cv[('sig', 0.1)]
        elif it.get('conv_manual') is not None:
            c1 = it['conv_manual']
        rej = 100 * s['obs_rej'] / (s['obs_rej'] + s['obs_ok']) if s.get('obs_ok') and s.get('obs_rej') is not None else None
        out.append(dict(nome=it['nome'], cen=r['cen'], dur=r['dur'], modo=r['modo'], sl=s.get('s_lat'), so=s.get('s_lon'), sh=s.get('s_h'),
                        s3=sigma3d(s['s_lat'], s['s_lon'], s['s_h']) if ok else None, sp=float(np.hypot(s['s_lat'], s['s_lon'])) if ok else None,
                        epocas=s.get('epocas'), sats=s.get('satelites'), rej=rej, rc=s.get('res_cod'), rf=s.get('res_fase'), c1=c1, c01=c01,
                        ant=f"{s.get('antena') or '—'} / {s['marco_cfa'] if s.get('marco_cfa') is not None else '—'}"))
    return out


def fmt(v, d=4):
    return '—' if v is None else (f'{v:.{d}f}' if isinstance(v, float) else str(v))


def tabela(L):
    return [[x['nome'], x['cen'] or '—', fmt(x['dur'], 1) if x['dur'] else '—', x['modo'] or '—', fmt(x['sl']), fmt(x['so']), fmt(x['sh']), fmt(x['s3']), fmt(x['sp']),
             fmt(x['epocas'], 0), fmt(x['sats'], 0), fmt(x['rej'], 1), fmt(x['rc'], 2), fmt(x['rf'], 1), fmt(x['c1'], 1), fmt(x['c01'], 1), x['ant']] for x in L]


def csv_ppp(L):
    return '\n'.join([','.join(CAB)] + [','.join(str(c).replace(',', ' ') for c in r) for r in tabela(L)])


# ------------------------------------------------------------------ gráficos
def _rotulo(m):
    return f'{int(m)} min' if m < 60 else (f'{m / 60:g} h')


def graficos(L, idx_item=None, itens=None):
    out = []
    ok = [x for x in L if x['s3'] is not None and x['dur']]
    modos = [m for m in ('L1/L2', 'L1') if any(x['modo'] == m for x in ok)] or [None]
    for modo in modos:
        sel = [x for x in ok if x['modo'] == modo]
        if not sel:
            continue
        durs = sorted({x['dur'] for x in sel}); pos = {d: i for i, d in enumerate(durs)}
        g = Figure(figsize=(11, 7.5)); ax = g.subplots(2, 2, sharex=True)
        for a, (k, tit) in zip(ax.ravel(), (('sl', 'σφ (latitude)'), ('so', 'σλ (longitude)'), ('sh', 'σh (altitude)'), ('s3', 'σ3D'))):
            for cen in sorted({x['cen'] for x in sel}, key=lambda c: (c is None, c)):
                pts = sorted((pos[x['dur']], x[k]) for x in sel if x['cen'] == cen and x[k] is not None)
                if pts:
                    a.plot(*zip(*pts), '-o', color=CORES[cen], label=NOMES_CEN[cen])
            a.set(title=tit, ylabel='sigma (m)'); a.grid(alpha=.3, which='both')
            vals = [x[k] for x in sel if x[k]]
            if vals and max(vals) / max(min(vals), 1e-9) > 20:
                a.set_yscale('log')
        ax[1, 0].set_xticks(range(len(durs))); ax[1, 0].set_xticklabels([_rotulo(d) for d in durs]); ax[1, 1].set_xticklabels([_rotulo(d) for d in durs])
        ax[1, 0].set_xlabel('Tempo de rastreio'); ax[1, 1].set_xlabel('Tempo de rastreio'); ax[0, 0].legend(fontsize=8)
        g.suptitle(f'Tempo de rastreio × sigma — {"modo " + modo if modo else "todos os arquivos"}'); g.tight_layout()
        out.append((f'Tempo de rastreio × sigmas ({modo or "todos"})', _png(g)))
    if all(any(x['modo'] == m for x in ok) for m in ('L1', 'L1/L2')):
        cens = sorted({x['cen'] for x in ok}, key=lambda c: (c is None, c)); durs = sorted({x['dur'] for x in ok})
        g = Figure(figsize=(4.2 * len(cens) + 1, 4.6)); axs = g.subplots(1, len(cens), sharey=True, squeeze=False)[0]
        for a, cen in zip(axs, cens):
            for j, (modo, cor) in enumerate((('L1/L2', '#1565c0'), ('L1', '#ef6c00'))):
                v = [next((x['s3'] for x in ok if x['cen'] == cen and x['dur'] == d and x['modo'] == modo), np.nan) for d in durs]
                a.bar(np.arange(len(durs)) + (j - .5) * .38, v, .38, color=cor, label=modo)
            a.set_xticks(range(len(durs))); a.set_xticklabels([_rotulo(d) for d in durs]); a.set_title(NOMES_CEN[cen], fontsize=9); a.grid(axis='y', alpha=.3, which='both'); a.set_yscale('log')
        axs[0].set_ylabel('σ3D (m)'); axs[0].legend(); g.suptitle('Solução L1 pura × dupla frequência (L1/L2): σ3D'); g.tight_layout()
        out.append(('L1 pura × L1/L2 (σ3D)', _png(g)))
    if itens is not None and idx_item is not None and itens[idx_item]['pos']:
        p = itens[idx_item]['pos']; cv, s3, d3 = convergencia(p); nome = itens[idx_item]['nome']
        g = Figure(figsize=(10, 8)); a1, a2 = g.subplots(2, 1, sharex=True)
        a1.plot(p['t'], p['slat'], label='σφ (SLAT)', color='#c62828'); a1.plot(p['t'], p['slon'], label='σλ (SLON)', color='#1565c0'); a1.plot(p['t'], p['sh'], label='σh (SHGT)', color='#2e7d32')
        a1.plot(p['t'], s3, 'k--', label='σ3D'); a1.set_yscale('log'); a1.axhline(1.0, color='gray', lw=.8); a1.axhline(0.1, color='gray', lw=.8, ls=':')
        a1.set(ylabel='desvio padrão (m)', title=f'Desvio padrão época a época — {nome}'); a1.legend(fontsize=8); a1.grid(alpha=.3, which='both')
        a2.plot(p['t'], p['dlat'] - p['dlat'][-1], color='#c62828', label='ΔLAT'); a2.plot(p['t'], p['dlon'] - p['dlon'][-1], color='#1565c0', label='ΔLON')
        a2.plot(p['t'], p['dh'] - p['dh'][-1], color='#2e7d32', label='ΔALT'); a2.axhline(0, color='k', lw=.6)
        a2.set(xlabel='Tempo de rastreio (min)', ylabel='diferença para a solução final (m)', title='Convergência das coordenadas'); a2.legend(fontsize=8); a2.grid(alpha=.3)
        g.tight_layout(); out.append((f'Convergência e desvio padrão — {nome}', _png(g)))
    return out