# -*- coding: utf-8 -*-
"""Edição de RINEX 2.xx de observação, equivalente a partes do teqc:
  - corte de sessões (teqc -st ... -e ...)
  - filtro de constelação (teqc -R -E -C -S)
  - seleção de observáveis (teqc -O.obs L1,C1)
  - correção do modelo e da altura da antena no cabeçalho.

Interface usada pelo app.py:
  nome_seguro(texto) -> str
  ler(txt) -> dict(cab, tipos, ep, versao)      ep = [dict(t, flag, clk, sats, obs)]
  np_median(lista) -> float
  sessoes(txt, pref, ext, inicio, durs, completa, l1, dupla, sist, tipos, antena, altura)
      -> gera (nome_arquivo, texto_rinex, estatisticas, descricao)

Épocas com flag de evento (2 a 6) são descartadas nas saídas.
"""
import math
import re
import statistics
import unicodedata
from datetime import datetime, timedelta

POR_LINHA = 5  # observáveis por linha no registro de dados (RINEX 2)


# ------------------------------------------------------------------ utilitários
def nome_seguro(texto):
    """Nome de arquivo sem acentos, espaços ou caracteres especiais (ex.: 'Cen 1' -> 'CEN_1')."""
    s = unicodedata.normalize('NFKD', str(texto or '')).encode('ascii', 'ignore').decode()
    s = re.sub(r'[^A-Za-z0-9_\-]+', '_', s).strip('_-')
    s = re.sub(r'_+', '_', s)
    return s.upper() or 'RINEX'


def np_median(valores):
    v = [x for x in valores if x is not None]
    return float(statistics.median(v)) if v else 0.0


def _rot(rotulo, conteudo=''):
    """Linha de cabeçalho: conteúdo nas colunas 1-60 e rótulo nas 61-80."""
    return conteudo[:60].ljust(60) + rotulo


def _sat(tok):
    """'G01', 'G 1', ' 1' -> 'G01'."""
    sis = tok[0] if tok[0] != ' ' else 'G'
    return f'{sis}{int(tok[1:]):02d}'


def _rotulo_dur(m):
    if m >= 60 and abs(m / 60 - round(m / 60)) < 1e-9:
        return f'{int(round(m / 60))}H'
    return (f'{m:g}'.replace('.', 'p')) + 'MIN'


def _texto_dur(m):
    if m >= 60 and abs(m / 60 - round(m / 60)) < 1e-9:
        return f'{int(round(m / 60))} h'
    return f'{m:g} min'.replace('.', ',')


# ------------------------------------------------------------------ leitura
def ler(txt):
    linhas = txt.replace('\r\n', '\n').replace('\r', '\n').split('\n')
    cab, tipos, i = [], [], 0
    while i < len(linhas):
        l = linhas[i].rstrip()
        cab.append(l)
        i += 1
        rot = l.ljust(80)[60:].strip()
        if rot == '# / TYPES OF OBSERV':
            tipos += l[6:60].split()
        elif rot == 'END OF HEADER':
            break
    else:
        raise ValueError('cabeçalho do RINEX sem END OF HEADER')
    versao = cab[0][:9].strip() if cab else ''
    try:
        if not 2 <= float(versao) < 3:
            raise ValueError(f'só RINEX 2.xx é suportado (o arquivo é versão {versao})')
    except ValueError as e:
        if 'suportado' in str(e):
            raise
    if not tipos:
        raise ValueError('cabeçalho sem "# / TYPES OF OBSERV"')

    nl = math.ceil(len(tipos) / POR_LINHA)
    ep = []
    while i < len(linhas):
        l = linhas[i].ljust(80)
        if not l.strip():
            i += 1
            continue
        try:
            flag = l[28].strip() or '0'
            n = int(l[29:32])
        except (ValueError, IndexError):
            raise ValueError(f'linha de época inválida (linha {i + 1}): {linhas[i].strip()[:40]}')
        if flag in '2345':          # registros especiais / mudança de cabeçalho
            i += 1 + n
            continue
        sats_txt = l[32:68]
        clk = l[68:80].strip()
        for k in range(1, math.ceil(n / 12)):
            sats_txt += linhas[i + k].ljust(80)[32:68]
        i += max(1, math.ceil(n / 12))
        if flag == '6':             # registros de cycle slip
            i += n * nl
            continue
        yy = int(l[1:3])
        t = datetime(yy + (2000 if yy < 80 else 1900), int(l[4:6]), int(l[7:9]), int(l[10:12]), int(l[13:15])) \
            + timedelta(seconds=float(l[15:26]))
        sats = [_sat(sats_txt[3 * k:3 * k + 3]) for k in range(n)]
        obs = {}
        for s in sats:
            campos = ''
            for k in range(nl):
                campos += (linhas[i + k] if i + k < len(linhas) else '').ljust(80)[:80]
            obs[s] = [campos[16 * j:16 * j + 16] for j in range(len(tipos))]
            i += nl
        ep.append(dict(t=t, flag=flag, clk=clk, sats=sats, obs=obs))
    if not ep:
        raise ValueError('nenhuma época de observação no arquivo')
    return dict(cab=cab, tipos=tipos, ep=ep, versao=versao)


# ------------------------------------------------------------------ escrita
def _linha_tempo(t, sistema):
    s = t.second + t.microsecond / 1e6
    return f'{t.year:6d}{t.month:6d}{t.day:6d}{t.hour:6d}{t.minute:6d}{s:13.7f}     {sistema:3s}'


def _antena(antena):
    partes = antena.split()
    if len(partes) >= 2 and len(partes[-1]) == 4 and partes[-1].isalnum():
        modelo, radome = ' '.join(partes[:-1]), partes[-1].upper()
    else:
        modelo, radome = ' '.join(partes), 'NONE'
    return modelo.upper()[:16].ljust(16) + radome[:4]


def _cabecalho(cab, tipos, t0, t1, sist_txt, antena, altura, comentario):
    out, tipos_feito = [], False
    for l in cab:
        l80 = l.ljust(80)
        rot = l80[60:].strip()
        if rot == '# / TYPES OF OBSERV':
            if not tipos_feito:
                for k in range(0, len(tipos), 9):
                    bloco = ''.join(f'{x:>6s}' for x in tipos[k:k + 9])
                    out.append(_rot('# / TYPES OF OBSERV', (f'{len(tipos):6d}' if k == 0 else ' ' * 6) + bloco))
                tipos_feito = True
            continue
        if rot in ('# OF SATELLITES', 'PRN / # OF OBS'):   # deixariam de valer após o corte
            continue
        if rot == 'RINEX VERSION / TYPE' and sist_txt:
            l80 = l80[:40] + {'G': 'G (GPS)', 'M': 'M (MIXED)'}.get(sist_txt, sist_txt).ljust(20) + l80[60:]
        elif rot == 'TIME OF FIRST OBS':
            l80 = _rot(rot, _linha_tempo(t0, l80[48:51].strip() or 'GPS'))
        elif rot == 'TIME OF LAST OBS':
            l80 = _rot(rot, _linha_tempo(t1, l80[48:51].strip() or 'GPS'))
        elif rot == 'ANT # / TYPE' and antena:
            l80 = _rot(rot, l80[:20] + _antena(antena))
        elif rot == 'ANTENNA: DELTA H/E/N' and altura is not None:
            l80 = _rot(rot, f'{altura:14.4f}' + l80[14:60])
        elif rot == 'END OF HEADER':
            out.append(_rot('COMMENT', comentario))
        out.append(l80.rstrip())
    return out


def _linha_epoca(t, flag, sats, clk):
    s = t.second + t.microsecond / 1e6
    ini = f' {t.year % 100:02d} {t.month:2d} {t.day:2d} {t.hour:2d} {t.minute:2d}{s:11.7f}  {flag}{len(sats):3d}'
    linhas = []
    for k in range(0, len(sats), 12):
        corpo = ''.join(sats[k:k + 12])
        if k == 0:
            l = ini + corpo
            if clk:
                l = l.ljust(68) + f'{float(clk):12.9f}'
        else:
            l = ' ' * 32 + corpo
        linhas.append(l)
    return linhas


def _gerar(R, t0, t1, sist, sel, comentario, antena, altura):
    """sist = None (todas) ou string de sistemas; sel = None (todas) ou lista de observáveis."""
    tipos = R['tipos'] if sel is None else sel
    idx = [R['tipos'].index(x) for x in tipos]
    corpo, n_ep, todos_sats, prim, ult = [], 0, set(), None, None
    for e in R['ep']:
        if not (t0 <= e['t'] <= t1):
            continue
        sats = []
        for s in e['sats']:
            if sist and s[0] not in sist:
                continue
            vals = [e['obs'][s][j] for j in idx]
            if sel is not None and all(not v[:14].strip() for v in vals):
                continue
            sats.append((s, vals))
        if not sats:
            continue
        corpo += _linha_epoca(e['t'], e['flag'], [s for s, _ in sats], e['clk'])
        for _, vals in sats:
            for k in range(0, len(vals), POR_LINHA):
                corpo.append(''.join(vals[k:k + POR_LINHA]).rstrip())
        todos_sats.update(s for s, _ in sats)
        n_ep += 1
        prim = prim or e['t']
        ult = e['t']
    if not n_ep:
        return None
    sist_txt = (sist if len(sist) == 1 else 'M') if sist else None
    cab = _cabecalho(R['cab'], tipos, prim, ult, sist_txt, antena, altura, comentario)
    est = dict(epocas=n_ep, inicio=prim, fim=ult, satelites=len(todos_sats), tipos=list(tipos))
    return '\n'.join(cab + corpo) + '\n', est


def sessoes(txt, pref, ext, inicio, durs, completa, l1, dupla, sist, tipos, antena=None, altura=None):
    R = ler(txt)
    t_ini = inicio or R['ep'][0]['t']
    if t_ini > R['ep'][-1]['t']:
        raise ValueError(f'o início escolhido ({t_ini:%d/%m/%Y %H:%M:%S}) é depois da última época do arquivo '
                         f'({R["ep"][-1]["t"]:%d/%m/%Y %H:%M:%S})')
    sist = ''.join(sorted(set(sist or 'G')))
    nome_sis = 'GPS' if sist == 'G' else sist
    disp = ' '.join(R['tipos'])

    if l1:
        faltam = [x for x in tipos if x not in R['tipos']]
        if faltam:
            raise ValueError(f'observável {", ".join(faltam)} não existe no arquivo (disponíveis: {disp})')
    if dupla:
        sel2 = [x for x in ('L1', 'L2', 'C1', 'P1', 'P2', 'C2') if x in R['tipos']]
        if 'L1' not in sel2 or 'L2' not in sel2:
            raise ValueError(f'o arquivo não tem L1 e L2 para a versão de dupla frequência (disponíveis: {disp})')

    for d in durs:
        if d <= 0:
            raise ValueError('as durações precisam ser maiores que zero')
        t_fim = t_ini + timedelta(minutes=d)
        lab, txt_d = _rotulo_dur(d), _texto_dur(d)
        aviso = ' (o arquivo termina antes do fim da sessão)' if t_fim > R['ep'][-1]['t'] + timedelta(seconds=1) else ''
        variantes = []
        if completa:
            variantes.append((f'{pref}_{lab}.{ext}', None, None, f'Sessão de {txt_d}, sem filtro'))
        if l1:
            variantes.append((f'{pref}_{lab}_{"".join(tipos[:1]) or "L1"}_{nome_sis}.{ext}', sist, list(tipos),
                              f'Sessão de {txt_d}, {nome_sis} {"+".join(tipos)}'))
        if dupla:
            variantes.append((f'{pref}_{lab}_L1L2_{nome_sis}.{ext}', sist, sel2,
                              f'Sessão de {txt_d}, {nome_sis} L1/L2'))
        for nome, s, sel, desc in variantes:
            com = f'rinex_edit: {t_ini:%H:%M:%S}+{d:g}min {s or "todos"} {",".join(sel) if sel else "todas obs"}'
            r = _gerar(R, t_ini, t_fim, s, sel, com, antena, altura)
            if r is None:
                raise ValueError(f'nenhuma época/satélite para {nome} a partir de {t_ini:%d/%m/%Y %H:%M:%S}')
            yield nome, r[0], r[1], desc + aviso