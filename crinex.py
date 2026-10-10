# -*- coding: utf-8 -*-
"""Descompactação Hatanaka (CRINEX 1.0 e 3.0 -> RINEX 2/3) em Python puro (sem executáveis externos)."""
import math


def _diff(old, new):
    """Diferença de texto: ' ' = mantém, '&' = vira branco, outro = substitui. Resultado tem tamanho máximo."""
    n = max(len(old), len(new)); o = old.ljust(n)
    return ''.join(' ' if c == '&' else (o[i] if c == ' ' else c) for i, c in enumerate(new.ljust(n)))


def _fmt_clock(v, dec, larg):
    t = _fmt(v, dec, larg).strip()          # estilo Fortran: sem zero antes do ponto (como o crx2rnx oficial)
    t = t[1:] if t.startswith('0.') else ('-' + t[2:] if t.startswith('-0.') else t)
    return t.rjust(larg)


def _fmt(v, dec, larg):
    s = abs(v); d = 10 ** dec
    return (('-' if v < 0 else '') + '%d.%0*d' % (s // d, dec, s % d)).rjust(larg)


class _Arco:
    __slots__ = ('m', 'n', 'd')

    def __init__(self, m, v):
        self.m, self.n, self.d = m, 0, [v] + [0] * m

    def passo(self, delta):
        o = min(self.n + 1, self.m); self.d[o] = delta
        for k in range(o - 1, -1, -1):
            self.d[k] += self.d[k + 1]
        self.n += 1
        return self.d[0]


def _valor(item, arco):
    """Devolve (valor_inteiro ou None, arco)."""
    if item.strip() == '':
        return None, None
    if '&' in item:
        m, v = item.split('&'); arco = _Arco(int(m), int(v)); return arco.d[0], arco
    if arco is None:
        raise ValueError('diferença sem inicialização')
    return arco.passo(int(item)), arco


def crx2rnx(dados):
    if isinstance(dados, bytes):
        dados = dados.decode('latin-1')
    L = dados.replace('\r', '').split('\n')
    v3 = L[0][:3].strip().startswith('3')
    i = 2; cab = []
    while True:
        cab.append(L[i]); i += 1
        if 'END OF HEADER' in cab[-1]:
            break
    ntip = {}
    if v3:
        for l in cab:
            if 'SYS / # / OBS TYPES' in l and l[0] != ' ':
                ntip[l[0]] = int(l[3:6])
    else:
        n = next(int(l[:6]) for l in cab if '# / TYPES OF OBSERV' in l)
    out = list(cab)
    INIT = '>' if v3 else '&'
    prev, estado, clock = '', {}, None
    while i < len(L):
        l = L[i]
        if l == '' and i >= len(L) - 1:
            break
        if l[:1] == INIT:
            prev, estado, clock = '', {}, None
        ep = _diff(prev, l); prev = ep
        ep = ep.rstrip()
        flag = int(ep[31 if v3 else 28]); ns = int(ep[32:35] if v3 else ep[29:32])
        i += 1
        if flag > 1:                                   # evento: linhas copiadas como estão
            out.append(ep[:35 if v3 else 32].rstrip())
            for _ in range(ns):
                out.append(L[i]); i += 1
            prev, estado, clock = '', {}, None
            continue
        sats = [ep[len(ep) - 3 * ns + 3 * k: len(ep) - 3 * ns + 3 * k + 3] for k in range(ns)]
        cl = L[i]; i += 1                              # linha do relógio do receptor
        ck_val = None
        if cl.strip():
            ck_val, clock = _valor(cl.strip(), clock)
        else:
            clock = None
        cur = {}
        linhas = []
        for s in sats:
            n_t = ntip.get(s[0], 0) if v3 else n
            ant = estado.get(s) or ([None] * n_t, ' ' * (2 * n_t))
            partes = L[i].split(' ', n_t); i += 1
            itens = (partes[:n_t] + [''] * n_t)[:n_t]; fl_novo = partes[n_t] if len(partes) > n_t else ''
            base = list(ant[1].ljust(2 * n_t))                      # observação com arco reiniciado ('n&valor'): descarta os flags anteriores
            for k in range(n_t):
                if '&' in itens[k]:
                    base[2 * k] = base[2 * k + 1] = ' '
            base = ''.join(base)
            fl = _diff(base, fl_novo).ljust(2 * n_t)[:2 * n_t] if fl_novo else base
            arcos, campos = [], []
            for k in range(n_t):
                val, arco = _valor(itens[k], ant[0][k])
                arcos.append(arco)
                campos.append(' ' * 16 if val is None else _fmt(val, 3, 14) + fl[2 * k] + fl[2 * k + 1])
            cur[s] = (arcos, fl)
            if v3:
                linhas.append(s + ''.join(campos))
            else:
                linhas += [''.join(campos[j:j + 5]) for j in range(0, n_t, 5)]
        estado = cur
        if v3:
            out.append(ep[:35].rstrip() if ck_val is None else ep[:35].ljust(41) + _fmt_clock(ck_val, 12, 15))
        else:
            cab_ep = ep[:32]
            sl = [''.join(sats[j:j + 12]) for j in range(0, ns, 12)] or ['']
            ck = '' if ck_val is None else _fmt_clock(ck_val, 9, 12)
            out.append((cab_ep + sl[0]).ljust(68) + ck if ck else cab_ep + sl[0])
            out += [' ' * 32 + t for t in sl[1:]]
        out += linhas
    n_cab = len(cab)                         # como o crx2rnx oficial: linhas de dados sem espaços no fim, cabeçalho intacto
    return '\n'.join(out[:n_cab] + [l.rstrip() for l in out[n_cab:]]) + '\n'