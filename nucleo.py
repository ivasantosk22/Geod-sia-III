# -*- coding: utf-8 -*-
"""Núcleo: leitura de RINEX de navegação (v2 GPS / v3 GPS+Galileo), SP3, cálculo da
posição (mesmas fórmulas do Orbita_galileo_Final.sce) e gráficos."""

import io
import os
import zipfile
import subprocess
from pathlib import Path
from dataclasses import dataclass
from datetime import datetime, timedelta
import numpy as np
from matplotlib.figure import Figure

GM, WE, PI = 3.986004418e14, 7.2921151467e-5, 3.1415926535898
GPS0 = datetime(1980, 1, 6)


def _n(t):
    t = t.strip()
    return float(t.replace('D', 'E').replace('d', 'e')) if t else 0.0


@dataclass
class Efemeride:
    sat: str
    toc: datetime
    clk: tuple      # a0, a1, a2
    m: np.ndarray   # matriz 7x4 (igual à matriz "inp" do Scilab)
    raw: list

    @property
    def fonte(self):
        if self.sat[0] == 'E':
            f = int(self.m[4, 1])
            nomes = [n for b, n in ((1, 'I/NAV E1-B'), (2, 'F/NAV E5a-I'), (4, 'I/NAV E5b-I')) if f & b]
            return f"IOD {int(self.m[0, 0])} · fonte {f} ({', '.join(nomes)})"
        return f"IODE {int(self.m[0, 0])}"


def ler_nav(texto):
    L = texto.splitlines()
    v3 = float(L[0][:9]) >= 3
    ini = next(i for i, l in enumerate(L) if 'END OF HEADER' in l) + 1
    corpo = [l for l in L[ini:] if l.strip()]
    if v3:
        blocos = []
        for l in corpo:
            if l[0] != ' ':
                blocos.append([l])
            elif blocos:
                blocos[-1].append(l)
    else:
        blocos = [corpo[i:i + 8] for i in range(0, len(corpo), 8)]
    out = []
    for b in blocos:
        p = b[0].rstrip()
        head, tail = p[:-57], p[-57:]            # 3 coeficientes do relógio = últimos 57 caracteres
        off = 4 if v3 else 3
        if v3:
            if p[0] not in 'GE':
                continue
            sat, f = head[:3].replace(' ', '0'), head[3:].split()
        else:
            t = head.split()
            sat, f = 'G%02d' % int(t[0]), t[1:7]
            f[0] = str(int(f[0]) + (2000 if int(f[0]) < 80 else 1900))
        toc = datetime(*[int(float(x)) for x in f])
        clk = tuple(_n(tail[19 * i:19 * (i + 1)]) for i in range(3))
        m = np.zeros((7, 4))
        for r, l in enumerate(b[1:8]):
            for c in range(4):
                m[r, c] = _n(l[off + 19 * c:off + 19 * (c + 1)])
        out.append(Efemeride(sat, toc, clk, m, [x.rstrip() for x in b]))
    return sorted(out, key=lambda e: (e.sat, e.toc))


def ler_sp3(texto):
    d, t = {}, None
    for l in texto.splitlines():
        if l.startswith('*'):
            t = datetime(*[int(float(x)) for x in l[1:].split()[:6]])
        elif l.startswith('P') and t:
            s = l[1:4].strip()
            s = 'G%02d' % int(s) if s.isdigit() else s.replace(' ', '0')
            x, y, z = (float(l[a:a + 14]) for a in (4, 18, 32))
            if x != 0 and x < 999999:
                d.setdefault(s, {})[t] = (x, y, z)
    return d


def calcular(e, ini, passo_min, n, gm_icd=False, aplicar_dts=True):
    """Mesmo cálculo do Scilab. Retorna (instantes, matriz n x 3 em km, avisos)."""
    m, (a0, a1, a2) = e.m, e.clk
    Crs, dn, M0 = m[0, 1:4]
    Cuc, ec, Cus, ra = m[1]
    Toe, Cic, O0, Cis = m[2]
    i0, Crc, w, Od = m[3]
    idot, trs = m[4, 0], m[6, 0]
    a = ra ** 2
    mu = 3.986005e14 if (gm_icd and e.sat[0] == 'G') else GM   # padrão da aula: GM WGS-84
    toc_rec = (e.toc - GPS0).total_seconds() % 604800          # toc do registro (slides 7-8)
    tempos, P, dmax = [], [], 0
    for k in range(n):
        t = ini + timedelta(minutes=k * passo_min)
        Toc = (t - GPS0).total_seconds() % 604800          # segundos da semana GPS
        x_ = (trs - toc_rec + 302400) % 604800 - 302400
        dts = a0 + a1 * x_ + a2 * x_ ** 2                        # dts = a0 + a1(trs-toc) + a2(trs-toc)^2
        Tg = Toc - dts if aplicar_dts else Toc
        dt = (Tg - Toe + 302400) % 604800 - 302400
        dmax = max(dmax, abs(dt))
        Mk = M0 + (np.sqrt(mu / a ** 3) + dn) * dt
        Ek = Mk
        for _ in range(8):
            Ek = Mk + ec * np.sin(Ek)
        V1 = (np.cos(Ek) - ec) / (1 - ec * np.cos(Ek))
        V2 = np.sqrt(1 - ec ** 2) * np.sin(Ek) / (1 - ec * np.cos(Ek))
        Vk = np.arctan2(V2, V1) % (2 * PI)                  # análise de quadrante do Scilab
        fi = Vk + w
        u = fi + Cuc * np.cos(2 * fi) + Cus * np.sin(2 * fi)
        r = a * (1 - ec * np.cos(Ek)) + Crc * np.cos(2 * fi) + Crs * np.sin(2 * fi)
        ik = i0 + idot * dt + Cic * np.cos(2 * fi) + Cis * np.sin(2 * fi)
        x, y = r * np.cos(u), r * np.sin(u)
        Ok = O0 + Od * dt - WE * Tg
        P.append([(x * np.cos(Ok) - y * np.sin(Ok) * np.cos(ik)) / 1000,
                  (x * np.sin(Ok) + y * np.cos(Ok) * np.cos(ik)) / 1000,
                  (y * np.sin(ik)) / 1000])
        tempos.append(t)
    av = []
    if dmax > 14400:
        av.append(f'Há épocas a {dmax / 3600:.1f} h do Toe: fora da validade típica (~4 h) da efeméride.')
    return tempos, np.array(P), av


def parametros(e):
    m, (a0, a1, a2) = e.m, e.clk
    L = [('Época do relógio (Toc)', e.toc.strftime('%d/%m/%Y %H:%M:%S'), ''),
         ('a0', a0, 's'), ('a1', a1, 's/s'), ('a2', a2, 's/s²'),
         ('IOD', m[0, 0], ''), ('Crs', m[0, 1], 'm'), ('Δn', m[0, 2], 'rad/s'), ('M0', m[0, 3], 'rad'),
         ('Cuc', m[1, 0], 'rad'), ('e', m[1, 1], ''), ('Cus', m[1, 2], 'rad'), ('√a', m[1, 3], 'm^½'),
         ('Toe', m[2, 0], 's'), ('Cic', m[2, 1], 'rad'), ('Ω0', m[2, 2], 'rad'), ('Cis', m[2, 3], 'rad'),
         ('i0', m[3, 0], 'rad'), ('Crc', m[3, 1], 'm'), ('ω', m[3, 2], 'rad'), ('Ω̇', m[3, 3], 'rad/s'),
         ('IDOT', m[4, 0], 'rad/s'), ('Semana (registro)', m[4, 2], ''), ('Tempo de transmissão', m[6, 0], 's')]
    return [(a, v if isinstance(v, str) else f'{v:.12g}', u) for a, v, u in L]


def _png(fig):
    b = io.BytesIO()
    fig.savefig(b, format='png', dpi=150)
    return b.getvalue()


def graficos(rot, P, Ps, terra, titulo):
    """Retorna [png1, png2 ou None, png3]."""
    cor = ['#1f77b4', '#ff7f0e', '#d62728']
    xs = np.arange(len(rot))
    f1 = Figure(figsize=(9, 5)); ax = f1.subplots()
    for j, c in enumerate('XYZ'):
        ax.plot(xs, P[:, j], '-o', color=cor[j], label=c)
    ax.set(xlabel='Hora (tempo GPS)', ylabel='Coordenadas (km)', title='Evolução de X(t), Y(t), Z(t)\n' + titulo)
    ax.set_xticks(xs); ax.set_xticklabels(rot, rotation=45); ax.grid(alpha=.3); ax.legend(); f1.tight_layout()

    g2 = None
    if Ps is not None and not np.isnan(Ps).all():
        d = (P - Ps) * 1000
        f2 = Figure(figsize=(9, 5)); ax = f2.subplots()
        for j, c in enumerate('XYZ'):
            ax.bar(xs + (j - 1) * .27, d[:, j], .27, color=cor[j], label='Δ' + c)
        ax.axhline(0, color='k', lw=.8)
        ax.set(xlabel='Hora (tempo GPS)', ylabel='Discrepância (m)', title='Discrepâncias (transmitida − SP3)\n' + titulo)
        ax.set_xticks(xs); ax.set_xticklabels(rot, rotation=45); ax.grid(axis='y', alpha=.3); ax.legend(); f2.tight_layout()
        g2 = _png(f2)

    f3 = Figure(figsize=(8, 7)); ax = f3.add_subplot(111, projection='3d')
    ax.plot(P[:, 0], P[:, 1], P[:, 2], '-o', color='crimson', lw=2.5, ms=6, label='Trajetória (transmitidas)')
    ax.plot(*P[0], 's', color='green', ms=10, label=f'Início ({rot[0]})')
    ax.plot(*P[-1], '^', color='black', ms=10, label=f'Fim ({rot[-1]})')
    for k in range(len(rot)):
        ax.text(*P[k], '  ' + rot[k], fontsize=7)
    if terra:
        u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 40), np.linspace(0, np.pi, 20))
        R = 6378.137
        ax.plot_surface(R * np.cos(u) * np.sin(v), R * np.sin(u) * np.sin(v), R * np.cos(v), color='lightsteelblue', alpha=.35, linewidth=0)
        ax.set_xlim(-26000, 26000); ax.set_ylim(-26000, 26000); ax.set_zlim(-26000, 26000)
    else:
        c = (P.max(0) + P.min(0)) / 2
        h = max(np.ptp(P, axis=0).max(), 1) * .575
        ax.set_xlim(c[0] - h, c[0] + h); ax.set_ylim(c[1] - h, c[1] + h); ax.set_zlim(c[2] - h, c[2] + h)
    ax.set_box_aspect((1, 1, 1)); ax.view_init(elev=20, azim=-60)
    ax.set_xlabel('X (km)', labelpad=10); ax.set_ylabel('Y (km)', labelpad=10); ax.set_zlabel('Z (km)', labelpad=10)
    ax.set_title('Órbita reconstruída em 3D\n' + titulo); ax.legend(loc='upper left', fontsize=8); f3.tight_layout()
    return [_png(f1), g2, _png(f3)]


# --- NOVA FUNÇÃO ADICIONADA PARA PROCESSAMENTO GNSS (CRX2RNX e TEQC) ---
def processar_dados_gnss(caminho_zip, pasta_saida=None):
    caminho_zip = Path(caminho_zip)
    if pasta_saida is None:
        pasta_saida = caminho_zip.parent
    else:
        pasta_saida = Path(pasta_saida)
        pasta_saida.mkdir(parents=True, exist_ok=True)

    # 1. Extração do ficheiro ZIP
    with zipfile.ZipFile(caminho_zip, 'r') as zip_ref:
        zip_ref.extractall(pasta_saida)
        arquivos_extraidos = zip_ref.namelist()

    # Identifica o ficheiro de observação Hatanaka (.25d, .24d, etc.)
    arquivo_hatanaka = next((pasta_saida / f for f in arquivos_extraidos if f.endswith('d')), None)

    if not arquivo_hatanaka:
        return "Erro: Nenhum ficheiro Hatanaka (.d) encontrado no zip."

    # 2. Conversão Hatanaka para RINEX (crx2rnx)
    try:
        subprocess.run(["crx2rnx", str(arquivo_hatanaka)], check=True, cwd=pasta_saida)
    except subprocess.CalledProcessError as e:
        return f"Erro na conversão crx2rnx: {e}"
    except FileNotFoundError:
        return "Erro: Executável 'crx2rnx' não encontrado no PATH do sistema."

    arquivo_rinex = str(arquivo_hatanaka)[:-1] + 'o'

    if not os.path.exists(arquivo_rinex):
        return "Erro: O ficheiro RINEX não foi gerado."

    # 3. Execução do TEQC
    try:
        subprocess.run(
            ["teqc", "+qc", Path(arquivo_rinex).name],
            cwd=pasta_saida,
            capture_output=True,
            text=True,
            check=True
        )
        return f"Processamento concluído com sucesso. Relatórios gerados em {pasta_saida}"
    
    except subprocess.CalledProcessError as e:
        return f"Erro ao executar TEQC: {e.stderr}"
    except FileNotFoundError:
        return "Erro: Executável 'teqc' não encontrado no PATH do sistema."