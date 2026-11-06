# -*- coding: utf-8 -*-
"""Aba "Sobre" (dados do autor). Importado pelo app.py antes de modelos.sincronizar():
acrescenta o link no menu de todas as telas e cria a tela sobre.html."""
from modelos import TEMPLATES


# ------------------------------------------------------------------ aba "Sobre" (dados do criador)
_ICONE_SOBRE = ('<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" '
                'stroke-linejoin="round"><circle cx="12" cy="8" r="3.6"/><path d="M4.5 20.5c1.2-3.8 4.1-5.6 7.5-5.6s6.3 1.8 7.5 5.6"/></svg>')
_LINK_SOBRE = ('<a class="tab{on}" href="/sobre" title="Sobre o sistema e o autor">' + _ICONE_SOBRE +
               '<div class="tt"><b>Sobre</b></div></a>')

for _n in list(TEMPLATES):
    if 'href="/sobre"' not in TEMPLATES[_n]:
        TEMPLATES[_n] = TEMPLATES[_n].replace('</nav>', _LINK_SOBRE.format(on='') + '</nav>', 1)

_base = TEMPLATES['tempo']
_cab = _base[:_base.index('<main>')]
_cab = (_cab.replace('<title>Sistemas de tempo</title>', '<title>Sobre o sistema</title>')
            .replace('<h1>Sistemas de tempo</h1>', '<h1>Sobre o sistema</h1>')
            .replace('<span>UTC, UT1, TAI, TGPS, TT, TCG, Data Juliana e semana GPS · Aula 7</span>',
                     '<span>Autoria, instituição e módulos</span>')
            .replace('<a class="tab on" href="/tempo"', '<a class="tab" href="/tempo"')
            .replace('<a class="tab" href="/sobre"', '<a class="tab on" href="/sobre"'))

TEMPLATES['sobre'] = _cab + '''<main>
<section class="card">
<div style="display:flex;gap:20px;align-items:center;flex-wrap:wrap">
<div style="width:84px;height:84px;border-radius:50%;background:#12384f;color:#fff;display:flex;align-items:center;justify-content:center;font-size:30px;font-weight:700;flex:none">IS</div>
<div><h2 style="margin:0">Ivailton Pinheiro dos Santos</h2>
<p style="margin:4px 0 0">Graduando em <b>Engenharia Cartográfica e de Agrimensura</b></p>
<p style="margin:2px 0 0">Universidade Federal Rural da Amazônia (<b>UFRA</b>) · Campus Belém (ICIBE) · Belém, Pará</p></div></div>
</section>
<section class="card"><h2>Sobre o sistema</h2>
<p>Sistema web desenvolvido em Python (Flask) para apoiar os trabalhos práticos de <b>Geodésia</b>, reunindo em um só lugar o processamento e a análise de dados GNSS e de parâmetros de orientação da Terra.</p>
<div class="rolar"><table class="tt3"><tr><th>Aba</th><th>O que faz</th></tr>
<tr><td>Efemérides</td><td>Posição de satélites GNSS por efemérides transmitidas e comparação com as precisas (SP3)</td></tr>
<tr><td>Observáveis</td><td>Combinações Lg, Lc e Pc e detecção de saltos de ciclo em dados da RBMC</td></tr>
<tr><td>Qualidade (TEQC)</td><td>Completude, multicaminho (MP1, MP2), SNR e saltos de ciclo, no estilo do teqc +qc</td></tr>
<tr><td>Editar RINEX</td><td>Corte de sessões, filtro de constelação e de observáveis para o IBGE-PPP</td></tr>
<tr><td>PPP (IBGE)</td><td>Leitura dos resultados do IBGE-PPP: σ3D, convergência e gráficos comparativos</td></tr>
<tr><td>Movimento do polo</td><td>Coordenadas xp, yp do IERS (Boletim B) em torno de um evento</td></tr>
<tr><td>Sistemas de tempo</td><td>Conversões entre UTC, UT1, TAI, TGPS, TT, Data Juliana e semana GPS</td></tr></table></div>
<p style="margin-top:12px">Tecnologias: Python, Flask, NumPy, Matplotlib. Os resultados devem ser conferidos antes do uso em relatórios.</p>
<p style="color:#5a6675;font-size:13px">© 2026 Ivailton Pinheiro dos Santos · UFRA</p></section>
</main></body></html>
'''