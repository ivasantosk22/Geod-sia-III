# -*- coding: utf-8 -*-
"""Sistema local (Flask): efemérides transmitidas x precisas (SP3).
Execute:  python app.py   -> abre http://127.0.0.1:5000"""
import base64, csv, gzip, io, os, re, threading, uuid, webbrowser, zipfile
from datetime import datetime
import numpy as np
from flask import Flask, render_template, request, redirect, url_for, send_file, abort, flash
import nucleo, observaveis
import controle_qualidade as teqc
import tempo as TEMPO, polo as POLO_M, rinex_edit as RE, ppp as PPP
from datetime import date as _date

BASE = os.path.dirname(os.path.abspath(__file__))
import modelos
import sobre                                   # acrescenta a aba "Sobre" (dados do autor) às telas
modelos.sincronizar(BASE)                      # recria telas ausentes, vazias ou antigas
app = Flask(__name__)
app.config['TEMPLATES_AUTO_RELOAD'] = True      # usa o arquivo .html novo sem precisar reiniciar o servidor
app.jinja_env.auto_reload = True


@app.after_request
def sem_cache(r):                               # o navegador não guarda páginas (evita tela em branco antiga)
    if r.mimetype == 'text/html':
        r.headers['Cache-Control'] = 'no-store'
    return r
app.secret_key = 'local'
app.config['MAX_CONTENT_LENGTH'] = 200 * 1024 * 1024
SESSOES = {}
CAB = ['Época', 'X (km)', 'Y (km)', 'Z (km)', 'X SP3 (km)', 'Y SP3 (km)', 'Z SP3 (km)', 'ΔX (m)', 'ΔY (m)', 'ΔZ (m)', '3D (m)']


def _achar(dados, pad):
    """Procura em .zip/.gz (inclusive aninhados) o arquivo cujo nome casa com pad."""
    if dados[:2] == b'\x1f\x8b':
        dados = gzip.decompress(dados)
    if dados[:2] != b'PK':
        return dados
    z = zipfile.ZipFile(io.BytesIO(dados))
    for n in z.namelist():
        if re.search(pad, n, re.I) or n.lower().endswith(('.zip', '.gz')):
            r = _achar(z.read(n), pad)
            if r is not None and not (r[:2] == b'PK'):
                return r
    return None


def _texto(dados, campo):
    """Aceita texto puro, .gz e .zip (escolhe o .sp3 ou o arquivo de navegação dentro do zip)."""
    if dados[:2] == b'\x1f\x8b':
        dados = gzip.decompress(dados)
    elif dados[:2] == b'PK':
        pad = r'\.sp3$' if campo == 'sp3' else r'\.(\d\d[nlp]|n|rnx)$'
        dados = _achar(dados, pad)
        if dados is None:
            raise ValueError('o .zip não contém ' + ('um .sp3' if campo == 'sp3' else 'arquivo de navegação (.26n)'))
    elif dados[:2] == b'\x1f\x9d':
        raise ValueError('arquivo .Z não é suportado: descompacte com o 7-Zip e envie o .sp3')
    if b'\x00' in dados[:2000]:
        raise ValueError('o arquivo não parece ser de texto (está compactado?)')
    return dados.decode('latin-1')


def _ler(campo, exemplo, nome_ex):
    if exemplo:
        return open(os.path.join(BASE, 'exemplos', nome_ex), encoding='latin-1').read()
    arq = request.files.get(campo)
    return _texto(arq.read(), campo) if arq and arq.filename else None


@app.route('/')
def index():
    return render_template('pagina.html', sid=None)


@app.route('/enviar', methods=['POST'])
def enviar():
    ex = 'exemplo' in request.form
    try:
        nav = _ler('nav', ex, 'transmitida.14l')
        sp3 = _ler('sp3', ex, 'tum18022.sp3')
        if not nav:
            raise ValueError('envie o arquivo de navegação')
        efs = nucleo.ler_nav(nav)
        if not efs:
            raise ValueError('nenhuma efeméride GPS/Galileo encontrada')
        sid = uuid.uuid4().hex[:8]
        SESSOES[sid] = {'efs': efs, 'sp3': nucleo.ler_sp3(sp3) if sp3 else {}, 'res': None}
        return redirect(url_for('painel', sid=sid))
    except Exception as ex_:
        flash(f'Não foi possível ler os arquivos: {ex_}')
        return redirect(url_for('index'))


def executar(S, e, ini, passo, n, terra, gm_icd=False, aplicar_dts=True):
    tempos, P, avisos = nucleo.calcular(e, ini, passo, n, gm_icd, aplicar_dts)
    dias = len({t.date() for t in tempos}) > 1
    rot = [t.strftime('%d/%m %H:%M' if dias else '%H:%M') for t in tempos]
    sp = S['sp3'].get(e.sat)
    Ps = np.array([sp.get(t, (np.nan,) * 3) for t in tempos]) if sp else None
    if Ps is None:
        avisos.append('Sem SP3 para este satélite: comparação não realizada.')
    elif np.isnan(Ps).any():
        ts = sorted(sp)
        pas = int((ts[1] - ts[0]).total_seconds() // 60) if len(ts) > 1 else 0
        avisos.append(f'{int(np.isnan(Ps[:, 0]).sum())} época(s) sem correspondência no SP3. O SP3 deste satélite cobre de '
                      f'{ts[0]:%d/%m/%Y %H:%M} a {ts[-1]:%d/%m/%Y %H:%M} (passo de {pas} min): confira se é o dia certo e se o passo é múltiplo dele.')
    D = (P - Ps) * 1000 if Ps is not None else np.full_like(P, np.nan)
    D3 = np.linalg.norm(D, axis=1)
    f = lambda v, d=6: '—' if np.isnan(v) else f'{v:.{d}f}'
    rows = []
    for k in range(len(tempos)):
        s = Ps[k] if Ps is not None else (np.nan,) * 3
        rows.append([rot[k], *[f(v) for v in P[k]], *[f(v) for v in s], *[f(v, 3) for v in D[k]], f(D3[k], 3)])
    ok = ~np.isnan(D3)
    est = [(nm, f(np.nanmean(c), 3), f(np.sqrt(np.nanmean(c ** 2)), 3), f(np.nanmax(np.abs(c)), 3))
           for nm, c in (('ΔX', D[:, 0]), ('ΔY', D[:, 1]), ('ΔZ', D[:, 2]), ('3D', D3))] if ok.any() else []
    resumo = ''
    if ok.sum() > 2:
        th = np.array([(t - tempos[0]).total_seconds() / 3600 for t in tempos])
        inc = np.polyfit(th[ok], D3[ok], 1)[0]
        resumo = (f'Discrepância 3D com RMS de {np.sqrt(np.mean(D3[ok] ** 2)):.2f} m; tendência '
                  f'{"crescente" if inc > 0 else "decrescente"} de {abs(inc):.3f} m/h no intervalo.')
    dd = (ini - nucleo.GPS0).days
    titulo = f'{e.sat} - {ini:%d/%m/%Y} (semana GPS {dd // 7}, dia {dd % 7})'
    pngs = nucleo.graficos(rot, P, Ps, terra, titulo)
    return dict(sat=e.sat, semana=dd // 7, dow=dd % 7, doy=ini.timetuple().tm_yday, ini=ini, passo=passo, n=n,
                avisos=avisos, gm='3,986005e14' if (gm_icd and e.sat[0] == 'G') else '3,986004418e14', dts=aplicar_dts, rows=rows, est=est, resumo=resumo, pngs=pngs, params=nucleo.parametros(e),
                raw='\n'.join(e.raw), b64=[base64.b64encode(p).decode() if p else None for p in pngs])


@app.route('/painel/<sid>', methods=['GET', 'POST'])
def painel(sid):
    S = SESSOES.get(sid) or abort(404)
    if request.method == 'POST':
        try:
            fm = request.form
            S['res'] = executar(S, S['efs'][int(fm['reg'])], datetime.fromisoformat(fm['inicio']),
                                float(fm['passo']), min(int(fm['n']), 500), 'terra' in fm, 'gm_icd' in fm, 'sem_dts' not in fm)
        except Exception as ex:
            flash(f'Erro no cálculo: {ex}')
    regs = [dict(i=i, sat=e.sat, toc=e.toc.strftime('%Y-%m-%dT%H:%M'),
                 rot=f'{e.toc:%d/%m/%Y %H:%M} · {e.fonte}') for i, e in enumerate(S['efs'])]
    return render_template('pagina.html', sid=sid, regs=regs, sats=sorted({r['sat'] for r in regs}),
                           tem_sp3=bool(S['sp3']), res=S['res'], cab=CAB, f=request.form)


@app.route('/baixar/<sid>/<what>')
def baixar(sid, what):
    res = (SESSOES.get(sid) or {}).get('res') or abort(404)
    if what == 'csv':
        b = io.StringIO(); w = csv.writer(b); w.writerow(CAB); w.writerows(res['rows'])
        return send_file(io.BytesIO(b.getvalue().encode('utf-8-sig')), mimetype='text/csv',
                         as_attachment=True, download_name=f"resultado_{res['sat']}.csv")
    k = {'g1': 0, 'g2': 1, 'g3': 2}.get(what)
    if k is None or not res['pngs'][k]:
        abort(404)
    return send_file(io.BytesIO(res['pngs'][k]), mimetype='image/png', as_attachment=True,
                     download_name=f"grafico{k + 1}_{res['sat']}.png")


OBS = {}


@app.route('/obs')
def obs_index():
    return render_template('obs.html', sid=None)


@app.route('/obs/enviar', methods=['POST'])
def obs_enviar():
    try:
        arqs = []
        for x in request.files.getlist('arquivos'):
            b = x.read() if x.filename else b''
            if x.filename and not (b[:2] not in (b'PK', b'\x1f\x8b') and b'NAV' in b[:600] and b'OBSERVATION DATA' not in b[:600]):   # ignora navegação
                arqs.append((b, x.filename))
        if not arqs:
            raise ValueError('envie ao menos um arquivo de observação')
        fs = [f for b, n_ in arqs for f in observaveis.carregar_todos(b, n_)]
    except Exception as ex:
        flash(f'Não foi possível ler os arquivos: {ex}')
        return redirect(url_for('obs_index'))
    sid = uuid.uuid4().hex[:8]
    OBS[sid] = {'fs': fs}
    return redirect(url_for('obs_painel', sid=sid))


@app.route('/obs/<sid>', methods=['GET', 'POST'])
def obs_painel(sid):
    S = OBS.get(sid) or abort(404)
    fm, fs = request.form, S['fs']
    try:
        limiar = float(fm.get('limiar', 0.15)); idx = min(int(fm.get('arq', 0)), len(fs) - 1)
        f = fs[idx]; ns = np.isfinite(f['A'][:, :, 0]).sum(0)
        sat = fm.get('sat') if fm.get('sat') in f['sats'] else f['sats'][int(np.argmax(ns))]
        figs = observaveis.graficos(fs, idx, sat, limiar)
        S['est'] = dict(limiar=limiar, idx=idx)
        est = sorted({x['estacao'] for x in fs}); comp = None
        if len(est) >= 2:
            ea = fm.get('estA') if fm.get('estA') in est else est[0]
            eb = fm.get('estB') if fm.get('estB') in est and fm.get('estB') != ea else next(e for e in est if e != ea)
            cab_c, rows_c, figs_c = observaveis.comparar(fs, ea, eb, limiar)
            comp = dict(cab=cab_c, rows=rows_c, ea=ea, eb=eb, figs=[(t, base64.b64encode(p).decode()) for t, p in figs_c],
                        pngs=[p for _, p in figs_c])
        res = dict(resumo=[observaveis.resumo(x, limiar) for x in fs],  figs=[(t, base64.b64encode(p).decode()) for t, p in figs],
                   pngs=[p for _, p in figs], sat=sat, limiar=limiar, idx=idx, obs=f['obs'], est=est, comp=comp)
    except Exception as ex:
        flash(f'Erro no processamento: {ex}')
        return redirect(url_for('obs_index'))
    S['res'] = res
    return render_template('obs.html', sid=sid, res=res, cab=observaveis.RES_CAB, arqs=[(i, x['nome'], x['estacao']) for i, x in enumerate(fs)],
                           sats=f['sats'])


@app.route('/obs/<sid>/baixar/<what>')
def obs_baixar(sid, what):
    S = OBS.get(sid) or abort(404); res = S.get('res') or abort(404)
    if what in ('resumo', 'dados'):
        txt = (','.join(observaveis.RES_CAB) + '\n' + '\n'.join(','.join(map(str, r)) for r in res['resumo'])) if what == 'resumo' \
            else observaveis.csv_dados(S['fs'][res['idx']], res['limiar'])
        return send_file(io.BytesIO(txt.encode('utf-8-sig')), mimetype='text/csv', as_attachment=True, download_name=f'obs_{what}.csv')
    if what == 'comp' and res.get('comp'):
        c = res['comp']; txt = ','.join(c['cab']) + '\n' + '\n'.join(','.join(map(str, r)) for r in c['rows'])
        return send_file(io.BytesIO(txt.encode('utf-8-sig')), mimetype='text/csv', as_attachment=True, download_name='obs_comparacao.csv')
    if what[:1] == 'c' and what[1:].isdigit() and res.get('comp'):
        k = int(what[1:]) - 1
        if k >= len(res['comp']['pngs']):
            abort(404)
        return send_file(io.BytesIO(res['comp']['pngs'][k]), mimetype='image/png', as_attachment=True, download_name=f'obs_comparacao{k + 1}.png')
    k = int(what[1:]) - 1 if what[:1] == 'g' and what[1:].isdigit() else None
    if k is None or k >= len(res['pngs']):
        abort(404)
    return send_file(io.BytesIO(res['pngs'][k]), mimetype='image/png', as_attachment=True, download_name=f'obs_grafico{k + 1}.png')


QC = {}


@app.route('/qc')
def qc_index():
    return render_template('qc.html', sid=None)


@app.route('/qc/enviar', methods=['POST'])
def qc_enviar():
    try:
        obs_, navs = [], []
        for x in request.files.getlist('arquivos'):
            if not x.filename:
                continue
            b = x.read()
            if b[:2] not in (b'PK', b'\x1f\x8b') and b'NAV' in b[:600] and b'OBSERVATION DATA' not in b[:600]:
                navs.append(_texto(b, 'nav'))                              # navegação enviada junto com a observação
            else:
                obs_.append((b, x.filename))
        nv = request.files.get('nav')
        if nv and nv.filename:
            navs.append(_texto(nv.read(), 'nav'))
        if not obs_:
            raise ValueError('envie ao menos um arquivo de observação (.o, .d, .zip ou .gz)')
        fs = [f for b, n_ in obs_ for f in teqc.carregar_todos(b, n_, navs)]
    except Exception as ex:
        flash(f'Não foi possível ler os arquivos: {ex}')
        return redirect(url_for('qc_index'))
    sid = uuid.uuid4().hex[:8]
    QC[sid] = {'fs': fs}
    return redirect(url_for('qc_painel', sid=sid))


@app.route('/qc/<sid>', methods=['GET', 'POST'])
def qc_painel(sid):
    S = QC.get(sid) or abort(404)
    fm, fs = request.form, S['fs']
    try:
        par = dict(mask=float(fm.get('mask', 10)), iod=float(fm.get('iod', 400)), mp=float(fm.get('mp', 5)), arco=int(fm.get('arco', 20)))
        idx = min(int(fm.get('arq', 0)), len(fs) - 1); f = fs[idx]
        Rs = [teqc.analisar(x, par['mask'], par['iod'], par['mp'], par['arco']) for x in fs]; R = Rs[idx]
        ov = S.setdefault('over', {})
        if request.method == 'POST' and 'rec_nome' in fm:                         # nomes digitados pelo usuário (valem para o arquivo escolhido)
            ov[idx] = dict(rec=fm.get('rec_nome', '').strip(), ant=fm.get('ant_nome', '').strip())

        def ident(i, ff):
            h, o_ = ff.get('hdr') or {}, ov.get(i, {})
            vazio = lambda t: (t or '').strip().lower() in ('', 'unknown', 'unknown ext', 'unknown ext.', 'none')
            rec, ant = o_.get('rec') or h.get('rec') or '', o_.get('ant') or h.get('ant') or ''
            d = h.get('delta')
            return dict(rec=rec or '—', ant=ant or '—', rec_num=h.get('rec_num') or '—', rec_ver=h.get('rec_ver') or '—', ant_num=h.get('ant_num') or '—',
                        observador=' / '.join(x for x in (h.get('observador'), h.get('agencia')) if x) or '—', pgm=h.get('pgm') or '—',
                        altura=f'{d[0]:.4f} m' if d else '—', vazio_rec=vazio(rec), vazio_ant=vazio(ant),
                        dig_rec=o_.get('rec', ''), dig_ant=o_.get('ant', ''))
        ns = np.isfinite(f['A'][:, :, 0]).sum(0)
        sat = fm.get('sat') if fm.get('sat') in f['sats'] else f['sats'][int(np.argmax(ns))]
        figs = teqc.graficos(f, R, sat, par['mask'])
        res = dict(resumo=[teqc.linha_res(r) for r in Rs], sats=teqc.linhas_sat(R), figs=[(t, base64.b64encode(p).decode()) for t, p in figs],
                   pngs=[p for _, p in figs], sat=sat, par=par, idx=idx, tem_el=R['tem_el'], nome=f"{f['estacao']} · {f['nome']}",
                   info=ident(idx, f), rec_ant=[(ident(i, x)['rec'], ident(i, x)['ant']) for i, x in enumerate(fs)])
    except Exception as ex:
        flash(f'Erro no processamento: {ex}')
        return redirect(url_for('qc_index'))
    S['res'] = res
    return render_template('qc.html', sid=sid, res=res, cab_res=teqc.CAB_RES, cab_sat=teqc.CAB_SAT,
                           arqs=[(i, x['nome'], x['estacao']) for i, x in enumerate(fs)], sats=f['sats'])


@app.route('/qc/<sid>/baixar/<what>')
def qc_baixar(sid, what):
    S = QC.get(sid) or abort(404); res = S.get('res') or abort(404)
    if what in ('resumo', 'sats'):
        cab, rows = ((teqc.CAB_RES + ['Receptor', 'Antena'], [r + list(ra) for r, ra in zip(res['resumo'], res['rec_ant'])]) if what == 'resumo' else (teqc.CAB_SAT, res['sats']))
        txt = ','.join(cab) + '\n' + '\n'.join(','.join(map(str, r)) for r in rows)
        return send_file(io.BytesIO(txt.encode('utf-8-sig')), mimetype='text/csv', as_attachment=True, download_name=f'qc_{what}.csv')
    k = int(what[1:]) - 1 if what[:1] == 'g' and what[1:].isdigit() else None
    if k is None or k >= len(res['pngs']):
        abort(404)
    return send_file(io.BytesIO(res['pngs'][k]), mimetype='image/png', as_attachment=True, download_name=f'qc_grafico{k + 1}.png')


# ------------------------------------------------------------------ Sistemas de tempo (Aula 7)
@app.route('/tempo', methods=['GET', 'POST'])
def tempo_pagina():
    ex = request.args.get('ex')
    if ex:
        e = TEMPO.EXEMPLO; v = dict(hl=e['hl'].strftime('%Y-%m-%dT%H:%M:%S'), fuso='3', dut1='0.4', tai='35')
    elif request.method == 'POST':
        v = dict(hl=request.form.get('hl', ''), fuso=request.form.get('fuso', '3'), dut1=request.form.get('dut1', '0'), tai=request.form.get('tai', ''))
    else:
        v = dict(hl=datetime.now().strftime('%Y-%m-%dT%H:%M:%S'), fuso='3', dut1='0', tai='')
    res = erro = None
    if ex or request.method == 'POST':
        try:
            hl = datetime.fromisoformat(v['hl']); fuso = float(v['fuso'].replace(',', '.')); info = []
            dut1 = float(v['dut1'].replace(',', '.') or 0); tai = int(v['tai']) if v['tai'].strip() else None
            bol = request.files.get('bol')
            if bol and bol.filename and not ex:
                dados = POLO_M.juntar([POLO_M.ler_boletim(bol.read().decode('latin-1'))])
                mjd = TEMPO.calcular(hl, fuso, 0, tai or 35)['mjd']
                u = POLO_M.ut1_utc_ms(dados, mjd)
                if u is None:
                    info.append('A data não está dentro do período do Boletim B enviado: DUT1 e TAI−UTC não foram alterados.')
                else:
                    dut1 = u / 1000; info.append(f'DUT1 = UT1−UTC = {dut1:+.6f} s, interpolado do Boletim B (valor das 0h UTC dos dias vizinhos).')
                    tb = POLO_M.tai_utc_boletim(dados, int(mjd))
                    if tb and tai is None:
                        tai = tb; info.append(f'TAI−UTC = {tb} s, da seção 4 do Boletim B.')
            r = TEMPO.calcular(hl, fuso, dut1, tai)
            if tai is None and r['utc'].date() > _date(2026, 12, 31):
                info.append('TAI−UTC = 37 s é o valor desde 2017. O Boletim C 72 (jul/2026) confirma que não haverá segundo intercalado até dez/2026; para datas posteriores, confirme no Boletim C.')
            tab = TEMPO.tabela(r); check = None
            if ex:
                m = {a_: c for a_, b_, c in tab}
                check = [(k, w, m[k], m[k] == w) for k, w in TEMPO.ESPERADO.items()]
            res = dict(tab=tab, info=info, check=check)
        except Exception as ex_:
            erro = f'Não foi possível calcular: {ex_}'
    return render_template('tempo.html', v=v, res=res, erro=erro, exemplo=bool(ex))


# ------------------------------------------------------------------ Movimento do polo (Aula 4)
POLO = {}


def _textos_boletins(arq):
    nome, b = arq.filename, arq.read()
    if b[:2] == b'PK':
        z = zipfile.ZipFile(io.BytesIO(b))
        return [(n, z.read(n).decode('latin-1')) for n in z.namelist() if n.lower().endswith(('.txt', '.dat', '.csv'))]
    return [(nome, b.decode('latin-1'))]


@app.route('/polo')
def polo_index():
    return render_template('polo.html', sid=None, eventos=POLO_M.EVENTOS)


@app.route('/polo/enviar', methods=['POST'])
def polo_enviar():
    try:
        fm = request.form; evento = _date.fromisoformat(fm['evento']); meses = int(fm.get('meses', 6))
        textos = [t for x in request.files.getlist('arquivos') if x.filename for t in _textos_boletins(x)]
        if not textos and 'baixar' in fm:
            nums = POLO_M.boletins_necessarios(POLO_M.soma_meses(evento, -meses), POLO_M.soma_meses(evento, meses))
            ok_, erros = POLO_M.baixar_boletins(nums)
            if erros:
                raise ValueError('não consegui baixar os boletins do IERS (' + erros[0] + '). Baixe-os manualmente no IERS Data Center e envie os .txt aqui')
            textos = [(f'bulletinb-{n}.txt', t) for n, t in ok_.items()]
        if not textos:
            raise ValueError('envie os boletins (.txt) ou marque a opção de baixar do IERS')
        bs = [POLO_M.ler_boletim(t) for _, t in textos]
        if not any(b['linhas'] for b in bs):
            raise ValueError('não encontrei valores de xp e yp nos arquivos (formato de Boletim B ou série C04)')
    except Exception as ex_:
        flash(f'Não foi possível ler os boletins: {ex_}')
        return redirect(url_for('polo_index'))
    sid = uuid.uuid4().hex[:8]
    POLO[sid] = dict(textos=textos, dados=POLO_M.juntar(bs), formatos=sorted({b['formato'] for b in bs}), evento=evento.isoformat(), meses=meses)
    return redirect(url_for('polo_painel', sid=sid))


@app.route('/polo/<sid>', methods=['GET', 'POST'])
def polo_painel(sid):
    S = POLO.get(sid) or abort(404); fm = request.form; D = S['dados']
    par = dict(evento=fm.get('evento', S['evento']), meses=int(fm.get('meses', S['meses'])), ma=int(fm.get('ma', 7)))
    res = erro = None
    try:
        ev = _date.fromisoformat(par['evento']); R = POLO_M.analisar(D, ev, par['meses'], par['ma']); figs = POLO_M.graficos(R)
        res = dict(R=R, est=R['est'], ind=R['ind'], figs=[(t, base64.b64encode(p).decode()) for t, p in figs], pngs=[p for _, p in figs])
        S['res'] = res; S['par'] = par
    except Exception as ex_:
        erro = str(ex_)
    ev0 = _date.fromisoformat(par['evento'])
    info = dict(n=len(D['mjd']), ini=POLO_M.mjd_data(D['mjd'][0]) if len(D['mjd']) else None, fim=POLO_M.mjd_data(D['mjd'][-1]) if len(D['mjd']) else None,
                arquivos=len(S['textos']), formatos=', '.join(S['formatos']),
                nec=POLO_M.boletins_necessarios(POLO_M.soma_meses(ev0, -par['meses']), POLO_M.soma_meses(ev0, par['meses'])))
    return render_template('polo.html', sid=sid, res=res, erro=erro, par=par, info=info, cab_est=POLO_M.CAB_EST, eventos=POLO_M.EVENTOS)


@app.route('/polo/<sid>/baixar/<what>')
def polo_baixar(sid, what):
    S = POLO.get(sid) or abort(404); res = S.get('res') or abort(404)
    if what == 'csv':
        return send_file(io.BytesIO(POLO_M.csv_polo(res['R']).encode('utf-8-sig')), mimetype='text/csv', as_attachment=True, download_name='movimento_polo.csv')
    if what == 'boletins':
        b = io.BytesIO()
        with zipfile.ZipFile(b, 'w', zipfile.ZIP_DEFLATED) as z:
            for n, t in S['textos']:
                z.writestr(os.path.basename(n) or 'boletim.txt', t)
        b.seek(0)
        return send_file(b, mimetype='application/zip', as_attachment=True, download_name='boletins_IERS.zip')
    k = int(what[1:]) - 1 if what[:1] == 'g' and what[1:].isdigit() else None
    if k is None or k >= len(res['pngs']):
        abort(404)
    return send_file(io.BytesIO(res['pngs'][k]), mimetype='image/png', as_attachment=True, download_name=f'polo_grafico{k + 1}.png')


# ------------------------------------------------------------------ Editar RINEX (equivalente ao teqc -st/-e, -R -E -C, -O.obs)
RIN = {}


@app.route('/rinex')
def rinex_index():
    return render_template('rinex.html', sid=None)


@app.route('/rinex/processar', methods=['POST'])
def rinex_processar():
    try:
        fm = request.form; arquivos = [(x.filename, x.read()) for x in request.files.getlist('arquivos') if x.filename]
        if not arquivos:
            raise ValueError('envie ao menos um arquivo de observação (.26o, .zip ou .gz)')
        durs = [float(v.replace(',', '.')) for v in re.split(r'[;, ]+', fm.get('duracoes', '5,30,60,120').strip()) if v]
        if not durs:
            raise ValueError('informe as durações em minutos (ex.: 5,30,60,120)')
        inicio = datetime.fromisoformat(fm['inicio']) if fm.get('inicio') else None
        sist = ''.join(c for c in fm.get('sistemas', 'G').upper() if c.isalpha()) or 'G'
        tipos = tuple(t.strip().upper() for t in fm.get('tipos_l1', 'L1,C1').split(',') if t.strip())
        antena = fm.get('antena', '').strip() or None; altura = float(fm['altura'].replace(',', '.')) if fm.get('altura', '').strip() else None
        saidas, entradas = {}, []
        its = [(n, it) for n, b in arquivos for it in observaveis._varrer(b, n)]
        if not its:
            raise ValueError('não encontrei RINEX de observação nos arquivos enviados')
        for n, it in its:
            txt = observaveis._decodifica(it['obs']); base = re.sub(r'\.[A-Za-z0-9]+$', '', os.path.basename(it['nome'] or n))
            pref = RE.nome_seguro(fm['prefixo']) if (fm.get('prefixo') and len(its) == 1) else RE.nome_seguro(re.sub(r'[_\-]?2H$', '', base, flags=re.I))
            R0 = RE.ler(txt); ep = R0['ep']
            m = re.search(r'\.(\d\d)[oOdD]$', it['nome'] or n); ext = (m[1] if m else f'{ep[0]["t"].year % 100:02d}') + 'o'
            entradas.append(dict(nome=base, n=len(ep), ini=ep[0]['t'], fim=ep[-1]['t'], tipos=' '.join(R0['tipos']),
                                 passo=RE.np_median([(b_['t'] - a_['t']).total_seconds() for a_, b_ in zip(ep, ep[1:])]) if len(ep) > 1 else 0))
            for nome, saida, est, desc in RE.sessoes(txt, pref, ext, inicio, durs, 'completa' in fm, 'l1' in fm, 'dupla' in fm, sist, tipos, antena, altura):
                saidas[nome] = (saida.encode('latin-1', 'replace'), est, desc)
        if not saidas:
            raise ValueError('escolha ao menos uma variante de saída')
    except Exception as ex:
        flash(f'Não foi possível processar: {ex}')
        return redirect(url_for('rinex_index'))
    sid = uuid.uuid4().hex[:8]; RIN[sid] = dict(arqs=saidas, entradas=entradas)
    return redirect(url_for('rinex_painel', sid=sid))


@app.route('/rinex/<sid>')
def rinex_painel(sid):
    S = RIN.get(sid) or abort(404)
    linhas = [dict(nome=n, desc=d, ep=e['epocas'], ini=e['inicio'], fim=e['fim'], sats=e['satelites'], tipos=' '.join(e['tipos']), kb=len(b) / 1024)
              for n, (b, e, d) in S['arqs'].items()]
    return render_template('rinex.html', sid=sid, linhas=linhas, entradas=S['entradas'])


@app.route('/rinex/<sid>/baixar/<nome>')
def rinex_baixar(sid, nome):
    S = RIN.get(sid) or abort(404)
    if nome == 'todos.zip':
        b = io.BytesIO()
        with zipfile.ZipFile(b, 'w', zipfile.ZIP_DEFLATED) as z:
            for n, (c, _, _) in S['arqs'].items():
                z.writestr(n, c)
        b.seek(0)
        return send_file(b, mimetype='application/zip', as_attachment=True, download_name='rinex_editados.zip')
    if nome not in S['arqs']:
        abort(404)
    return send_file(io.BytesIO(S['arqs'][nome][0]), mimetype='text/plain', as_attachment=True, download_name=nome)


# ------------------------------------------------------------------ PPP (IBGE): leitura dos resultados
PPP_S = {}


@app.route('/ppp')
def ppp_index():
    return render_template('ppp.html', sid=None)


@app.route('/ppp/enviar', methods=['POST'])
def ppp_enviar():
    try:
        arqs = [(x.filename, x.read()) for x in request.files.getlist('arquivos') if x.filename]
        itens = PPP.carregar(arqs) if arqs else []
        if request.form.get('manual', '').strip():
            itens += PPP.manual(request.form['manual'])
        if not itens:
            raise ValueError('envie os .zip do IBGE-PPP (ou os .sum/.pos) ou preencha a tabela manual')
    except Exception as ex:
        flash(f'Não foi possível ler: {ex}')
        return redirect(url_for('ppp_index'))
    sid = uuid.uuid4().hex[:8]; PPP_S[sid] = dict(itens=itens)
    return redirect(url_for('ppp_painel', sid=sid))


@app.route('/ppp/<sid>', methods=['GET', 'POST'])
def ppp_painel(sid):
    S = PPP_S.get(sid) or abort(404); itens = S['itens']; fm = request.form
    if request.method == 'POST':
        for i, it in enumerate(itens):
            try:
                it['rot'] = dict(cen=int(fm[f'cen{i}']) if fm.get(f'cen{i}') else None, dur=float(fm[f'dur{i}'].replace(',', '.')) if fm.get(f'dur{i}') else None,
                                 modo=fm.get(f'modo{i}') or None)
            except ValueError:
                flash(f'duração inválida na linha {i + 1}')
    L = PPP.linhas(itens)
    com_pos = [i for i, it in enumerate(itens) if it['pos']]
    curva = int(fm.get('curva', com_pos[0] if com_pos else 0)) if itens else 0
    figs = PPP.graficos(L, curva if com_pos else None, itens)
    S['res'] = dict(L=L, pngs=[p for _, p in figs])
    return render_template('ppp.html', sid=sid, itens=itens, tab=PPP.tabela(L), cab=PPP.CAB, figs=[(t, base64.b64encode(p).decode()) for t, p in figs],
                           com_pos=com_pos, curva=curva)


@app.route('/ppp/<sid>/baixar/<what>')
def ppp_baixar(sid, what):
    S = PPP_S.get(sid) or abort(404); res = S.get('res') or abort(404)
    if what == 'csv':
        return send_file(io.BytesIO(PPP.csv_ppp(res['L']).encode('utf-8-sig')), mimetype='text/csv', as_attachment=True, download_name='ppp_resultados.csv')
    k = int(what[1:]) - 1 if what[:1] == 'g' and what[1:].isdigit() else None
    if k is None or k >= len(res['pngs']):
        abort(404)
    return send_file(io.BytesIO(res['pngs'][k]), mimetype='image/png', as_attachment=True, download_name=f'ppp_grafico{k + 1}.png')


# ------------------------------------------------------------------ Sobre (dados do autor)
@app.route('/sobre')
def sobre_pagina():
    return render_template('sobre.html')


@app.route('/diagnostico')
def diagnostico():
    import sys, importlib
    L = [f'Python {sys.version.split()[0]}', f'Pasta do sistema: {BASE}', '']
    for n in modelos.TEMPLATES:
        p = os.path.join(BASE, 'templates', n + '.html'); tam = os.path.getsize(p) if os.path.exists(p) else -1
        L.append(f'templates/{n}.html: {tam} bytes  ' + ('OK' if tam > 1000 else 'PROBLEMA (vazio ou ausente)'))
    L.append('')
    for m in ('nucleo', 'observaveis', 'controle_qualidade', 'crinex', 'tempo', 'polo', 'rinex_edit', 'ppp', 'sobre'):
        try:
            L.append(f'{m}.py: OK ({os.path.basename(importlib.import_module(m).__file__)})')
        except Exception as e:
            L.append(f'{m}.py: PROBLEMA ({e})')
    return '<pre style="font:14px monospace;padding:20px">' + '\n'.join(L) + '</pre>'


if __name__ == '__main__':
    threading.Timer(1, lambda: webbrowser.open('http://127.0.0.1:5000')).start()
    app.run(host='127.0.0.1', port=5000)