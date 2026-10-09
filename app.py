# -*- coding: utf-8 -*-
"""Sistema local (Flask): efemérides transmitidas x precisas (SP3).
Execute:  python app.py   -> abre http://127.0.0.1:5000"""
import base64, csv, gzip, io, os, re, threading, uuid, webbrowser, zipfile
from datetime import datetime
import numpy as np
from flask import Flask, render_template, request, redirect, url_for, send_file, abort, flash
import nucleo, observaveis, teqc

BASE = os.path.dirname(os.path.abspath(__file__))
app = Flask(__name__)
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
        arqs = [x for x in request.files.getlist('arquivos') if x.filename]
        if not arqs:
            raise ValueError('envie ao menos um arquivo')
        fs = [observaveis.carregar(x.read(), x.filename) for x in arqs]
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
        arqs = [x for x in request.files.getlist('arquivos') if x.filename]
        if not arqs:
            raise ValueError('envie ao menos um arquivo de observação')
        nv = request.files.get('nav')
        nav_txt = _texto(nv.read(), 'nav') if nv and nv.filename else None
        fs = [teqc.carregar(x.read(), x.filename, nav_txt) for x in arqs]
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
        ns = np.isfinite(f['A'][:, :, 0]).sum(0)
        sat = fm.get('sat') if fm.get('sat') in f['sats'] else f['sats'][int(np.argmax(ns))]
        figs = teqc.graficos(f, R, sat, par['mask'])
        res = dict(resumo=[teqc.linha_res(r) for r in Rs], sats=teqc.linhas_sat(R), figs=[(t, base64.b64encode(p).decode()) for t, p in figs],
                   pngs=[p for _, p in figs], sat=sat, par=par, idx=idx, tem_el=R['tem_el'], nome=f"{f['estacao']} · {f['nome']}",
                   info=dict(rec=f.get('rec') or '—', ant=f.get('ant') or '—'))
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
        cab, rows = (teqc.CAB_RES, res['resumo']) if what == 'resumo' else (teqc.CAB_SAT, res['sats'])
        txt = ','.join(cab) + '\n' + '\n'.join(','.join(map(str, r)) for r in rows)
        return send_file(io.BytesIO(txt.encode('utf-8-sig')), mimetype='text/csv', as_attachment=True, download_name=f'qc_{what}.csv')
    k = int(what[1:]) - 1 if what[:1] == 'g' and what[1:].isdigit() else None
    if k is None or k >= len(res['pngs']):
        abort(404)
    return send_file(io.BytesIO(res['pngs'][k]), mimetype='image/png', as_attachment=True, download_name=f'qc_grafico{k + 1}.png')


if __name__ == '__main__':
    threading.Timer(1, lambda: webbrowser.open('http://127.0.0.1:5000')).start()
    app.run(host='127.0.0.1', port=5000)