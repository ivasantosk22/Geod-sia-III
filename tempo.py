# -*- coding: utf-8 -*-
"""Sistemas de tempo (UTC, UT1, TAI, TGPS, TT, TCG), Data Juliana, MJD e semana GPS, conforme a Aula 7."""
from datetime import date, datetime, timedelta

LG = 6.969290134e-10                     # TCG - TT = Lg * (MJD - 43144,0) * 86400 s  (MJD do TAI)
# TAI - UTC (segundos intercalados, IERS Boletim C): (a partir de, valor)
SALTOS = [(date(1972, 1, 1), 10), (date(1972, 7, 1), 11), (date(1973, 1, 1), 12), (date(1974, 1, 1), 13), (date(1975, 1, 1), 14),
          (date(1976, 1, 1), 15), (date(1977, 1, 1), 16), (date(1978, 1, 1), 17), (date(1979, 1, 1), 18), (date(1980, 1, 1), 19),
          (date(1981, 7, 1), 20), (date(1982, 7, 1), 21), (date(1983, 7, 1), 22), (date(1985, 7, 1), 23), (date(1988, 1, 1), 24),
          (date(1990, 1, 1), 25), (date(1991, 1, 1), 26), (date(1992, 7, 1), 27), (date(1993, 7, 1), 28), (date(1994, 7, 1), 29),
          (date(1996, 1, 1), 30), (date(1997, 7, 1), 31), (date(1999, 1, 1), 32), (date(2006, 1, 1), 33), (date(2009, 1, 1), 34),
          (date(2012, 7, 1), 35), (date(2015, 7, 1), 36), (date(2017, 1, 1), 37)]
ULTIMO_SALTO = SALTOS[-1][0]
DIAS = ['Domingo', 'Segunda-feira', 'Terça-feira', 'Quarta-feira', 'Quinta-feira', 'Sexta-feira', 'Sábado']


def tai_utc(d):
    """TAI-UTC (s) na data d, ou None se anterior a 1972."""
    v = None
    for ini, s in SALTOS:
        if d >= ini:
            v = s
    return v


def jd_data(Y, M, D):
    """Data Juliana às 0h (Leick, 1995), com divisões inteiras como na Aula 7."""
    a = int(7 * (Y + int((M + 9) / 12)) / 4)
    b = int(3 * (int((Y + int((M - 9) / 7)) / 100) + 1) / 4)       # (M-9)/7 é divisão INTEIRA (trunca), como na fórmula original
    return 367 * Y - a - b + int(275 * M / 9) + D + 1721028.5


def _seg(dt):
    return dt.hour * 3600 + dt.minute * 60 + dt.second + dt.microsecond / 1e6


def fmt_dt(dt):
    s = f'{dt.second + dt.microsecond / 1e6:09.6f}'.rstrip('0').rstrip('.')
    s = s if '.' in s or len(s) > 2 else s.zfill(2)
    return f'{dt:%d/%m/%Y} {dt.hour:02d}:{dt.minute:02d}:{s.zfill(2) if "." not in s else s}'


def calcular(hl, fuso_h=3.0, dut1=0.0, tai_utc_s=None):
    """hl: hora local (datetime). UTC = HL + fuso_h (Belém: +3). Retorna dict com tudo."""
    utc = hl + timedelta(hours=fuso_h)
    saltos = tai_utc(utc.date()) if tai_utc_s is None else tai_utc_s
    if saltos is None:
        raise ValueError('TAI-UTC não está definido antes de 1972: informe o valor manualmente')
    ut1 = utc + timedelta(seconds=dut1)
    tai = utc + timedelta(seconds=saltos)
    tgps = tai - timedelta(seconds=19)
    tt = tai + timedelta(seconds=32.184)
    jd0 = jd_data(utc.year, utc.month, utc.day)
    jd = jd0 + _seg(utc) / 86400
    mjd0, mjd = jd0 - 2400000.5, jd - 2400000.5
    mjd_tai = mjd + saltos / 86400
    tcg_tt = LG * (mjd_tai - 43144.0) * 86400
    tcg = tt + timedelta(seconds=tcg_tt)
    jdg = jd_data(tgps.year, tgps.month, tgps.day)                 # semana/dia GPS pela data em tempo GPS
    n = int(round(jdg - 2444244.5))
    semana, dow = n // 7, n % 7
    return dict(hl=hl, utc=utc, ut1=ut1, tai=tai, tgps=tgps, tt=tt, tcg=tcg, tcg_tt=tcg_tt, jd0=jd0, jd=jd, mjd0=mjd0, mjd=mjd,
                semana=semana, dow=dow, doy=utc.timetuple().tm_yday, sow=dow * 86400 + _seg(tgps), saltos=saltos, dut1=dut1, fuso=fuso_h,
                sp3_curto=f'igs{semana}{dow}.sp3', sp3_longo=f'IGS0OPSFIN_{utc.year}{utc.timetuple().tm_yday:03d}0000_01D_15M_ORB.SP3')


def tabela(r):
    """Linhas (sistema, equação, resultado) no formato da tabela da aula."""
    f = fmt_dt
    return [
        ('Hora local', 'dado de entrada', f(r['hl'])),
        ('UTC', f'UTC = HL + {r["fuso"]:g} h (Belém: +3 h)', f(r['utc'])),
        ('UT1', f'UT1 = UTC + DUT1  (DUT1 = {r["dut1"]:+g} s)', f(r['ut1'])),
        ('TAI', f'TAI = UTC + {r["saltos"]} s (TAI−UTC)', f(r['tai'])),
        ('TGPS', 'TGPS = TAI − 19 s', f(r['tgps'])),
        ('TT = TDT', 'TT = TAI + 32,184 s', f(r['tt'])),
        ('TCG', f'TCG = TT + Lg·(MJD_TAI − 43144,0)·86400  (= TT + {r["tcg_tt"]:.6f} s)', f(r['tcg'])),
        ('JD (0h UTC)', 'Leick (1995): fórmula da aula', f'{r["jd0"]:.1f}'),
        ('JD (com a hora)', 'JD = JD(0h) + hora UTC/24', f'{r["jd"]:.6f}'),
        ('MJD (0h UTC)', 'MJD = JD − 2400000,5', f'{r["mjd0"]:.0f}'),
        ('MJD (com a hora)', 'MJD = JD − 2400000,5', f'{r["mjd"]:.6f}'),
        ('Semana GPS', 'int[(JD − 2444244,5)/7]', str(r['semana'])),
        ('Dia da semana GPS', '0 = domingo … 6 = sábado', f'{r["dow"]}  ({DIAS[r["dow"]]})'),
        ('Dia do ano (DOY)', 'data em UTC', f'{r["doy"]:03d}'),
        ('Segundos da semana GPS', 'dia·86400 + segundos do dia em TGPS', f'{r["sow"]:.3f}'),
        ('Arquivo SP3 (IGS final)', 'nomes curto e longo', f'{r["sp3_curto"]}  |  {r["sp3_longo"]}.gz'),
    ]


# Exemplo resolvido da aula: Belém, 05/07/2012 09:10:20, DUT1 = +0,4 s, TAI−UTC = 35 s
EXEMPLO = dict(hl=datetime(2012, 7, 5, 9, 10, 20), fuso=3.0, dut1=0.4, tai=35)
ESPERADO = {'UTC': '05/07/2012 12:10:20', 'UT1': '05/07/2012 12:10:20.4', 'TAI': '05/07/2012 12:10:55', 'TGPS': '05/07/2012 12:10:36',
            'TT = TDT': '05/07/2012 12:11:27.184', 'JD (0h UTC)': '2456113.5', 'MJD (0h UTC)': '56113', 'Semana GPS': '1695', 'Dia da semana GPS': '4  (Quinta-feira)'}