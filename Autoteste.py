# -*- coding: utf-8 -*-
"""Autoteste do sistema no SEU computador.
Uso:  python autoteste.py dados_BEPA.zip     (ou qualquer .zip da RBMC)
Confere: arquivos/versões corretos, leitura do Hatanaka e o controle de qualidade."""
import importlib, os, sys, time


def ok(msg, cond):
    print(('[OK]    ' if cond else '[FALHA] ') + msg)
    return cond


def main(arq=None):
    print('Python', sys.version.split()[0], '| pasta:', os.getcwd())
    todos = True
    for m in ('nucleo', 'crinex', 'observaveis', 'controle_qualidade', 'app'):
        try:
            importlib.import_module(m); todos &= ok(f'módulo {m}', True)
        except Exception as e:
            todos &= ok(f'módulo {m}: {e}', False)
    if not todos:
        return print('\nCorrija os arquivos acima (versão antiga ou faltando) e rode de novo.')
    import observaveis as O, controle_qualidade as Q
    todos &= ok('observaveis.carregar_todos existe', hasattr(O, 'carregar_todos'))
    todos &= ok('controle_qualidade.carregar_todos existe', hasattr(Q, 'carregar_todos'))
    if arq:
        t = time.time(); fs = Q.carregar_todos(open(arq, 'rb').read(), os.path.basename(arq))
        todos &= ok(f'{len(fs)} arquivo(s) de observação lido(s) em {time.time() - t:.0f}s', len(fs) > 0)
        for f in fs:
            r = Q.analisar(f)['res']
            print(f"   {f['estacao']} {f['ep'][0]:%d/%m/%Y}: {r['nep']} épocas, {r['nsat']} sats, completude {r['compl']:.1f}%, MP1 {r['mp1']:.2f} m, MP2 {r['mp2']:.2f} m, saltos {r['slips']}")
        if 'BEPA' in os.path.basename(arq).upper() and len(fs) == 1:      # valores de referência obtidos nos testes
            r = Q.analisar(fs[0])['res']
            ref = (5760, 31, 88.1, 0.38, 0.33, 29); got = (r['nep'], r['nsat'], round(r['compl'], 1), round(r['mp1'], 2), round(r['mp2'], 2), r['slips'])
            todos &= ok(f'BEPA bate com os valores de referência {ref}', got == ref)
            if got != ref:
                print('   obtido:', got)
    print('\n' + ('TUDO CERTO.' if todos else 'HÁ FALHAS: me envie esta saída.'))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)