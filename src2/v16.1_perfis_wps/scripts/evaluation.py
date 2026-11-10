"""Resumos de pares com pesos explícitos; sem supor réplicas por par."""
import numpy as np
import pandas as pd


def concordance(same, cross, same_weights=None, cross_weights=None):
    s, c = np.asarray(same, float), np.asarray(cross, float)
    if len(s) == 0 or len(c) == 0 or not np.isfinite(s).all() or not np.isfinite(c).all():
        raise ValueError('Escores finitos de ambos os grupos são necessários.')
    sw = np.ones(len(s)) if same_weights is None else np.asarray(same_weights, float)
    cw = np.ones(len(c)) if cross_weights is None else np.asarray(cross_weights, float)
    order = np.greater(s[:, None], c[None, :]).astype(float) + .5*np.equal(s[:, None], c[None, :])
    return float((sw[:, None]*cw[None, :]*order).sum()/(sw.sum()*cw.sum()))


def group_weights(groups):
    g = pd.Series(np.asarray(groups))
    return (1/g.map(g.value_counts())).to_numpy(float)


def summarize(frame):
    same = frame[frame.mesmo_padrao]
    cross = frame[~frame.mesmo_padrao]
    sw, cw = group_weights(same.grupo), group_weights(cross.grupo)
    return dict(pares_intra=len(same), pares_entre=len(cross),
                grupos_intra=same.grupo.nunique(), grupos_entre=cross.grupo.nunique(),
                A_pares=concordance(same.similaridade, cross.similaridade),
                A_balanceada=concordance(same.similaridade, cross.similaridade, sw, cw),
                media_intra=float(same.valor.mean()), media_entre=float(cross.valor.mean()),
                media_intra_balanceada=float(np.average(same.valor, weights=sw)),
                media_entre_balanceada=float(np.average(cross.valor, weights=cw)))


def summarize_pattern(frame, pattern):
    same = frame[frame.mesmo_padrao & (frame.padrao_x == pattern)]
    cross = frame[~frame.mesmo_padrao & ((frame.padrao_x == pattern) | (frame.padrao_y == pattern))]
    if len(same) == 0:
        return dict(pares_intra=0, pares_entre=len(cross), A_pares=np.nan, A_balanceada=np.nan)
    return dict(pares_intra=len(same), pares_entre=len(cross),
                A_pares=concordance(same.similaridade, cross.similaridade),
                A_balanceada=concordance(same.similaridade, cross.similaridade,
                                        cross_weights=group_weights(cross.grupo)))
