# -*- coding: utf-8 -*-

import argparse
import os
import re
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import mathtext
from sympy import lambdify, Symbol, I
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations,
                                        implicit_multiplication_application,
                                        convert_xor)

# ========== 第一部分：把 LaTeX 变成可计算的函数 ==========

_FUNCS = ['sinh', 'cosh', 'tanh', 'arcsin', 'arccos', 'arctan',
          'sin', 'cos', 'tan', 'exp', 'ln', 'log']


def _match_brace(s, i):
    """s[i] 是 '{'，返回配对的 '}' 的下标。"""
    depth = 0
    for j in range(i, len(s)):
        if s[j] == '{':
            depth += 1
        elif s[j] == '}':
            depth -= 1
            if depth == 0:
                return j
    raise ValueError("LaTeX 花括号不匹配")


def _grab(s, i):
    """从下标 i 取一个参数：{...} 返回内部，否则返回单个字符。"""
    if i < len(s) and s[i] == '{':
        j = _match_brace(s, i)
        return s[i+1:j], j + 1
    return s[i], i + 1


def latex_to_math(tex):
    """把常见 LaTeX 数学式转换成标准中缀表达式字符串。"""
    s = tex.strip().replace(r'\left', '').replace(r'\right', '')
    s = s.replace(r'\cdot', '*').replace(r'\times', '*')
    s = s.replace(r'\,', '').replace(' ', '')

    # \frac{A}{B} -> ((A)/(B))，支持嵌套
    while r'\frac' in s:
        k = s.index(r'\frac'); i = k + 5
        A, i = _grab(s, i)
        B, i = _grab(s, i)
        s = s[:k] + '((' + A + ')/(' + B + '))' + s[i:]

    # \sqrt{A} -> sqrt((A))
    while r'\sqrt' in s:
        k = s.index(r'\sqrt'); i = k + 5
        A, i = _grab(s, i)
        s = s[:k] + 'sqrt((' + A + '))' + s[i:]

    # e^{A} / e^A -> exp((A))
    if 'e^' in s:
        out = ''; i = 0
        while i < len(s):
            if (s[i] == 'e' and i+1 < len(s) and s[i+1] == '^'
                    and (i == 0 or not (s[i-1].isalnum() or s[i-1] in ')}'))):
                A, j = _grab(s, i+2)
                out += 'exp((' + A + '))'; i = j
            else:
                out += s[i]; i += 1
        s = out

    # 虚数单位 i 紧贴函数时插入乘号：i\sin -> i*\sin
    s = re.sub(r'(?<![A-Za-z])i(?=\\[a-zA-Z])', 'i*', s)

    # 函数名去掉反斜杠；ln 视作自然对数
    for fn in _FUNCS:
        s = s.replace('\\' + fn, fn)
    s = s.replace('ln', 'log')

    # ^{A} -> **(A)
    out = ''; i = 0
    while i < len(s):
        if s[i] == '^':
            A, j = _grab(s, i+1)
            out += '**(' + A + ')'; i = j
        else:
            out += s[i]; i += 1
    s = out

    return s.replace('{', '(').replace('}', ')')


def make_function(tex):
    """LaTeX 字符串 -> 可对 numpy 复数数组计算的函数 f。"""
    math = latex_to_math(tex)
    tr = standard_transformations + (convert_xor,
                                     implicit_multiplication_application)
    expr = parse_expr(math, transformations=tr,
                      local_dict={'z': Symbol('z'), 'I': I, 'i': I, 'e': np.e})
    f = lambdify(Symbol('z'), expr, modules=['numpy'])

    def safe_f(z):
        with np.errstate(all='ignore'):
            w = f(z)
        return np.broadcast_to(np.asarray(w, dtype=complex), np.shape(z))
    return safe_f, math


# ========== 第二部分：画网格 + 自动调整取景 ==========

def _auto_limits(values, lo=2.0, hi=98.0, pad=0.08):
    """根据数据 2%/98% 分位数定出主体范围，忽略奇点处极端值。"""
    v = values[np.isfinite(values)]
    if v.size == 0:
        return -1.0, 1.0
    a, b = np.percentile(v, [lo, hi])
    if b - a < 1e-9:
        a, b = a - 1, b + 1
    m = (b - a) * pad
    return a - m, b + m


def _make_grid_lines(grid):
    """返回若干 (z数组, 颜色)，构成初始网格。"""
    lines = []
    if grid == 'rect':
        xr, yr, n, m = (-2, 2), (-2, 2), 21, 400
        for x in np.linspace(*xr, n):
            lines.append((x + 1j*np.linspace(*yr, m), "#1f6feb"))
        for y in np.linspace(*yr, n):
            lines.append((np.linspace(*xr, m) + 1j*y, "#d1242f"))
    else:
        r_min, r_max, n_c, n_ray, m = 0.3, 2.0, 9, 24, 400
        theta = np.linspace(0, 2*np.pi, m)
        for r in np.linspace(r_min, r_max, n_c):
            lines.append((r*np.exp(1j*theta), "#1f6feb"))
        radii = np.linspace(r_min, r_max, m)
        for a in np.linspace(0, 2*np.pi, n_ray, endpoint=False):
            lines.append((radii*np.exp(1j*a), "#d1242f"))
    return lines


def _title_text(tex):
    """优先用数学体标题；LaTeX 不被 mathtext 支持时退回纯文本。"""
    cand = f"$w = {tex}$"
    try:
        mathtext.MathTextParser('agg').parse(cand)
        return cand
    except Exception:
        return f"w = {tex}"


def visualize(f, tex, grid, outpath):
    lines = _make_grid_lines(grid)
    ws = [(z, f(z), color) for z, color in lines]
    all_w = np.concatenate([w for _, w, _ in ws])
    xlo, xhi = _auto_limits(all_w.real)
    ylo, yhi = _auto_limits(all_w.imag)
    cx, cy = (xlo+xhi)/2, (ylo+yhi)/2
    half = max(xhi-xlo, yhi-ylo) / 2
    cap = half * 3

    fig, (az, aw) = plt.subplots(1, 2, figsize=(12, 6))
    for z, w, color in ws:
        az.plot(z.real, z.imag, color=color, lw=0.8)
        wr, wi = w.real.copy(), w.imag.copy()
        bad = ~np.isfinite(wr) | ~np.isfinite(wi) | \
              (np.abs(wr-cx) > cap) | (np.abs(wi-cy) > cap)
        wr[bad] = np.nan; wi[bad] = np.nan
        aw.plot(wr, wi, color=color, lw=0.8)

    for ax, name in [(az, "z-plane (original grid)"),
                     (aw, "w-plane (after f)")]:
        ax.set_aspect("equal")
        ax.axhline(0, color="gray", lw=0.5, zorder=0)
        ax.axvline(0, color="gray", lw=0.5, zorder=0)
        ax.set_title(name, fontsize=12)
        ax.grid(alpha=0.2)
    aw.set_xlim(cx-half, cx+half)
    aw.set_ylim(cy-half, cy+half)

    fig.suptitle(_title_text(tex), fontsize=15, fontweight="bold")
    fig.tight_layout()
    fig.savefig(outpath, dpi=130, bbox_inches="tight")
    plt.close(fig)


# ========== 第三部分：终端交互 ==========

def _safe_name(tex):
    name = re.sub(r'[^A-Za-z0-9]+', '_', tex).strip('_')
    return ("map_" + name)[:40] if name else "map_result"


def main(argv=None):
    p = argparse.ArgumentParser(description="解析函数对网格的共形变换可视化")
    p.add_argument("--latex", default=None, help="函数 LaTeX（自变量 z），非交互用")
    p.add_argument("--grid", default=None, choices=["rect", "polar"],
                   help="初始网格：rect(直角) / polar(极坐标)")
    p.add_argument("-o", "--out", default=None, help="输出文件名")
    args = p.parse_args(argv)

    tex = args.latex
    if tex is None:
        print(r"请输入函数的 LaTeX 代码（自变量用 z），例如：z^2 + \frac{1}{z}")
        tex = input("f(z) = ").strip()
    if not tex:
        print("没有输入公式，已退出。"); sys.exit(1)

    grid = args.grid
    if grid is None:
        ans = input("初始网格  [1] 直角坐标   [2] 极坐标   （默认 1）：").strip()
        grid = "polar" if ans == "2" else "rect"

    try:
        f, math = make_function(tex)
        f(np.array([1+1j, 0.5-0.3j]))
    except Exception as ex:
        print(f"公式无法解析：{ex}")
        print(r"提示：自变量用 z；分式用 \frac{分子}{分母}；幂用 ^。")
        sys.exit(1)

    here = os.path.dirname(os.path.abspath(__file__))
    out_name = args.out or (_safe_name(tex) + ".png")
    outpath = os.path.join(here, out_name)

    visualize(f, tex, grid, outpath)
    print(f"完成！（解析为 {math}）")
    print(f"图片已保存到：{outpath}")


if __name__ == "__main__":
    main()
