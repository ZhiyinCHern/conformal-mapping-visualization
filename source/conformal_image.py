# -*- coding: utf-8 -*-

import argparse
import os
import re
import sys
import numpy as np
from PIL import Image
from scipy.ndimage import map_coordinates
from sympy import symbols, solve, lambdify, I
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations,
                                        implicit_multiplication_application,
                                        convert_xor)

_Z, _W = symbols('z w')
_TR = standard_transformations + (convert_xor, implicit_multiplication_application)
_FUNCS = ['sinh', 'cosh', 'tanh', 'arcsin', 'arccos', 'arctan',
          'sin', 'cos', 'tan', 'exp', 'ln', 'log']


# ========== 第一部分：LaTeX -> 表达式 / 函数 ==========

def _match_brace(s, i):
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
    if i < len(s) and s[i] == '{':
        j = _match_brace(s, i)
        return s[i+1:j], j + 1
    return s[i], i + 1


def latex_to_math(tex):
    """把常见 LaTeX 数学式转换成标准中缀表达式字符串。"""
    s = tex.strip().replace(r'\left', '').replace(r'\right', '')
    s = s.replace(r'\cdot', '*').replace(r'\times', '*')
    s = s.replace(r'\,', '').replace(' ', '')
    while r'\frac' in s:
        k = s.index(r'\frac'); i = k + 5
        A, i = _grab(s, i); B, i = _grab(s, i)
        s = s[:k] + '((' + A + ')/(' + B + '))' + s[i:]
    while r'\sqrt' in s:
        k = s.index(r'\sqrt'); i = k + 5
        A, i = _grab(s, i)
        s = s[:k] + 'sqrt((' + A + '))' + s[i:]
    if 'e^' in s:
        out = ''; i = 0
        while i < len(s):
            if (s[i] == 'e' and i+1 < len(s) and s[i+1] == '^'
                    and (i == 0 or not (s[i-1].isalnum() or s[i-1] in ')}'))):
                A, j = _grab(s, i+2); out += 'exp((' + A + '))'; i = j
            else:
                out += s[i]; i += 1
        s = out
    s = re.sub(r'(?<![A-Za-z])i(?=\\[a-zA-Z])', 'i*', s)
    for fn in _FUNCS:
        s = s.replace('\\' + fn, fn)
    s = s.replace('ln', 'log')
    out = ''; i = 0
    while i < len(s):
        if s[i] == '^':
            A, j = _grab(s, i+1); out += '**(' + A + ')'; i = j
        else:
            out += s[i]; i += 1
    return out.replace('{', '(').replace('}', ')')


def _expr(tex):
    return parse_expr(latex_to_math(tex), transformations=_TR,
                      local_dict={'z': _Z, 'I': I, 'i': I, 'e': np.e})


def _to_callable(expr, var):
    g = lambdify(var, expr, modules=['numpy'])

    def safe(x):
        with np.errstate(all='ignore'):
            v = g(x)
        return np.broadcast_to(np.asarray(v, dtype=complex), np.shape(x))
    return safe


def build_forward_and_inverse(tex):
    """返回 (f 正向函数, g 逆函数, inverted标志)。
    inverted=False 表示求逆失败、g 退回为 f（方向可能相反）。"""
    fe = _expr(tex)
    f = _to_callable(fe, _Z)
    try:
        sols = solve(fe - _W, _Z)
    except Exception:
        sols = []
    if not sols:
        return f, f, False
    # 多值函数会有多支逆，挑一支"最主"的（在几个测试点上模最小）
    tps = np.array([0.5+0.3j, -0.4+0.5j, 0.6-0.2j])

    def score(s):
        gg = _to_callable(s, _W)
        return np.nansum(np.abs(gg(tps)))
    best = min(sols, key=score)
    return f, _to_callable(best, _W), True


# ========== 第二部分：照片变换（直接采样，锐利） ==========

def _auto_limits(v, lo=2.0, hi=98.0, pad=0.08):
    v = v[np.isfinite(v)]
    if v.size == 0:
        return -1.0, 1.0
    a, b = np.percentile(v, [lo, hi])
    if b - a < 1e-9:
        a, b = a - 1, b + 1
    m = (b - a) * pad
    return a - m, b + m


def transform_photo(image, f, g, out_size=600, src_half=2.0, fill="edge"):
    """
    把照片按 w=f(z) 正向变形。
      f, g    : 正向函数与其闭式逆函数（g 用于对每个输出点反查原图坐标）
      src_half: 照片较短边对应复平面 [-src_half, src_half]（整张照片都参与）
      fill    : 越界填充 "edge"(边缘色,默认) / "white" / "black"
    """
    src = np.asarray(image.convert("RGB"))
    H, W = src.shape[:2]
    scale = (min(W, H) / 2) / src_half
    hx, hy = src_half * W / min(W, H), src_half * H / min(W, H)

    # 1) 用 f 在照片区域粗采样，自动决定输出取景范围（忽略奇点极端值）
    ZX, ZY = np.meshgrid(np.linspace(-hx, hx, 160), np.linspace(-hy, hy, 160))
    Wv = f(ZX + 1j*ZY)
    xlo, xhi = _auto_limits(Wv.real)
    ylo, yhi = _auto_limits(Wv.imag)
    cx, cy = (xlo+xhi)/2, (ylo+yhi)/2
    half = max(xhi-xlo, yhi-ylo) / 2
    line = np.linspace(-half, half, out_size)
    OX, OY = np.meshgrid(line + cx, line + cy)

    # 2) 闭式反查每个输出点的原图坐标 z = g(w)，直接在原图采样（锐利）
    Z = g(OX + 1j*OY)
    col = W/2 + Z.real*scale
    row = H/2 - Z.imag*scale
    mode = "nearest" if fill == "edge" else "constant"
    cval = 255 if fill == "white" else 0
    out = np.zeros((out_size, out_size, 3), np.uint8)
    for ch in range(3):
        out[..., ch] = map_coordinates(src[..., ch], [row, col], order=1,
                                       mode=mode, cval=cval)
    return Image.fromarray(out)


# ========== 第三部分：终端交互 ==========

def _clean_path(p):
    return p.strip().strip('"').strip("'").strip()


def _safe_name(tex):
    name = re.sub(r'[^A-Za-z0-9]+', '_', tex).strip('_')
    return name[:30] if name else "result"


def main(argv=None):
    p = argparse.ArgumentParser(description="解析函数对照片的共形变换")
    p.add_argument("--latex", default=None, help="函数 LaTeX（自变量 z）")
    p.add_argument("--image", default=None, help="图片路径")
    p.add_argument("-o", "--out", default=None, help="输出文件名")
    p.add_argument("--fill", default="edge", choices=["edge", "white", "black"])
    args = p.parse_args(argv)

    tex = args.latex
    if tex is None:
        print(r"请输入函数的 LaTeX 代码（自变量用 z），例如：\frac{1}{z}")
        tex = input("f(z) = ").strip()
    if not tex:
        print("没有输入公式，已退出。"); sys.exit(1)

    img_path = args.image
    if img_path is None:
        img_path = input("图片路径：")
    img_path = _clean_path(img_path)
    if not os.path.exists(img_path):
        print(f"找不到图片：{img_path}"); sys.exit(1)

    try:
        f, g, inverted = build_forward_and_inverse(tex)
        f(np.array([1+1j, 0.5-0.3j]))
        g(np.array([1+1j, 0.5-0.3j]))
    except Exception as ex:
        print(f"公式无法解析：{ex}")
        print(r"提示：自变量用 z；分式用 \frac{分子}{分母}；幂用 ^。")
        sys.exit(1)
    if not inverted:
        print("注意：未能求出该函数的闭式逆，变形方向可能与预期相反。")

    print("处理中……")
    img = Image.open(img_path)
    out_img = transform_photo(img, f, g, fill=args.fill)

    here = os.path.dirname(os.path.abspath(__file__))
    stem = os.path.splitext(os.path.basename(img_path))[0]
    out_name = args.out or f"{stem}_{_safe_name(tex)}.png"
    outpath = os.path.join(here, out_name)
    out_img.save(outpath)
    print("完成！")
    print(f"图片已保存到：{outpath}")


if __name__ == "__main__":
    main()
