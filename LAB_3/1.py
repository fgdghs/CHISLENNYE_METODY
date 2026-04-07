import sympy as sp
import numpy as np

x, h_sym, x0_sym = sp.symbols("x h x0")
n = 5  
k = 2  
h_val = 0.1
x0_val = 1.5

f_syms = sp.symbols(f'f0:{n+1}')
uzly_x_sym = [x0_sym + i * h_sym for i in range(n + 1)]

L_n_x = 0
for i in range(n + 1):
    chislitel = 1
    znamenatel = 1
    for j in range(n + 1):
        if i != j:
            chislitel *= (x - uzly_x_sym[j])
            znamenatel *= (uzly_x_sym[i] - uzly_x_sym[j])
    
    L_n_x += f_syms[i] * (chislitel / znamenatel)

L_n_deriv2_x = sp.diff(L_n_x, x, k)
formula_at_x0 = sp.simplify(L_n_deriv2_x.subs(x, x0_sym))

print("="*70)
print("АНАЛИТИЧЕСКАЯ ФОРМУЛА L''(x0):")
print("="*70)
sp.pprint(formula_at_x0, use_unicode=True)
print("="*70 + "\n")

f_func = lambda v: v**2 + np.log(v) - 4
f_2_tochno = lambda v: 2 - (1 / v**2)
f_6_tochno = lambda v: -120 * (v**-6)

uzly_x_num = np.array([x0_val + i * h_val for i in range(n + 1)])
y_vals = f_func(uzly_x_num)

subs_dict = {h_sym: h_val}
subs_dict.update({f_syms[i]: y_vals[i] for i in range(n+1)})

priblizh_val = float(formula_at_x0.subs(subs_dict))
tochnoe_val = f_2_tochno(x0_val)
fakt_oshibka = tochnoe_val - priblizh_val

x_num = sp.symbols('x_num')
omega = 1
for xi in uzly_x_num:
    omega *= (x_num - xi)

omega_diff_k = sp.diff(omega, x_num, k).subs(x_num, x0_val)
teor_const = float(omega_diff_k) / sp.factorial(n + 1)

xi_vals = np.linspace(uzly_x_num[0], uzly_x_num[-1], 100)
proizv_6_vals = f_6_tochno(xi_vals)

R_min = teor_const * np.min(proizv_6_vals)
R_max = teor_const * np.max(proizv_6_vals)

if R_min > R_max: R_min, R_max = R_max, R_min

print(f"Приближенное L''(x0): {priblizh_val:.16f}")
print(f"Точное f''(x0):      {tochnoe_val:.16f}")
print(f"Фактическая ошибка:  {fakt_oshibka:.16e}")
print("-" * 45)
print(f"Диапазон погрешности: [{R_min:.16e}, {R_max:.16e}]")
print(f"Попадает в диапазон: {R_min <= fakt_oshibka <= R_max}")