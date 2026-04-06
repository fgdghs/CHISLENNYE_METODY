import sympy as sp
import numpy as np

# Этап 1
t = sp.symbols("t")
n = 5  # Степень многочлена
k = 2  # Порядок производной
m = 0  # Индекс узла x0

uzly_t = list(range(n + 1))
koeff_A = []

for i in uzly_t:
    # Числитель
    chislitel = 1
    for j in uzly_t:
        if i != j:
            chislitel *= t - j

    # Знаменатель
    znamenatel = 1
    for j in uzly_t:
        if i != j:
            znamenatel *= i - j

    l_i = chislitel / znamenatel
    sp.pprint(l_i, use_unicode=True)
    # Дифференцируем  k раз и подставляем точку m
    proizv_sym = sp.diff(l_i, t, k)
    sp.pprint(proizv_sym, use_unicode=True)
    znachenie = proizv_sym.subs(t, m)
    koeff_A.append(float(znachenie))

koeff_A = np.array(koeff_A)

# Этап 2
h = 0.1  # Шаг сетки
x0 = 1.5  # Начальная точка

# Исходная функция и её производные
f = lambda x: x**2 + np.log(x) - 4
f_2_tochno = lambda x: 2 - (1 / x**2)
f_6_tochno = lambda x: -120 * (x**-6)

uzly_x = np.array([x0 + i * h for i in range(n + 1)])
y_vals = f(uzly_x)

# Расчет значения производной по формуле Лагранжа
priblizh_val = np.sum(koeff_A * y_vals) / (h**k)
tochnoe_val = f_2_tochno(x0)
fakt_oshibka = tochnoe_val - priblizh_val

# Этап 3
# Вычисление константы omega
t_sym = sp.symbols("t_sym")
omega = 1
for i in range(n + 1):
    omega *= t_sym - i

omega_diff_k = sp.diff(omega, t_sym, k).subs(t_sym, m)
teor_const = float(omega_diff_k) / sp.factorial(n + 1)

# Поиск экстремумов (n+1)-й производной на отрезке
xi_vals = np.linspace(uzly_x[0], uzly_x[-1], 100)
proizv_6_vals = f_6_tochno(xi_vals)

R_min = teor_const * (h ** (n + 1 - k)) * np.min(proizv_6_vals)
R_max = teor_const * (h ** (n + 1 - k)) * np.max(proizv_6_vals)

if R_min > R_max:
    R_min, R_max = R_max, R_min


print(f"--- Результаты ---")
print(f"Приближенное L''(x0): {priblizh_val:.10f}")
print(f"Точное f''(x0):      {tochnoe_val:.10f}")
print(f"Фактическая ошибка:  {fakt_oshibka:.2e}")
print("-" * 45)
print(f"Теоретический диапазон погрешности:")
print(f"min Rn,k: {R_min:.2e}")
print(f"max Rn,k: {R_max:.2e}")
print(f"Условие min < R < max: {R_min <= fakt_oshibka <= R_max}")
