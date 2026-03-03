import math


def f(x):
    return x**2 + math.log(x) - 4


def f_vtoraya(x):
    return 2.0 - 1.0 / (x * x)


def f_tretya(x):
    return 2.0 / x**3


# Исход данные
a, b = 1.5, 2.0
n = 10  # число интервалов
h = (b - a) / n
x_star = 1.52


# 1
print("=" * 70)
print("ЭТАП 1. ТАБЛИЦА ЗНАЧЕНИЙ ФУНКЦИИ")
print("-" * 70)
x_uzly = [round(a + i * h, 2) for i in range(n + 1)]
y_uzly = [f(x) for x in x_uzly]

print(f"{'i':3} | {'x_i':6} | {'f(x_i)':12}")
print("-" * 70)
for i, (xi, yi) in enumerate(zip(x_uzly, y_uzly)):
    print(f"{i:3d} | {xi:6.2f} | {yi:12.8f}")
print("=" * 70)

# Определяю индекс i для линейной интерполяции
i_lin = 0
for j in range(len(x_uzly) - 1):
    if x_uzly[j] <= x_star <= x_uzly[j + 1]:
        i_lin = j
        break

x_i = x_uzly[i_lin]
x_i1 = x_uzly[i_lin + 1]
y_i = y_uzly[i_lin]
y_i1 = y_uzly[i_lin + 1]

print(f"\nЛИНЕЙНАЯ ИНТЕРПОЛЯЦИЯ: интервал [{x_i:.2f}, {x_i1:.2f}], индекс i={i_lin}")

# 2.
# МН Лагранжа
L1 = y_i * (x_star - x_i1) / (x_i - x_i1) + y_i1 * (x_star - x_i) / (x_i1 - x_i)
print(f"\nЭТАП 2. L1({x_star}) = {L1:.10f}")

# 3
omega2 = (x_star - x_i) * (x_star - x_i1)
# Значения второй производной на концах интервала
f2_lev = f_vtoraya(x_i)
f2_prav = f_vtoraya(x_i1)

R1_lev = f2_lev * omega2 / 2.0
R1_prav = f2_prav * omega2 / 2.0
min_R1 = min(R1_lev, R1_prav)
max_R1 = max(R1_lev, R1_prav)

# Точное значение функции и фактическая погрешность
y_tochn = f(x_star)
R1_fakt = L1 - y_tochn

print("\nЭТАП 3. ОЦЕНКА ПОГРЕШНОСТИ ЛИНЕЙНОЙ ИНТЕРПОЛЯЦИИ")
print(f"   ω₂(x*) = {omega2:.8f}")
print(f"   f'' на концах: {f2_lev:.8f}, {f2_prav:.8f}")
print(f"   R1 левый = {R1_lev:.8f}, R1 правый = {R1_prav:.8f}")
print(f"   min R1 = {min_R1:.8f}, max R1 = {max_R1:.8f}")
print(f"   Фактическая погрешность R1_fakt = L1 - f(x*) = {R1_fakt:.8f}")
# 4
print("\nЭТАП 4. ОЦЕНКА R и ответ на вопрос пункта 2")
if min_R1 < R1_fakt < max_R1:
    print("   --> Неравенство min R1 < R1_fakt < max R1 ВЫПОЛНЯЕТСЯ")
else:
    print("   --> Неравенство min R1 < R1_fakt < max R1 НЕ ВЫПОЛНЯЕТСЯ")

# Вопрос пункта 2
if abs(L1 - y_tochn) <= 1e-4:
    print("\nЛинейная интерполяция обеспечивает точность 1e-4.")
else:
    print("\nЛинейная интерполяция НЕ обеспечивает точность 1e-4.")
print("-" * 70)


# 5. Квадратичная интерполяция
i_kv = i_lin + 1  # для выполнения условия x_{i-1} < x* < x_i

if i_kv > n - 1:
    # Если x* в последнем интервале, придётся использовать i_kv = n-1 (условие не выполнится)
    i_kv = n - 1
    print("\n Условие x_{i-1}<x*<x_i не выполнено")
# Узлы
x_im1 = x_uzly[i_kv - 1]
x_i0 = x_uzly[i_kv]
x_ip1 = x_uzly[i_kv + 1]
y_im1 = y_uzly[i_kv - 1]
y_i0 = y_uzly[i_kv]
y_ip1 = y_uzly[i_kv + 1]

print(
    f"\nКВАДРАТИЧНАЯ ИНТЕРПОЛЯЦИЯ: узлы x{i_kv-1}={x_im1:.2f}, x{i_kv}={x_i0:.2f}, x{i_kv+1}={x_ip1:.2f}"
)
print(
    f"Условие x_im1 < x* < x_i0 выполнено? {x_im1 < x_star < x_i0 }".replace(
        "True", "ДА"
    ).replace("False", "НЕТ")
)

# Вычисление L2
ch1 = y_im1 * (x_star - x_i0) * (x_star - x_ip1) / ((x_im1 - x_i0) * (x_im1 - x_ip1))
ch2 = y_i0 * (x_star - x_im1) * (x_star - x_ip1) / ((x_i0 - x_im1) * (x_i0 - x_ip1))
ch3 = y_ip1 * (x_star - x_im1) * (x_star - x_i0) / ((x_ip1 - x_im1) * (x_ip1 - x_i0))
L2 = ch1 + ch2 + ch3

print(f"\nЭТАП 5. L2({x_star}) = {L2:.10f}")

# 6
omega3 = (x_star - x_im1) * (x_star - x_i0) * (x_star - x_ip1)
# Значения третьей производной на концах интервала [x_im1, x_ip1]
f3_lev = f_tretya(x_im1)
f3_prav = f_tretya(x_ip1)
R2_lev = f3_lev * omega3 / 6.0
R2_prav = f3_prav * omega3 / 6.0
min_R2 = min(R2_lev, R2_prav)
max_R2 = max(R2_lev, R2_prav)

R2_fakt = L2 - y_tochn

print("\nЭТАП 6. ОЦЕНКА ПОГРЕШНОСТИ КВАДРАТИЧНОЙ ИНТЕРПОЛЯЦИИ")
print(f"   ω₃(x*) = {omega3:.8f}")
print(f"   f''' на концах: {f3_lev:.8f}, {f3_prav:.8f}")
print(f"   R2 левый = {R2_lev:.8f}, R2 правый = {R2_prav:.8f}")
print(f"   min R2 = {min_R2:.8f}, max R2 = {max_R2:.8f}")
print(f"   Фактическая погрешность R2_fakt = L2 - f(x*) = {R2_fakt:.8f}")
# 7
print("\nЭТАП 7.R2 и ответ на вопрос пункта 3")
if min_R2 < R2_fakt < max_R2:
    print("    Неравенство min R2 < R2_fakt < max R2 ВЫПОЛНЯЕТСЯ")
else:
    print("   Неравенство min R2 < R2_fakt < max R2 НЕ ВЫПОЛНЯЕТСЯ")


# Вопрос пункта 3
if abs(L2 - y_tochn) <= 1e-5:
    print("Квадратичная интерполяция обеспечивает точность 1e-5.")
else:
    print("Квадратичная интерполяция НЕ обеспечивает точность 1e-5.")
print("=" * 70)


#  8
print("\nЭТАП 8. ИНТЕРПОЛЯЦИЯ НЬЮТОНА")
print("Таблица разделённых разностей (узлы x_im1, x_i0, x_ip1):")

# Разделённые разности первого порядка
f_im1_i0 = (y_i0 - y_im1) / (x_i0 - x_im1)
f_i0_ip1 = (y_ip1 - y_i0) / (x_ip1 - x_i0)

# Разделённая разность второго порядка
f_im1_i0_ip1 = (f_i0_ip1 - f_im1_i0) / (x_ip1 - x_im1)

print(f"   f[x_im1, x_i0] = {f_im1_i0:.8f}")
print(f"   f[x_i0, x_ip1] = {f_i0_ip1:.8f}")
print(f"   f[x_im1, x_i0, x_ip1] = {f_im1_i0_ip1:.8f}")

# Линейная интерполяция Ньютона по узлам x_i0, x_ip1
N1 = y_i0 + f_i0_ip1 * (x_star - x_i0)
print(f"\n   Линейный многочлен Ньютона (узлы x_i, x_i+1): N1 = {N1:.10f}")
print(f"   Сравнение с L1: L1 = {L1:.10f}, разность = {abs(N1-L1):.2e}")

# Квадратичная интерполяция Ньютона по трём узлам
N2 = (
    y_im1
    + f_im1_i0 * (x_star - x_im1)
    + f_im1_i0_ip1 * (x_star - x_im1) * (x_star - x_i0)
)
print(f"\n   Квадратичный многочлен Ньютона: N2 = {N2:.10f}")
print(f"   Сравнение с L2: L2 = {L2:.10f}, разность = {abs(N2-L2):.2e}")
