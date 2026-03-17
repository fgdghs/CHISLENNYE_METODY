import math
import sys


f_output = open("interpolation_results.txt", "w", encoding="utf-8")
sys.stdout = f_output


def f(x):
    return x**2 + math.log(x) - 4


def f_deriv_11(x):
    return math.factorial(10) / (x**11)


def f_deriv_3(x):
    return 2 / (x**3)


a, b = 1.5, 2.0
n_shagov = 10
h = 0.05
# Для более удобного вывода создал список словарей, который будет хранить информацию о каждой точке.
x_points_to_eval = [
    {"name": "x**", "val": 1.52, "LN": "NewTon 1", "deg": 10},
    {"name": "x***", "val": 1.97, "LN": "Gaus 1", "deg": 2},
    {"name": "x****", "val": 1.77, "LN": "Gaus 1", "deg": 10},
]

x_uzly = [round(a + i * h, 2) for i in range(n_shagov + 1)]
y_uzly = [f(x) for x in x_uzly]

# ЭТАП 1. ТАБЛИЦА КОНЕЧНЫХ РАЗНОСТЕЙ
n_points = len(x_uzly)
dy = [[0.0] * n_points for _ in range(n_points)]

for i in range(n_points):
    dy[i][0] = y_uzly[i]

for j in range(1, n_points):
    for i in range(n_points - j):
        dy[i][j] = dy[i + 1][j - 1] - dy[i][j - 1]

print("=" * 160)
print(f"{'ЭТАП 1. ТАБЛИЦА КОНЕЧНЫХ РАЗНОСТЕЙ ':^160}")
print("-" * 160)

header = f"{'xi':>5} | {'f(xi)':>10}"
for j in range(1, n_points):
    header += f" | {'dy^' + str(j):>10}"
print(header)
print("-" * 160)

for line in range(2 * n_points - 1):
    row_str = ""
    if line % 2 == 0:
        idx = line // 2
        row_str += f"{x_uzly[idx]:5.2f} | {dy[idx][0]:10.6f}"
    else:
        row_str += f"{'':>5} | {'':>10}"

    for j in range(1, n_points):
        if (line - j) % 2 == 0:
            idx = (line - j) // 2
            if 0 <= idx < n_points - j:
                val = dy[idx][j]
                fmt = f"{val:10.2e}" if abs(val) < 1e-6 and val != 0 else f"{val:10.6f}"
                row_str += f" | {fmt}"
            else:
                row_str += f" | {'':>10}"
        else:
            row_str += f" | {'':>10}"
    print(row_str)

# ЭТАП 2. ВЫЧИСЛЕНИЕ Ln
print("\n" + "=" * 160)
print(f"{'ЭТАП 2. ВЫЧИСЛЕНИЕ ИНТЕРПОЛЯЦИОННЫХ МНОГОЧЛЕНОВ L(z)':^160}")
print("-" * 160)

# 1. x** (1-я формула Ньютона так как точка расположена в начале, а формула использует узлы от x0 до x10 (разности, идущие «вперёд»))
t2 = (x_points_to_eval[0]["val"] - x_uzly[0]) / h  # из x = x0 + t*h
L_n2 = dy[0][0]
p = 1.0
for k in range(1, 11):
    p *= (t2 - k + 1) / k
    L_n2 += p * dy[0][k]  # вычисления рекурсивный с t умножаем на разности

# 2. x*** (1-я формула Гаусса, так как точка 1.97 ближе к x9 = 1.95, t = 0.4)
idx0_x3 = 9
t3 = (x_points_to_eval[1]["val"] - x_uzly[idx0_x3]) / h
L_n3 = dy[idx0_x3][0]
p = 1.0
for k in range(1, 3):
    if k % 2 != 0:
        p *= (t3 + (k // 2)) / k
    else:
        p *= (t3 - (k // 2)) / k
    L_n3 += p * dy[idx0_x3 - (k // 2)][k]

# 3. x**** (1-я формула Гаусса, x0 = 1.75, так как точка расположена в середине, Формула Гаусса использует узлы,
# симметрично расположенные относительно центра (x5), и разности, идущие в обе стороны от этого центра)
idx0 = 5
t4 = (x_points_to_eval[2]["val"] - x_uzly[idx0]) / h
L_n4 = dy[idx0][0]
p = 1.0
for k in range(1, 11):
    if k % 2 != 0:
        p *= (t4 + (k // 2)) / k
    else:
        p *= (t4 - (k // 2)) / k
    L_n4 += p * dy[idx0 - (k // 2)][k]  # // деление без остатка

x_points_to_eval[0]["L"] = L_n2
x_points_to_eval[1]["L"] = L_n3
x_points_to_eval[2]["L"] = L_n4

for pt in x_points_to_eval:
    print(f"L{pt['deg']}({pt['name']} = {pt['val']}) = {pt['L']:.16f}  ({pt['LN']})")

# ЭТАП 3. ОЦЕНКА ПОГРЕШНОСТИ
print("\n" + "=" * 160)
print(
    f"{'ЭТАП 3 и 4. ОЦЕНКА ПОГРЕШНОСТИ И ПРОВЕРКА НЕРАВЕНСТВА min_Rn <= Rn(z) <= max_Rn':^160}"
)
print("-" * 160)

f11_a = f_deriv_11(a)
f11_b = f_deriv_11(b)

fact11 = math.factorial(11)
fact3 = math.factorial(3)

for pt in x_points_to_eval:
    z = pt["val"]
    Ln_z = pt["L"]
    fz = f(z)
    deg = pt["deg"]

    # Вычисление omega(z)
    omega = 1.0
    if pt["name"] == "x***":
        # Для Гаусса 2-й степени от x9 используются узлы x8, x9, x10
        relevant_nodes = x_uzly[8:11]
    else:
        relevant_nodes = x_uzly[:deg + 1]
        
    for xi in relevant_nodes:
        omega *= z - xi

    if deg == 10:
        val1 = (f11_a * omega) / fact11
        val2 = (f11_b * omega) / fact11
    else:
        # Для L2 используем 3-ю производную
        val1 = (f_deriv_3(a) * omega) / fact3
        val2 = (f_deriv_3(b) * omega) / fact3

    # Min and Max Rn
    min_Rn = min(val1, val2)
    max_Rn = max(val1, val2)

    # Разность Ln(z) - f(z)
    Rn_z = fz - Ln_z

    print(f"\nАНАЛИЗ ДЛЯ ТОЧКИ {pt['name']} ({z}):")
    print(
        f"   L{deg} = {pt['L']:.16f}  ({pt['LN']}), F({pt['val']}) = {f(pt['val']):.16f}, Разница = {abs(pt['L'] - f(pt['val'])):.16e})"
    )
    print(f"   omega_{deg+1}(z) = {omega:.16e}")
    print(
        f"   Min Rn = {min_Rn:.16e} , Max Rn = {max_Rn:.16e}, Max Rn - Min Rn = {abs((max_Rn - min_Rn)):.16e}"
    )
    print(f"   Фактическое Rn(z)     : {Rn_z:.16e}")

    # ЭТАП 4. ПРОВЕРКА НЕРАВЕНСТВА min_Rn <= Rn(z) <= max_Rn

    if (min_Rn) <= Rn_z <= (max_Rn):
        print(f"   РЕЗУЛЬТАТ: Неравенство (min_Rn) <= Rn_z <= (max_Rn) ВЫПОЛНЯЕТСЯ")
    else:
        print(f"   РЕЗУЛЬТАТ: НЕ ВЫПОЛНЯЕТСЯ ")

print("-" * 160)

f_output.close()
sys.stdout = sys.__stdout__

print("READY" * 100)