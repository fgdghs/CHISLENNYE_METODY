import numpy as np

# Я определяю функцию из варианта №3 и её вторую производную для краевых условий
def f(x): 
    return x**2 + np.log(x) - 4

def f_vtoraya_proizvodnaya(x): 
    return 2.0 - 1.0 / (x**2)

# ВАШ МЕТОД ПРОГОНКИ
def metod_progonki(A, B, C, F):
    n = len(B)
    alfa = np.zeros(n)
    beta = np.zeros(n)

    # ПРЯМОЙ ХОД ПРОГОНКИ
    alfa[0] = -C[0] / B[0]
    beta[0] = F[0] / B[0]

    for i in range(1, n - 1):
        znamenatel = B[i] + A[i] * alfa[i - 1]
        alfa[i] = -C[i] / znamenatel
        beta[i] = (F[i] - A[i] * beta[i - 1]) / znamenatel

    # ОБРАТНЫЙ ХОД ПРОГОНКИ
    reshenie_U = np.zeros(n)
    znamenatel_n = B[n - 1] + A[n - 1] * alfa[n - 2]
    reshenie_U[n - 1] = (F[n - 1] - A[n - 1] * beta[n - 2]) / znamenatel_n

    for i in range(n - 2, -1, -1):
        reshenie_U[i] = alfa[i] * reshenie_U[i + 1] + beta[i]

    return reshenie_U

# Я подготавливаю сетку и значения функции
x_tochki = np.linspace(1.5, 2.0, 5)
y = f(x_tochki)
N = len(x_tochki) - 1
h = np.diff(x_tochki)  # Шаги сетки h_i

# Краевые условия 2-го типа по методичке: M0 = A, MN = B
M_0 = f_vtoraya_proizvodnaya(x_tochki[0])
M_N = f_vtoraya_proizvodnaya(x_tochki[-1])

# Инициализация массивов для системы из (N-1) уравнений
n_vnutr = N - 1
koeff_A = np.zeros(n_vnutr)
koeff_B = np.zeros(n_vnutr)
koeff_C = np.zeros(n_vnutr)
pravaya_chast_F = np.zeros(n_vnutr)

# --- ПЕРВАЯ СТРОКА СИСТЕМЫ (i = 1) ВРУЧНУЮ ---
# Соответствует уравнению: 2*M1 + mu_1*M2 = d_1 - lambda_1*M0
h1, h2 = h[0], h[1]
lyambda_1 = h1 / (h1 + h2)  #
myu_1 = h2 / (h1 + h2)      #
d_1 = (6.0 / (h1 + h2)) * ((y[2] - y[1])/h2 - (y[1] - y[0])/h1) #

koeff_B[0] = 2.0
koeff_C[0] = myu_1
pravaya_chast_F[0] = d_1 - lyambda_1 * M_0


# --- СРЕДНИЕ СТРОКИ СИСТЕМЫ (от i = 2 до N-2) ЦИКЛОМ ---
# Заполняются строго по формуле: lambda_i*M_{i-1} + 2*M_i + mu_i*M_{i+1} = d_i
for i in range(2, N - 1):
    idx = i - 1
    hi, hi1 = h[i-1], h[i]
    
    koeff_A[idx] = hi / (hi + hi1)        # lyambda_i
    koeff_B[idx] = 2.0
    koeff_C[idx] = hi1 / (hi + hi1)       # myu_i
    pravaya_chast_F[idx] = (6.0 / (hi + hi1)) * ((y[i+1] - y[i])/hi1 - (y[i] - y[i-1])/hi) #


# --- ПОСЛЕДНЯЯ СТРОКА СИСТЕМЫ (i = N-1) ВРУЧНУЮ ---
# Соответствует уравнению: lambda_{N-1}*M_{N-2} + 2*M_{N-1} = d_{N-1} - mu_{N-1}*MN
idx_last = n_vnutr - 1
h_n_1, h_n = h[N-2], h[N-1]
lyambda_n_1 = h_n_1 / (h_n_1 + h_n)
myu_n_1 = h_n / (h_n_1 + h_n)
d_n_1 = (6.0 / (h_n_1 + h_n)) * ((y[N] - y[N-1])/h_n - (y[N-1] - y[N-2])/h_n_1)

koeff_A[idx_last] = lyambda_n_1
koeff_B[idx_last] = 2.0
pravaya_chast_F[idx_last] = d_n_1 - myu_n_1 * M_N


# Я нахожу внутренние моменты, используя ВАШ метод прогонки
M_vnutrennie = metod_progonki(koeff_A, koeff_B, koeff_C, pravaya_chast_F)

# Я собираю итоговый массив всех моментов M_i
M = np.concatenate(([M_0], M_vnutrennie, [M_N]))

# Вывод результатов
print("Моменты M_i (вторые производные в узлах):")
for i, val in enumerate(M):
    print(f"M_{i} = {val:.6f}")

# Я рассчитываю значение S''(x) в контрольной точке x**** = 1.77
x_test = 1.77
for i in range(N):
    if x_tochki[i] <= x_test <= x_tochki[i+1]:
        h_segment = h[i]
        # Формула из методички: S''(x) = Mi + (Mi+1 - Mi)/hi+1 * (x - xi)
        S_vtoraya_proizv = M[i] + ((M[i+1] - M[i]) / h_segment) * (x_test - x_tochki[i])
        print(f"\nТочка проверки x = {x_test}")
        print(f"S''({x_test}) по сплайну = {S_vtoraya_proizv:.6f}")
        print(f"Точное f''({x_test}) = {f_vtoraya_proizvodnaya(x_test):.6f}")
        break