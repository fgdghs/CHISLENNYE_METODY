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
y_znacheniya = f(x_tochki)
N = len(x_tochki) - 1

# Я вычисляю шаги сетки h_i
h = np.diff(x_tochki)

# Краевые условия 2-го типа по методичке: M0 = A, MN = B
M_0 = f_vtoraya_proizvodnaya(x_tochki[0])
M_N = f_vtoraya_proizvodnaya(x_tochki[-1])

# Я инициализирую массивы для коэффициентов СЛАУ (размер N-1)
n_vnutr = N - 1
koeff_A = np.zeros(n_vnutr) # поддиагональ
koeff_B = np.zeros(n_vnutr) # главная диагональ
koeff_C = np.zeros(n_vnutr) # наддиагональ
pravaya_chast_F = np.zeros(n_vnutr)

# Я заполняю коэффициенты системы строго по формулам из методички
for i in range(1, N):
    idx = i - 1
    h_i = h[i-1]
    h_i_sled = h[i]
    
    # Коэффициенты лямбда и мю
    lyambda_i = h_i / (h_i + h_i_sled)
    myu_i = h_i_sled / (h_i + h_i_sled)
    
    # Главная диагональ всегда 2
    koeff_B[idx] = 2.0
    
    # Заполнение поддиагонали и наддиагонали
    if idx > 0:
        koeff_A[idx] = lyambda_i
    if idx < n_vnutr - 1:
        koeff_C[idx] = myu_i
        
    # Правая часть d_i
    d_i = (6.0 / (h_i + h_i_sled)) * (
        (y_znacheniya[i+1] - y_znacheniya[i]) / h_i_sled - 
        (y_znacheniya[i] - y_znacheniya[i-1]) / h_i
    )
    
    # Учет краевых условий 2-го типа в правой части
    if idx == 0:
        d_i -= lyambda_i * M_0
    if idx == n_vnutr - 1:
        d_i -= myu_i * M_N
        
    pravaya_chast_F[idx] = d_i

# Я нахожу внутренние моменты, используя ВАШ метод прогонки
M_vnutrennie = metod_progonki(koeff_A, koeff_B, koeff_C, pravaya_chast_F)

# Я собираю итоговый массив всех моментов M_i
M = np.zeros(N + 1)
M[0] = M_0
M[N] = M_N
M[1:N] = M_vnutrennie

# Вывод результатов
print("Моменты M_i (вторые производные в узлах):")
for i, val in enumerate(M):
    print(f"M_{i} = {val:.6f}")

# Я рассчитываю значение S''(x) в контрольной точке x**** = 1.77
x_test = 1.77
for i in range(N):
    if x_tochki[i] <= x_test <= x_tochki[i+1]:
        h_segment = h[i]
        # Формула расчета второй производной сплайна
        S_vtoraya_proizv = M[i] + ((M[i+1] - M[i]) / h_segment) * (x_test - x_tochki[i])
        print(f"\nТочка проверки x = {x_test}")
        print(f"S''({x_test}) по сплайну = {S_vtoraya_proizv:.6f}")
        print(f"Точное f''({x_test}) = {f_vtoraya_proizvodnaya(x_test):.6f}")
        break