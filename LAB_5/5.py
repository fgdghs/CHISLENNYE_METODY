import numpy as np

# Я задаю размер системы
razmer = 5

# Я генерирую случайные элементы для диагоналей
A = np.random.uniform(1, 10, razmer)  # поддиагональ
A[0] = 0.0  # по условию (1) первый элемент поддиагонали A_1 = 0

C = np.random.uniform(1, 10, razmer)  # наддиагональ
C[razmer - 1] = 0.0  # по условию (1) последний элемент C_N = 0

# Я генерирую главную диагональ со строгим диагональным преобладанием
B = np.abs(A) + np.abs(C) + np.random.uniform(1, 10, razmer)

# Я создаю случайный вектор правой части
prav_chast = np.random.uniform(10, 50, razmer)


def metod_progonki(A, B, C, F):
    n = len(B)
    alfa = np.zeros(n)
    beta = np.zeros(n)

    # ПРЯМОЙ ХОД ПРОГОНКИ
    # Я инициализирую первые прогоночные коэффициенты (формула 2)
    alfa[0] = -C[0] / B[0]
    beta[0] = F[0] / B[0]

    # Я вычисляю остальные прогоночные коэффициенты (формулы 5)
    for i in range(1, n - 1):
        znamenatel = B[i] + A[i] * alfa[i - 1]
        alfa[i] = -C[i] / znamenatel
        beta[i] = (F[i] - A[i] * beta[i - 1]) / znamenatel

    # ОБРАТНЫЙ ХОД ПРОГОНКИ
    reshenie_U = np.zeros(n)

    # Я нахожу значение U_N для последнего уравнения
    znamenatel_n = B[n - 1] + A[n - 1] * alfa[n - 2]
    reshenie_U[n - 1] = (F[n - 1] - A[n - 1] * beta[n - 2]) / znamenatel_n

    # Я последовательно нахожу остальные значения с конца (формула 6)
    for i in range(n - 2, -1, -1):
        reshenie_U[i] = alfa[i] * reshenie_U[i + 1] + beta[i]

    return reshenie_U


# Я получаю искомое решение системы
reshenie = metod_progonki(A, B, C, prav_chast)

# ПРОВЕРКА РЕЗУЛЬТАТА
# Я собираю полную плотную матрицу из диагоналей, чтобы умножить ее на вектор решения
matrica_A = np.diag(B) + np.diag(A[1:], k=-1) + np.diag(C[:-1], k=1)

# Я вычисляю вектор невязки и его норму |A*x - f|
raznost = np.dot(matrica_A, reshenie) - prav_chast
nevyazka = np.linalg.norm(raznost)

print(
    "Сгенерированная матрица системы (с диагональным преобладанием):\n",
    np.round(matrica_A, 2),
)
print("\nВектор правой части F:\n", np.round(prav_chast, 2))
print("\nНайденное решение U:\n", np.round(reshenie, 4))
print("\nПроверка (норма вектора невязки |A*x - f|):", nevyazka)
