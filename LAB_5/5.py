import numpy as np

razmer = 5

A = np.random.uniform(1, 10, razmer)  # поддиагональ
A[0] = 0.0  # по условию  первый элемент поддиагонали A_1 = 0

C = np.random.uniform(1, 10, razmer)  # наддиагональ
C[razmer - 1] = 0.0  # по условию  последний элемент C_N = 0

# Главная диагональ с диагональным преобладанием (Элемент на диагонали строго больше суммы своих соседей по строке)
B = np.abs(A) + np.abs(C) + np.random.uniform(1, 10, razmer)

# случайный вектор правой части
prav_chast = np.random.uniform(10, 50, razmer)


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

    # Остальные значения последовательно
    for i in range(n - 2, -1, -1):
        reshenie_U[i] = alfa[i] * reshenie_U[i + 1] + beta[i]

    return reshenie_U


reshenie = metod_progonki(A, B, C, prav_chast)

matrica_A = (
    np.diag(B) + np.diag(A[1:], k=-1) + np.diag(C[:-1], k=1)
)  # k - сколько отступать от главной диагонали

raznost = np.dot(matrica_A, reshenie) - prav_chast
pogreshost = np.linalg.norm(raznost)

print(
    "Сгенерированная матрица системы (с диагональным преобладанием):\n",
    np.round(matrica_A, 2),
)
print("\nВектор правой части F:\n", np.round(prav_chast, 2))
print("\nНайденное решение U:\n", np.round(reshenie, 4))
print("\nПроверка |A*x - f|:", pogreshost)
