import numpy as np

def f(x):
    return 2 * np.log(x) - 1/x

def df(x):
    return 2/x + 1/x**2

def d2f(x):
    return -2/x**2 - 2/x**3

def solve_combined_method(a, b, eps=1e-5):
    # Проверка условия сходимости для метода касательных (f(x) * f''(x) > 0)
    # На [1, 2] f''(x) < 0. Значит выбираем конец, где f(x) < 0
    if f(a) * d2f(a) > 0:
        xn = a # приближение касательными
        xn_chord = b # приближение хордами
    else:
        xn = b
        xn_chord = a

    print(f"{'Итерация':<10} | {'x_касат (xn)':<15} | {'x_хорд (xn_wave)':<15} | {'Разность':<15}")
    print("-" * 65)

    iteration = 0
    while abs(xn - xn_chord) > eps:
        # Метод касательных (Ньютона)
        xn = xn - f(xn) / df(xn)
        
        # Метод хорд
        xn_chord = xn_chord - (f(xn_chord) * (xn - xn_chord)) / (f(xn) - f(xn_chord))
        
        iteration += 1
        print(f"{iteration:<10} | {xn:<15.6f} | {xn_chord:<15.6f} | {abs(xn - xn_chord):<15.6e}")

    return (xn + xn_chord) / 2

# Начальный отрезок [1, 2]
root = solve_combined_method(1, 2)
print("-" * 65)
print(f"Найденный корень: {root:.6f}")
print(f"Проверка f(root): {f(root):.2e}")