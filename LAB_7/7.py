import numpy as np

def f(x):
    return 2 * np.log(x) - 1/x

def df(x):
    return 2/x + 1/x**2

def d2f(x):
    return -2/x**2 - 2/x**3

# Поиск интервалов где находится корень
def find_intervals(start, end, step):
    intervals = []
    curr = start
    while curr < end:
        nxt = curr + step
        # Теорема Больцано-Коши
        if f(curr) * f(nxt) < 0:
            intervals.append((curr, nxt))
        curr = nxt
    return intervals

def solve_combined_method(a, b, eps):
    if f(a) * d2f(a) > 0:
        xn = a           
        xn_chord = b    
    else:
        xn = b
        xn_chord = a

    print(f"\nУточнение корня на отрезке [{a}, {b}]")
    print(f"{'Итерация':<7} | {'x_кас (xn)':<18} | {'x_хорд (xn_chord)':<18} | {'Разность':<12}")
    print("-" * 70)

    iteration = 0
    while abs(xn - xn_chord) > eps:
        iteration += 1
        
        f_xn = f(xn)
        f_ch = f(xn_chord)
            
        xn = xn - f_xn / df(xn)

        xn_chord = xn_chord - (f_ch * (xn - xn_chord)) / (f(xn) - f_ch)
        
        print(f"{iteration:<7} | {xn:<18.12f} | {xn_chord:<18.12f} | {abs(xn - xn_chord):<12.2e}")

    return (xn + xn_chord) / 2

intervals = find_intervals(0.1, 5, 0.1)

print(f"Найдено интервалов : {len(intervals)}")

for interval in intervals:
    a, b = interval
    res = solve_combined_method(a, b, eps=1e-15)
        
    print("-" * 70)
    print(f"Итоговый корень: {res:.14f}")
    print(f"Значение f(x):   {f(res):.2e}")