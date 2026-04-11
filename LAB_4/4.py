import math


def f(x):
    return x**2 + math.log(x) - 4


def pravye_pryamougolniki(a, b, n):
    shag = (b - a) / n
    summa = 0
    for i in range(1, n + 1):
        x = a + i * shag
        summa += f(x)
    return summa * shag


def vychislit_s_utochneniem(a, b, eps):
    n = 2
    i_n = pravye_pryamougolniki(a, b, n)

    print(f"{'n':>10} | {'I_n':>15} | {'Delta':>15}")
    print("-" * 60)

    while True:
        n_novoe = n * 2
        i_2n = pravye_pryamougolniki(a, b, n_novoe)
        delta = abs(i_2n - i_n)

        print(f"{n_novoe:>10} | {i_2n:>15.16f} | {delta:>15.16f}")

        if delta <= eps:
            return i_2n, n_novoe

        i_n = i_2n
        n = n_novoe


nizhnyaya_granica = 1.5
verhnyaya_granica = 2.0
tochnost = 0.001

itog, final_n = vychislit_s_utochneniem(nizhnyaya_granica, verhnyaya_granica, tochnost)

print("-" * 60)
print(f"Result: {itog:.8f}")
print(f"Final n: {final_n}")
