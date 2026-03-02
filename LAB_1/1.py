import math


# --- ШАГ 1: ОПРЕДЕЛЕНИЕ ФУНКЦИИ И ПРОИЗВОДНЫХ ---
def f(x):
    return x**2 + math.log(x) - 4


# Вторая производная f''(x) для оценки погрешности R1
def f_vtoraya_proizv(x):
    return 2 - (1 / x**2)


# Данные твоего варианта (№3)
nachalo, konec = 1.5, 2.0
kolichestvo_uzlov = 10
shag = (konec - nachalo) / kolichestvo_uzlov

# Точка для расчета (x* и x** одинаковы в табл. 1 для 3-го варианта)
x_celevoy = 1.52

# --- ШАГ 2: ПОСТРОЕНИЕ ТАБЛИЦЫ ЗНАЧЕНИЙ (Пункт 1 задания) ---
uzly_x = [round(nachalo + i * shag, 2) for i in range(kolichestvo_uzlov + 1)]
uzly_y = [f(x) for x in uzly_x]

print(f"--- ТАБЛИЦА ЗНАЧЕНИЙ ФУНКЦИИ ---")
for i, (ux, uy) in enumerate(zip(uzly_x, uzly_y)):
    print(f"x{i}: {ux:.2f} | f(x{i}): {uy:.6f}")

# --- ШАГ 3: АВТОМАТИЧЕСКИЙ ВЫБОР РАБОЧЕГО ИНТЕРВАЛА ---
# Вместо того чтобы писать i=0 руками, мы просим программу найти
# индекс j такой, чтобы наш x_celevoy лежал между узлами x[j] и x[j+1].
indeks = 0
for j in range(len(uzly_x) - 1):
    if uzly_x[j] <= x_celevoy <= uzly_x[j + 1]:
        indeks = j
        break

# Теперь присваиваем значения найденных узлов
xi, xi1 = uzly_x[indeks], uzly_x[indeks + 1]
yi, yi1 = uzly_y[indeks], uzly_y[indeks + 1]

print(f"\nВыбран рабочий интервал: [{xi}, {xi1}] (индекс i={indeks})")

# --- ШАГ 4: РАСЧЕТ ПО ФОРМУЛЕ ЛАГРАНЖА (L1) ---
# Линейная интерполяция через весовые коэффициенты узлов.
lagranzh1 = yi * (x_celevoy - xi1) / (xi - xi1) + yi1 * (x_celevoy - xi) / (xi1 - xi)

# --- ШАГ 5: РАСЧЕТ ПО ФОРМУЛЕ НЬЮТОНА (N1) ---
# Линейная интерполяция через разделенные разности.
razdel_raznost = (yi1 - yi) / (xi1 - xi)
nyuton1 = yi + razdel_raznost * (x_celevoy - xi)

# --- ШАГ 6: ОЦЕНКА ПОГРЕШНОСТИ (Пункт 2 задания) ---
y_tochnoe = f(x_celevoy)
pogreshnost_fakt = lagranzh1 - y_tochnoe

# Оцениваем R1 через вторую производную на краях интервала [xi, xi1]
f2_znacheniya = [f_vtoraya_proizv(xi), f_vtoraya_proizv(xi1)]
omega2 = (x_celevoy - xi) * (x_celevoy - xi1)

# Вычисляем границы теоретического остаточного члена
r1_granica1 = (f2_znacheniya[0] * omega2) / 2
r1_granica2 = (f2_znacheniya[1] * omega2) / 2

# --- ВЫВОД РЕЗУЛЬТАТОВ ---
print(f"\n--- РЕЗУЛЬТАТЫ ДЛЯ x = {x_celevoy} ---")
print(f"Интерполяция Лагранжа L1: {lagranzh1:.8f}")
print(f"Интерполяция Ньютона N1:   {nyuton1:.8f}")
print(f"Точное значение f(x):      {y_tochnoe:.8f}")
print(f"Фактическая погрешность:   {abs(pogreshnost_fakt):.8e}")

print(f"\n--- ОЦЕНКА ПОГРЕШНОСТИ ---")
teor_min, teor_max = min(r1_granica1, r1_granica2), max(r1_granica1, r1_granica2)
print(f"Теоретический интервал R1: ({teor_min:.8f}, {teor_max:.8f})")
print(
    f"R1 фактическое ({pogreshnost_fakt:.8f}) попадает в интервал? "
    f"{'Да' if teor_min < pogreshnost_fakt < teor_max else 'Нет'}"
)

# Проверка на соответствие точности 10^-4
dopustimo = abs(pogreshnost_fakt) <= 1e-4
print(
    f"\nДопустима ли линейная интерполяция (ошибка <= 10^-4)? {'Да' if dopustimo else 'Нет'}"
)
