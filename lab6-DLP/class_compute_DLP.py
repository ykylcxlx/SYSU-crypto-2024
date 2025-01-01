from sympy import symbols, discrete_log
from concurrent.futures import ThreadPoolExecutor, as_completed

def modular_exponentiation(base, exponent, modulus):
    result = 1
    base = base % modulus

    while exponent > 0:
        if (exponent % 2) == 1:
            result = (result * base) % modulus
        exponent = exponent >> 1
        base = (base * base) % modulus
    
    return result

def compute_dlp(i):
    p = 3768901521908407201157691198029711972876087647970824596533
    n = 9993115456385501509
    alpha = modular_exponentiation(3107382411142271813235322646657672922264748410711464860476,
                                    2 * 2 * 23 * 8783 * 2419781956425763 * 192888768642311611 * i, p)
    beta = 2120553873612439845419858696451540936395844505496867133711

    x = discrete_log(p, beta, alpha)

    return i, x if x < n else None

stuid = [
    22336041, 22336090, 22336156, 22336167,
    22336198, 22336206, 22336232, 22336239,
    22336243, 22336275, 22336276, 22336280,
    22336304, 22336314, 22336327, 21306154,
    22302124, 21306375
]

# 使用 ThreadPoolExecutor 来并行计算
results = []
with ThreadPoolExecutor(max_workers=10) as executor:
    futures = {executor.submit(compute_dlp, i): i for i in stuid}
    for future in as_completed(futures):
        result = future.result()
        if result is not None:
            print(f"stuid: {result[0]}")
            print(f"DLP: {result[1]}")
        else:
            print(f"No solution found for stuid: {futures[future]}")