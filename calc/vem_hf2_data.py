# Data

code_2 = 2
code_3 = 2
code_4 = 2

# Parameters

a_list = [1.2, 1.7, 2.1, 2.6]
m_0_list = [15, 20, 25, 30]
b_list = [5, 6, 7, 8]
d_list = [25, 35, 45, 55]
E_list = [170, 190, 210, 230]
ρ_list = [6000, 6500, 7000, 7500]

def get_data(code2, code3, code4):
    data = {
        'a': a_list[code2 - 1],
        'm_0': m_0_list[code2 - 1],
        'b': b_list[code3 - 1],
        'd_mm': d_list[code3 - 1],
        'E_GPa': E_list[code4 - 1],
        'ρ': ρ_list[code4 - 1]
    }
    # Change to SI

    data['E'] = data['E_GPa']*10**9
    data['d'] = data['d_mm']/1000
    return data
