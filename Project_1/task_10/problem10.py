import numpy as np

# Time-arrays for the general algorithm
time_n10 = np.array([0, 0, 0, 0, 0, 0])
time_n100 = np.array([5.5e-5, 0.000334, 0, 7.4e-5, 0, 0])
time_n1000 = np.array([0.001002, 0.001005, 0.000503, 0.000783, 0.000502, 0.000501])
time_n10000 = np.array([0.005946, 0.005819, 0.007009, 0.0055, 0.0045, 0.004502])
time_n100000 = np.array([0.058995, 0.063889, 0.056993, 0.048521, 0.056489, 0.059116])
time_n1000000 = np.array([0.606497, 0.586001, 0.577265, 0.503477, 0.587996, 0.586002])

time_array_general = np.array([time_n10, time_n100, time_n1000, time_n10000, time_n100000, time_n1000000])

print('\n')
print('For the general algorithm:')
print('\t\t\tmean\t\tstd')
for n, time_n in enumerate(time_array_general, start=1):
    time_n_mean = np.mean(time_n)
    time_n_std = np.std(time_n)
    print(f'Time for n={10**n:.0e}:\t{time_n_mean:.3e} +- {time_n_std:.3e} s')

print('')

# Time-arrays for the special algorithm
time_n10 = np.array([0, 0, 0, 0, 0, 0])
time_n100 = np.array([0, 0, 0, 0, 0, 0])
time_n1000 = np.array([0.000499, 0.000343, 0.000502, 0.000611, 0.000666, 0.000581])
time_n10000 = np.array([0.00304, 0.003, 0.002998, 0.003498, 0.003725, 0.002389])
time_n100000 = np.array([0.032, 0.032, 0.026998, 0.031727, 0.031999, 0.027001])
time_n1000000 = np.array([0.319499, 0.3195, 0.271499, 0.318999, 0.3185, 0.271])

time_array_special = np.array([time_n10, time_n100, time_n1000, time_n10000, time_n100000, time_n1000000])

print('For the special algorithm:')
print('\t\t\tmean\t\tstd')
for n, time_n in enumerate(time_array_special, start=1):
    time_n_mean = np.mean(time_n)
    time_n_std = np.std(time_n)
    print(f'Time for n={10**n:.0e}:\t{time_n_mean:.3e} +- {time_n_std:.3e} s')
    
print('')

print('Ratio Special / General:')
print('\t\t\tmean')
for n, (special, general) in enumerate(zip(time_array_special[1:], time_array_general[1:]), start=2):
    special_mean = np.mean(special)
    general_mean = np.mean(general)
    print(f'Relation for n={10**n:.0e}:\t{special_mean / general_mean}')
