import math
import numpy as np
W = np.array([[1, 1], [1, -1]],int)
p = (12, 8, 2, 1, 13, 4, 15, 6, 7, 0, 10, 5, 3, 14, 9, 11)
#####Треугольник Паскаля#####
def anf(function):
    buff = []
    iterations = len(function)
    ANF = [ function[0] ]
    for _ in range(iterations-1):
        for i in range(len(function)-1):
            buff.append(function[i]^function[i+1])
        print(buff)
        ANF.append(buff[0])
        function = buff
        buff = []
    return ANF
#############################
#####Матрица Адамара-Сильвестра#####
def hadamardMatrix(dimension):
    if dimension == 1:
        return W
    else:
        return np.hstack((np.vstack((hadamardMatrix(dimension-1),hadamardMatrix(dimension-1))),np.vstack((hadamardMatrix(dimension-1),-hadamardMatrix(dimension-1)))))
####################################
#####Получение массива значений функций#####
def functionValue(function):
    values = [ [] for _ in range(len(bin(max(function))[2:]))]
    for value in function:
        for function_number in range(len(bin(max(function))[2:])):
            values[function_number].append(value >> function_number & 1)
    return values
###################################################
#####Список Фурье образов для всех разрядных функций#####
def fourierSpectrum(function):
    values = functionValue(function)
    spectrum = []
    M = hadamardMatrix(int(math.log2(len(function))))
    for functions in range(len(bin(max(function))[2:])):
        spectrum.append(np.matmul(np.array(values[functions]).T,M).tolist())
    return spectrum
#########################################################
#####Список Уолш образов для всех разрядных функций#####
def walschSpectrum(function):
    values = functionValue(function)
    unit_array = np.array([1 for _ in range(len(function))])
    M = hadamardMatrix(int(math.log2(len(function))))
    spectrum = []
    for functions in range(len(bin(max(function))[2:])):
        spectrum.append((np.matmul((unit_array - 2 * np.array(values[functions])).T,M)).tolist())
    return [list(map(int, i)) for i in spectrum]
#########################################################
#####Подсчет вероятностей для строгого лавинного эффекта#####
def SAE(function):
        function = functionValue(function)
        result = []
        dim = (int)(math.log2(len(function[0])))
        num = 0
        for funcs in function:
                for var in range(len(funcs)): 
                        for i in range(dim):
                                if (funcs[var] != funcs[var ^ (1 << i)]):
                                        num += 1
                result.append(num/(2 ** dim * dim))
                num = 0
        return result
#############################################################
#####Подсчет вероятностей для критерия независимости битов#####
def BIC(function):
        p = (3, 5, 6, 9, 10, 12)
        result = [0.0,0.0,0.0,0.0,0.0,0.0]
        for i in range(16):
                for j in range(4):
                        for n,k in enumerate(p):
                                if ((function[i] ^ function[i^(1 << j)]) & k == k):
                                        result[n] += 1
        for i in range(6):
                result[i] = result[i]/64
        return result
###############################################################
#####Подсчет вероятности для лавинного эффекта#####
def AE(function):
        dim = (int)(math.log2(len(function)))
        num = []
        for var in range(2 ** dim): 
                for i in range(dim):
                        num.append(bin((function[var] ^ function[var ^ (1 << i)])).count('1'))
        return (sum(num)/(2 ** dim * dim * dim))
###################################################
