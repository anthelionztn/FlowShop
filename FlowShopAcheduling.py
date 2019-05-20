# -*- coding: utf-8 -*-
import copy
import math
import random
import matplotlib.pyplot as plt
import numpy as np

pack_size = 200
chrom_length = 8  # 相当于产品的个数
prob_mutate = 0.05  # 变异率
fit_value = []  # 适应性分数
epoch_num = 100  # 进化世代数
mean_value = [] # 记录每一代个体耗时的均值
min_value = [] # 记录每一代个体耗时的最小值

# 用户填写产品工艺耗时矩阵,这里用8种产品，每种产品有7个工艺环节为例
processTimeMatrix = [[10, 5, 7, 9, 25, 7, 8],
                     [12, 6, 3, 23, 11, 5, 9],
                     [8, 8, 5, 15, 12, 7, 10],
                     [4, 20, 6, 11, 14, 7, 13],
                     [9, 5, 12, 14, 5, 9, 10],
                     [9, 7, 16, 14, 5, 7, 10],
                     [9, 7, 16, 14, 5, 7, 10],
                     [9, 7, 16, 14, 5, 7, 10]]


# 随机生成种群，使用整数编码
def geneEncoding(pack_size, chrom_length):
    pack = []
    lst = []
    for i in range(chrom_length):
        lst.append(i)

    for j in range(pack_size):
        pack.append(random.sample(lst, len(lst)))
    return pack


def decoding(weightList, processTimeMatrix):
    newprocessTimeMatrix = []
    for i in range(len(processTimeMatrix)):
        newprocessTimeMatrix.append(processTimeMatrix[weightList[i]])
    return newprocessTimeMatrix


def calobjValue(pack):
    maxTimeValue = []
    temp = []
    for k in range(pack_size):
        temp = decoding(pack[k], processTimeMatrix)
        timeCostMatrix = [[0 for i in range(len(temp[0]))] for i in range(len(temp))]
        for i in range(len(temp)):
            for j in range(len(temp[i])):
                if i == 0 and j == 0:
                    timeCostMatrix[i][j] = temp[0][0]
                elif i == 0 and j != 0:
                    timeCostMatrix[i][j] = timeCostMatrix[i][j - 1] + temp[i][j]
                elif i != 0 and j == 0:
                    timeCostMatrix[i][j] = timeCostMatrix[i - 1][j] + temp[i][j]
                else:
                    timeCostMatrix[i][j] = max(timeCostMatrix[i][j - 1], timeCostMatrix[i - 1][j]) + temp[i][j]
        maxTimeValue.append(timeCostMatrix[-1][-1])
    #print maxTimeValue
    return maxTimeValue


def calfitValue(maxTimeValue):
    temp = 0
    lst = []
    for i in range(len(maxTimeValue)):
        lst.append(1 / float(maxTimeValue[i]))
    fit_value = map(costFunction, lst)
    sumValue = sum(fit_value)
    for j in range(len(fit_value)):
        fit_value[j] /= sumValue
        temp += fit_value[j]
        fit_value[j] = temp
    return fit_value


def costFunction(x):
    return math.pow(x, 9)


# 择优复制，适应性分数越高的个体被复制的几率越高，返回被复制个体在种群中的索引
def selection(cumValue):
    newPackIndex = []
    for i in range(pack_size):
        temp = random.random()
        j = 0
        while cumValue[j] < temp:
            j += 1
        newPackIndex.append(j)
    return newPackIndex


# 按照被复制个体在种群中的索引复制生成新一代种群
def geneNewPack(oldPack, newPackIndex):
    # print oldPack #为了观察染色体变异情况
    newPack = []
    for i in range(pack_size):
        newPack.append(mutate(oldPack[newPackIndex[i]]))
    return newPack


# 染色体按一定概率随机变异
def mutate(oldGene):
    temp = random.random()
    newGene = copy.copy(oldGene)
    if temp >= prob_mutate:
        return oldGene
    else:
        inverseBit = random.randrange(0, chrom_length - 2, 1)
        # print 'inverseBit',inverseBit
        bitBox = newGene[inverseBit]
        newGene[inverseBit] = newGene[inverseBit + 1]
        newGene[inverseBit + 1] = bitBox
        return newGene


# 杂交两个个体
def exchange(a, b):
    exchangePoint1 = 1
    exchangePoint2 = 1
    while exchangePoint1 == exchangePoint2:
        exchangePoint1 = random.randrange(0, chrom_length / 2, 1)
        exchangePoint2 = random.randrange(chrom_length / 2, chrom_length - 1, 1)
    aNewL = a[:exchangePoint1 + 1]
    aNewR = a[exchangePoint2 + 1:]
    aNewM = []
    for i in range(chrom_length):
        if b[i] not in (aNewL + aNewR):
            aNewM.append(b[i])
    aNew = aNewL + aNewM + aNewR
    bNewL = a[exchangePoint1 + 1:exchangePoint2 + 1]
    bNewR = []
    for j in range(chrom_length):
        if b[j] not in aNewM:
            bNewR.append(b[j])
    bNew = bNewL + bNewR
    return aNew, bNew


# 杂交过程
def crossover(newPack):
    for i in range(int(pack_size / 2)):
        tpl = exchange(newPack[i], newPack[pack_size - i - 1])

        newPack[i] = mutate(tpl[0])
        newPack[pack_size - i - 1] = mutate(tpl[1])
    return newPack


# 世代循环
def epoch(pack, T):
    for i in range(T):
        res = calobjValue(pack)
        #print 'Mean value of Epoch', i, ':', np.mean(res)
        mean_value.append(np.mean(res))
        min_value.append(np.min(res))
        fit = calfitValue(res)
        sel = selection(fit)
        newPack = geneNewPack(pack, sel)
        pack = crossover(newPack)
    return pack

first_pack = geneEncoding(pack_size, chrom_length)  # 随机生成第一世代种群
last_pack = epoch(first_pack, epoch_num)
last_pack_time = calobjValue(last_pack)
bestValue = min(last_pack_time)
print 'the optimal time:', bestValue
print 'the optimal solution:', last_pack[last_pack_time.index(bestValue)]  # 打印结果
plt.plot(mean_value)
plt.show()