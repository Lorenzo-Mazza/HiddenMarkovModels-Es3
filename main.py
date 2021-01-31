import math


def matrix_fill(list_1d):
    row, column = int(list_1d[0]), int(list_1d[1])
    matrix = [[0.0 for x in range(column)] for y in range(row)]
    counter = 2
    for row_counter in range(row):
        for column_counter in range(column):
            matrix[row_counter][column_counter] = list_1d[counter]
            counter += 1
    return matrix

def list_fill(matrix):
    row, column = len(matrix),len(matrix[0])
    list_1d = [0.0 for x in range(row*column+2)]
    list_1d[0], list_1d[1]= row, column
    counter = 2
    for i in range(row):
        for j in range(column):
            list_1d[counter]= matrix[i][j]
            if list_1d[counter]< 1e-6:
                list_1d[counter]=0
            counter += 1
    return list_1d


def computeAlpha():
    for t in range(T):
        c[t] = 0
        for i in range(Nstates):
            if t==0:
                alpha[t][i]= b[i][obs[t + 1]] * pi[i + 2]
            else:
                sigma = 0
                for j in range(Nstates):
                    sigma += a[j][i] * alpha[t - 1][j]
                alpha[t][i]= b[i][obs[t + 1]] * sigma
            c[t] += alpha[t][i]
        c[t] = 1 / c[t]
        for i in range(Nstates):
            alpha[t][i] = c[t] * alpha[t][i]
            alpha[t][i] = round(alpha[t][i],8)

def computeBeta():
    for t in range(T)[::-1]:
        for i in range(Nstates):
            if t==T-1:
                beta[t][i] = 1
            else:
                beta[t][i] = 0
                for j in range(Nstates):
                    beta[t][i] += a[i][j] \
                                  * b[j][obs[t + 2]] * beta[t + 1][j]
            beta[t][i] = c[t]*beta[t][i]
            beta[t][i] = round(beta[t][i],8)



def computeGammas():
    for t in range(T-1):
        for i in range(Nstates):
            gamma[t][i] = 0
            for j in range(Nstates):
                digamma[t][i][j] = alpha[t][i] * a[i][j] * b[j][obs[t + 2]]* beta[t+1][j]
                gamma[t][i]+= digamma[t][i][j]
    for i in range(Nstates):
        gamma[T-1][i]= alpha[T-1][i]
        gamma[T-1][i] = round(gamma[T-1][i], 8)


def re_estimateParam():
    for i in range(Nstates):
        pi[i+2]= gamma[0][i]
        denomA, denomB = 0,0
        for t in range(T-1):
            denomA+= gamma[t][i]
        for j in range(Nstates):
            numerA= 0
            for t in range(T-1):
                numerA+= digamma[t][i][j]
            a[i][j]= numerA/denomA
        for t in range(T):
            denomB+= gamma[t][i]
        for j in range(M):
            numerB= 0
            for t in range(T):
                if obs[t+1]==j:
                    numerB+= gamma[t][i]
            b[i][j] = numerB/denomB


def computeLogProb():
    logProbability= 0
    for t in range(T):
        logProbability+= math.log(c[t])
    return -logProbability

def estimateModel():
    computeAlpha()
    computeBeta()
    computeGammas()
    re_estimateParam()


inp = input().strip().split(" ")
state_list = list(map(float, inp))
Nstates = int(state_list[0])
a = matrix_fill(state_list)
inp = input().strip().split(" ")
emissions_list = list(map(float, inp))
b = matrix_fill(emissions_list)
M = len(b[0])
inp = input().strip().split(" ")
pi = list(map(float, inp))
inp = input().strip().split(" ")
obs = list(map(int, inp))
T = obs[0]
alpha = [[0.0 for j in range(Nstates)] for i in range(T)]
beta = [[0.0 for j in range(Nstates)] for i in range(T)]
c = [0.0 for i in range(T)]
gamma = [[0.0 for j in range(Nstates)] for i in range(T)]
digamma = [[[0.0 for x in range(Nstates)] for y in range(Nstates)] for i in range(T)]
estimateModel()
logProb= computeLogProb()
oldLogProb = logProb-1
iterations = 1
MAX_ITER = 10
while (iterations < MAX_ITER) and (logProb > oldLogProb):
    threshold = logProb-oldLogProb
    if threshold<0.01:
        break
    oldLogProb = logProb
    estimateModel()
    logProb= computeLogProb()
    iterations += 1

if iterations== MAX_ITER:
    print("The algorithm doesn't converge")
else:
    a_list = list_fill(a)
    b_list = list_fill(b)
    print(' '.join(map(str, a_list)))
    print(' '.join(map(str, b_list)))


