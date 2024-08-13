import copy
import numpy as np

q = 2
directions = (0, 3)

# Convenience function
def dimension(simplex):
    return len(simplex) - 1
def main():
    global q
    global directions
    i =directions[0]
    j = directions[1]
    # Example A.2.1
    B = ([
             [1 , 1 , 1 , 0 , 0],
             [0 , 1 , 1 , 1 , 0],
             [1 , 0 , 1 , 1 , 0],
             [1 , 1 , 0 , 1 , 0],
             [0 , 1 , 1 , 0 , 1],
             [0 , 1 , 0 , 1 , 1],
             [0 , 0 , 1 , 1 , 1],
             [1 , 1 , 1 , 1 , 0],
             [0 , 1 , 1 , 1 , 1]
    ])
    B = np.array(B)
    print("Q_hat:")
    print(compute_Q_original(B, q, i, j))
    print("Q:")
    print(compute_Q_new(B,q,i,j))

def compute_Q_original(B, q, i, j):
    number_simplices = B.shape[0]
    inclusion = argmax(B @ B.T) - np.diag([1]*number_simplices)
    d_iB = face(B, i)
    d_jB = face(B, j)
    q_near = threshold(d_iB @ d_jB.T, q)
    Q = threshold(q_near + inclusion, 0)
    for x in range(number_simplices):
        Q[x,x] = 1
    return Q

def compute_Q_new(B, q, i, j):
    number_simplices = B.shape[0]
    inclusion = argmax(B @ B.T) - np.diag([1] * number_simplices)

    B_Q1 = B[np.sum(B, axis=1) == q+2]

    C_S = threshold(B_Q1 @ B.T, q+1)

    q_plus_one_near = threshold(face(B_Q1, i) @ face(B_Q1, j).T, q) #+ np.diag([1] * len(simplices[q+1]))
    for x in range(B_Q1.shape[0]):
        q_plus_one_near[x,x] = 1
    q_near_new_ = threshold(C_S.T @ q_plus_one_near @ C_S, 0)
    for x in range(B.shape[0]):
        q_near_new_[x,x] = 1
    Q = threshold(inclusion + q_near_new_, 0)
    return Q

def threshold(M, p):
    def greater_p(x):
        return 1 if x > p else 0
    return np.vectorize(greater_p)(M)
def face(B, i):
    result = copy.deepcopy(B)
    for r in range(B.shape[0]):
        index = 0
        overflow = True
        for c in range(B.shape[1]):
            if result[r][c] != 0:
                if index == i:
                    overflow = False
                    result[r][c] = 0
                index += 1
        if overflow:
            for c in range(B.shape[1]-1, 0, -1):
                if result[r][c] != 0:
                    result[r][c] = 0
                    break
    return result


def argmax(M):
    result = np.zeros_like(M)
    for row in range(0, M.shape[0]):
        max = np.max(M[row, :])
        for index, e in enumerate(M[row, :]):
            if e == max:
                result[row, index] = 1
    return result





if __name__ == '__main__':
    main()
