import numpy as np
import pandas as pd
import myproject as c_proj
import sys
from traceback import print_exc


def main():
    # starts reading the arguments here:
    n = len(sys.argv)
    if n < 5 or n > 6:
        print("Invalid Input!")
        sys.exit(1)
    max_iter = 300
    try:
        K = int(sys.argv[1])
        if n == 6:
            max_iter = int(sys.argv[2])
        epsilon = float(sys.argv[-3])
        filename_1 = sys.argv[-2]
        filename_2 = sys.argv[-1]
    except Exception:
        print("Invalid Input!")
        sys.exit(1)

    assert max_iter > 0, "Invalid Input!"
    np.set_printoptions(suppress=True)

    df, d, N = merge(filename_1, filename_2, K)

    df.drop([0], axis=1, inplace=True)
    df = df.to_numpy()

    indexes = kmeansPP(df, K, N)

    str_index = ""
    for i in range(K - 1):
        str_index += str(indexes[i]) + ","
    str_index += str(indexes[-1])
    try:
        a = c_proj.KMeans(K, N, d, max_iter, epsilon, df.flatten().tolist(), str_index)
        result = np.reshape(a, (-1, d)).tolist()
        print(str(indexes.tolist())[1:-1].replace(" ", ""))
        for row in result:
            print(str([round(x, 4) for x in row])[1:-1].replace(" ", ""))
    except Exception:
        print_exc()
        print("Invalid result!")


def merge(file1, file2, k):
    file_1 = open(file1, "r")
    text1 = [line.strip().split(",") for line in file_1]
    file_1.close()

    file_2 = open(file2, "r")
    text2 = [line.strip().split(",") for line in file_2]
    file_2.close()

    dimension1 = len(text1[0]) - 1
    dimension2 = len(text2[0]) - 1
    d = dimension1 + dimension2

    df1 = pd.DataFrame(data=text1)
    df2 = pd.DataFrame(data=text2)

    merged = df1.merge(df2, on=0, how='inner')
    N = merged.shape[0]

    assert N > k, "Invalid Input!"

    merged = merged.astype(float)
    merged = merged.sort_values(merged.columns[0], ascending=True)
    return merged, d, N


def kmeansPP(df, k, N):
    np.random.seed(0)
    first_cent = np.random.randint(N)
    centIndex = np.array([first_cent])

    # first run:
    minDist = np.sum(np.power((df - df[first_cent]), 2), axis=1)
    probability = np.divide(minDist, np.sum(minDist))
    centIndex = np.append(centIndex, np.random.choice(N, p=probability))

    # rest:
    for j in range(1, k - 1):
        dists = np.sum(np.power((df - df[centIndex[j]]), 2), axis=1)
        minDist = np.minimum(minDist, dists)
        probability = np.divide(minDist, np.sum(minDist))
        centIndex = np.append(centIndex, np.random.choice(N, p=probability))
    return centIndex


if __name__ == "__main__":
    main()