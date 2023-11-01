# -*- coding:UTF-8 -*-
import sys
import time
import math
import collections
from binary_heap import *
import random
from random import choice

class Graph:

    def __init__(self, dataset):
        self.dataset = dataset
        self.adj_list, self.deg = self.graph_list()
        self.n = len(self.adj_list)
        self.m = 0
        for u in self.adj_list:
            self.m += len(self.adj_list[u])
        self.m = self.m / 2
        print("number of nodes:" + str(self.n))
        print("number of edges:" + str(self.m))

    def graph_list(self):
        starttime = time.time()
        adj_list = {}
        print("........................")
        print(self.dataset + " is loading...")
        file = open(self.dataset)
        while 1:
            lines = file.readlines(50000)
            if not lines:
                break
            for line in lines:
                line = line.split()
                if str(line[0]) == "#":
                    continue
                if int(line[0]) == int(line[1]):
                    continue
                doubleedge = [(int(line[0]), int(line[1])), (int(line[1]), int(line[0]))]
                for (f_id, t_id) in doubleedge:
                    if f_id in adj_list:
                        adj_list[f_id].add(t_id)
                    else:
                        adj_list[f_id] = {t_id}
        file.close()
        deg = {}
        for u in adj_list:
            deg[u] = len(adj_list[u])
        endtime = time.time()
        print("load_time(s)" + str(endtime - starttime))
        return adj_list, deg

    # Compute conductance
    def compu_conductance(self, S):
        d = {}
        dout = {}
        sum_d = 0
        sum_dout = 0
        for u in S:
            d[u] = len(self.adj_list[u])
            sum_d += d[u]
            dout[u] = d[u] - len(set(self.adj_list[u]) & set(S))
            sum_dout += dout[u]
        if 2 * self.m != sum_d:
            condu = sum_dout / min(sum_d, 2 * self.m - sum_d)
            return condu
        else:
            return 0  # 最好的情况condu = 0.0

    # ISCCI
    def local_greedy_dr(self, T, l, h):
        starttime = time.time()
        greedy_best_set = set()
        greedy_best_condu = float('inf')
        f = {}
        dr = {}
        g = dict.fromkeys(self.adj_list.keys(), 0)
        for i in range(T):
            C = set()
            count = 0
            mymin_heap = MinHeap([])
            mymin_heap.insert([QID, 1])
            g_s_upper = 0
            g_s_lower = 0
            volC = 0
            best_g = 0
            best_index = 0
            D = []
            while volC < self.m and count < h and len(mymin_heap.heap) != 0:
                u = mymin_heap.heap[0][0]
                g[u] = mymin_heap.heap_dict[u]
                mymin_heap.remove()

                C.add(u)
                D.append(u)
                count += 1
                volC += self.deg[u]
                g_s_upper += 2 * len(set(self.adj_list[u]) & C)
                g_s_lower += self.deg[u]
                if g_s_upper / g_s_lower >= best_g and volC <= self.m and l <= count <= h:
                    best_index = count
                    best_g = g_s_upper / g_s_lower

                for v in self.adj_list[u]:
                    if v not in C:
                        dr[v] = self.deg[v] / len(set(self.adj_list[v]) & set(C))
                        f[v] = dr[v] + g[v]
                        idxs = [i[0] for i in mymin_heap.heap]
                        if v not in idxs:
                            mymin_heap.insert([v, f[v]])
                        elif v in idxs:
                            if f[v] >= mymin_heap.heap_dict[v]:
                                mymin_heap.increase_key(v, f[v])
                            else:
                                mymin_heap.decrease_key(v, f[v])

            greedy_set = set()
            for x in range(best_index):
                greedy_set.add(D[x])
            greedy_condu = self.compu_conductance(greedy_set)

            if greedy_condu <= greedy_best_condu:
                greedy_best_condu = greedy_condu
                greedy_best_set = greedy_set.copy()

        endtime = time.time()
        greedy_time = endtime - starttime
        return greedy_time, greedy_best_set, greedy_best_condu

    def local_greedy_(self, l, h):
        starttime = time.time()
        C = []
        volC = 0
        count = 0
        cr = {}
        mymax_heap = MaxHeap([])
        mymax_heap.insert([QID, 1])
        a = 0
        b = 0
        greedy_set = set()
        greedy_condu = float('inf')
        while volC < self.m and count < h and len(mymax_heap.heap) != 0:
            u = mymax_heap.remove()[0]

            C.append(u)
            count += 1
            volC += self.deg[u]
            a += self.deg[u] - 2 * len(set(self.adj_list[u]) & set(C))  # a
            b += self.deg[u]  # b
            if volC <= self.m and l <= count <= h:
                temp_condu = self.compu_conductance(C)
                if temp_condu <= greedy_condu:
                    greedy_set = C.copy()
                    greedy_condu = temp_condu
            for u in C:
                for v in self.adj_list[u]:
                    if v not in C:
                        d1 = self.deg[v] - 2 * len(set(self.adj_list[v]) & set(C))
                        d2 = self.deg[v]
                        cr[v] = (a * d2 - b * d1) / (b * (b + d2))
                        idxs = [i[0] for i in mymax_heap.heap]
                        if v not in idxs:
                            mymax_heap.insert([v, cr[v]])
                        elif v in idxs:
                            if cr[v] >= mymax_heap.heap_dict[v]:
                                mymax_heap.increase_key(v, cr[v])
                            else:
                                mymax_heap.decrease_key(v, cr[v])

        endtime = time.time()
        greedy_time = endtime - starttime
        return greedy_time, greedy_set, greedy_condu

    # ISCCP
    def local_greedy_cr(self, T, l, h):
        starttime = time.time()
        time1, greedy_best_set, greedy_best_condu = self.local_greedy_(l, h)
        for i in range(T - 1):
            C = []
            volC = 0
            count = 0
            a = 0
            b = 0
            cr = {}

            w = random.choice(greedy_best_set[1:])
            for u in greedy_best_set:
                if u == w:
                    break
                volC += self.deg[u]
                a += self.deg[u] - 2 * len(set(self.adj_list[u]) & set(C))  # a
                b += self.deg[u]  # b
                count += 1
                C.append(u)

            mymax_heap = MaxHeap([])
            for u in C:
                for v in self.adj_list[u]:
                    if v not in C and v != w:
                        d1 = self.deg[v] - 2 * len(set(self.adj_list[v]) & set(C))
                        d2 = self.deg[v]
                        cr[v] = (a * d2 - b * d1) / (b * (b + d2))
                        idxs = [i[0] for i in mymax_heap.heap]
                        if v not in idxs:
                            mymax_heap.insert([v, cr[v]])

            while volC < self.m and count < h and len(mymax_heap.heap) != 0:
                u = mymax_heap.remove()[0]

                C.append(u)
                count += 1
                volC += self.deg[u]
                a += self.deg[u] - 2 * len(set(self.adj_list[u]) & set(C))  # a
                b += self.deg[u]  # b
                if volC <= self.m and l <= count <= h:
                    temp_condu = self.compu_conductance(C)
                    if temp_condu <= greedy_best_condu:
                        greedy_best_set = C.copy()
                        greedy_best_condu = temp_condu

                for v in self.adj_list[u]:
                    if v not in C:
                        d1 = self.deg[v] - 2 * len(set(self.adj_list[v]) & set(C))
                        d2 = self.deg[v]
                        cr[v] = (a * d2 - b * d1) / (b * (b + d2))
                        idxs = [i[0] for i in mymax_heap.heap]
                        if v not in idxs:
                            mymax_heap.insert([v, cr[v]])
                        elif v in idxs:
                            if cr[v] >= mymax_heap.heap_dict[v]:
                                mymax_heap.increase_key(v, cr[v])
                            else:
                                mymax_heap.decrease_key(v, cr[v])

        greedy_best_set = set(greedy_best_set)
        endtime = time.time()
        greedy_time = endtime - starttime
        return greedy_time, greedy_best_set, greedy_best_condu

    def sweep_cut(self, pi, l, h):
        starttime = time.time()
        for u in pi:
            pi[u] = pi[u] / self.deg[u]
        pi = sorted(pi.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)
        S = set()
        volS = 0
        cutS = 0
        best_condu, best_index, count = float('inf'), 0, 0
        best_set = set()
        for x in pi:
            u = x[0]
            S.add(u)
            count += 1
            cutS = cutS - 2 * len(set(self.adj_list[u]) & S) + self.deg[u]
            volS = volS + self.deg[u]
            if min(volS, 2 * self.m - volS) != 0 and cutS / min(volS,
                                                                2 * self.m - volS) < best_condu and l <= count <= h:
                best_condu = cutS / min(volS, 2 * self.m - volS)
                best_index = count
        for x in range(best_index):
            best_set.add(pi[x][0])
        if len(best_set) > self.n / 2:
            best_set = set(self.adj_list) - set(best_set)
        endtime = time.time()
        sweep_time = endtime - starttime
        return sweep_time, best_set, best_condu

    def fpush(self, alpha, epsilon):  # NIBBLE_PPR
        starttime = time.time()
        pi, r = {}, {}
        r[QID] = 1
        q = collections.deque()
        q.append(QID)
        while q:
            u = q.popleft()
            for v in self.adj_list[u]:
                if v not in r:
                    r[v] = 0
                update = (1 - alpha) * r[u] / self.deg[u]
                r[v] = r[v] + update  # unweighted graph
                if (r[v] - update) / self.deg[v] < epsilon and r[v] / self.deg[v] >= epsilon:
                    q.append(v)
            if u not in pi:
                pi[u] = 0
            pi[u] = pi[u] + alpha * r[u]
            r[u] = 0
        endtime = time.time()
        return pi, endtime - starttime

    def hk_relax(self, t, epsilon):
        N = int(2 * t * math.log(1 / epsilon)) + 1
        starttime = time.time()
        psis = {}
        psis[N] = 1.
        for i in range(N - 1, -1, -1):
            psis[i] = psis[i + 1] * t / (float(i + 1.)) + 1.  # eq(6) of hk-relax
        x = {}  # Store x, r as dictionaries
        r = {}  # initialize residual
        Q = collections.deque()  # initialize queue

        r[(QID, 0)] = 1.
        Q.append((QID, 0))
        while len(Q) > 0:
            (v, j) = Q.popleft()  # v has r[(v,j)] ...
            rvj = r[(v, j)]
            # perform the hk-relax step
            if v not in x: x[v] = 0.
            x[v] += rvj
            r[(v, j)] = 0.
            mass = (t * rvj / (float(j) + 1.)) / self.deg[v]
            for u in self.adj_list[v]:  # for neighbors of v
                next = (u, j + 1)  # in the next block
                if j + 1 == N:
                    if u not in x: x[u] = 0.
                    x[u] += rvj / self.deg[v]
                    continue
                if next not in r: r[next] = 0.
                thresh = math.exp(t) * epsilon * self.deg[u]
                thresh = thresh / (2 * N * psis[j + 1])
                if r[next] < thresh and r[next] + mass >= thresh:
                    Q.append(next)  # add u to queue
                r[next] = r[next] + mass
        endtime = time.time()
        return x, endtime - starttime

    def maintain_connected(self, temp):
        q, visited = [QID], {QID}
        while q:
            v = q.pop()
            for u in self.adj_list[v]:
                if u in temp and u not in visited:
                    q.append(u)
                    visited.add(u)
        return visited

    def Global_CSM(self, C, l, h):
        deg = {}
        mymin_heap = MinHeap([])
        for u in C:
            deg[u] = len(set(self.adj_list[u]) & C)
            mymin_heap.insert([u, deg[u]])
        best = mymin_heap.peek()[1]
        D, best_index = [], 0
        temp = C.copy()
        while mymin_heap:
            u = mymin_heap.remove()[0]
            temp.remove(u)
            if str(u) == str(QID):
                break
            D.append(u)
            for v in self.adj_list[u]:
                if v in temp:
                    deg[v] -= 1
                    mymin_heap.decrease_key(v, deg[v])
            if mymin_heap.peek()[1] > best:
                best = mymin_heap.peek()[1]
                best_index = len(D)
        R = C - set(D[:best_index])
        result = self.maintain_connected(R)
        while len(result) > h:
            if best_index >= len(D):
                break
            v = D[best_index]
            R.add(v)
            result = self.maintain_connected(R)
            best_index += 1

        if l <= len(result) <= h:
            return result

        if len(result) < l or len(result) > h:
            result = set()
            return result

    def Local_CSM(self, l, h):
        starttime = time.time()
        best, C, D = 0, set(), {QID}
        inter_degree, inter_degree_min_heap = {}, MinHeap([])
        Q = MaxHeap([])
        Q.insert([QID, 0])
        Q_with_C = {}
        inter_degree_max_min = 0
        min_conductance = float('inf')
        S = set()

        while Q.heap:
            v = Q.remove()[0]
            C.add(v)
            temp = set(self.adj_list[v]) & C
            inter_degree[v] = len(temp)
            inter_degree_min_heap.insert([v, inter_degree[v]])
            for u in temp:
                inter_degree[u] = inter_degree[u] + 1
                inter_degree_min_heap.increase_key(u, inter_degree[u])
            for w in self.adj_list[v]:
                if w in Q.heap_dict:
                    Q_with_C[w] = Q_with_C[w] + 1
                    Q.increase_key(w, Q_with_C[w])
            if (inter_degree_min_heap.peek())[1] > best:
                best = (inter_degree_min_heap.peek())[1]
                if (best == len(self.adj_list[QID]) or best == math.floor((1 + math.sqrt(9 + 8 * (self.m - self.n))) / 2)) and l <= len(C) <= h:
                    min_conductance = self.compu_conductance(C)
                    endtime = time.time()
                    return C, endtime - starttime , min_conductance

            if l <= len(C) <= h:
                if inter_degree_max_min < (inter_degree_min_heap.peek())[1]:
                    inter_degree_max_min = (inter_degree_min_heap.peek())[1]
                    S = C.copy()
                if len(C) == h:
                    if len(S) == 0:
                        min_conductance = 0
                    else:
                        min_conductance = self.compu_conductance(S)
                    endtime = time.time()
                    return S, endtime - starttime, min_conductance

            for u in self.adj_list[v]:
                if u not in D:
                    D.add(u)
                    if len(self.adj_list[u]) >= best:
                        Q_with_C[u] = len(set(self.adj_list[u]) & C)
                        Q.insert([u, Q_with_C[u]])
        H = self.Global_CSM(C, l, h)
        if len(H) == 0:
            min_conductance = 0
            endtime = time.time()
            return H, endtime - starttime, min_conductance
        else:
            min_conductance = self.compu_conductance(H)
            endtime = time.time()
            return H, endtime - starttime, min_conductance

    def calculata_modularity(self, S):
        S = set(S)
        ind = 0
        outd = 0
        degree = {}
        cut = {}
        for u in S:
            ind += len(set(self.adj_list[u] & S))
            degree[u] = len(self.adj_list[u])
            cut[u] = degree[u] - len(set(self.adj_list[u]) & set(S))
            outd += cut[u]
        ind = ind / 2
        M = ind / outd
        return M


if __name__ == "__main__":
    dataset = sys.argv[1]   # Network
    T = int(sys.argv[2])    # iterations
    QID = int(sys.argv[3])  # query vid
    l = int(sys.argv[4])    # size LB
    h = int(sys.argv[5])    # size UB
    G = Graph(dataset)

    total_nodes = len(G.adj_list)
    total_degree = sum(len(neighbors) for neighbors in G.adj_list.values())
    average_degree = total_degree / total_nodes
    seed_nodes_density = [node for node, neighbors in G.adj_list.items() if len(neighbors) > average_degree]

    selected_QID = set()
    while len(selected_QID) < 50:
        QID = int(choice(list(seed_nodes_density)))
        selected_QID.add(QID)

    alpha = 0.01
    t = 5
    epsilon = 1 / G.m

    i = 0
    while i < 1:
        i += 1
        if i > 1:
            l += 3
            h += 3

        avg_csm_time = 0
        avg_csm_condu = 0
        avg_csm_size = 0
        avg_csm_M = 0

        avg_HK_time = 0
        avg_HK_condu = 0
        avg_HK_size = 0
        avg_HK_M = 0

        avg_push_time = 0
        avg_push_condu = 0
        avg_push_size = 0
        avg_push_M = 0

        avg_cr_time = 0
        avg_cr_condu = 0
        avg_cr_size = 0
        avg_cr_M = 0

        avg_dr_time = 0
        avg_dr_condu = 0
        avg_dr_size = 0
        avg_dr_M = 0

        num = 0
        for QID in selected_QID:
            csm_set, time_csm, csm_condu = G.Local_CSM(l, h)
            avg_csm_time += time_csm
            avg_csm_condu += csm_condu
            avg_csm_size += len(csm_set)
            if len(csm_set) != 0:
                num += 1
                csm_M = G.calculata_modularity(csm_set)
                avg_csm_M += csm_M

            pii, HK_time = G.hk_relax(t, epsilon)
            sweep_time_hk, HK_result, HK_condu = G.sweep_cut(pii, l, h)
            avg_HK_time += HK_time + sweep_time_hk
            avg_HK_condu += HK_condu
            avg_HK_size += len(HK_result)
            HK_M = G.calculata_modularity(HK_result)
            avg_HK_M += HK_M

            pi, fpush_time = G.fpush(alpha, epsilon)
            sweep_time_ppr, ppr_result, ppr_condu = G.sweep_cut(pi, l, h)
            avg_push_time += fpush_time + sweep_time_ppr
            avg_push_condu += ppr_condu
            avg_push_size += len(ppr_result)
            push_M = G.calculata_modularity(ppr_result)
            avg_push_M += push_M

            greedy_time_cr, greedy_set_cr, greedy_condu_cr = G.local_greedy_cr(T, l, h)
            avg_cr_time += greedy_time_cr
            avg_cr_condu += greedy_condu_cr
            avg_cr_size += len(greedy_set_cr)
            crr_M = G.calculata_modularity(greedy_set_cr)
            avg_cr_M += crr_M

            greedy_time_dr, greedy_set_dr, greedy_condu_dr = G.local_greedy_dr(T, l, h)
            avg_dr_time += greedy_time_dr
            avg_dr_condu += greedy_condu_dr
            avg_dr_size += len(greedy_set_dr)
            dr_M = G.calculata_modularity(greedy_set_dr)
            avg_dr_M += dr_M

        print("avg_csm_time:" + str(avg_csm_time / len(selected_QID)))
        print("avg_csm_condu:" + str(avg_csm_condu / num))
        print("avg_csm_size:" + str(avg_csm_size / num))
        print("avg_csm_M:" + str(avg_csm_M / num) + '\n')

        print("avg_HK_time:" + str(avg_HK_time / len(selected_QID)))
        print("avg_HK_condu:" + str(avg_HK_condu / len(selected_QID)))
        print("avg_HK_size:" + str(avg_HK_size / len(selected_QID)))
        print("avg_HK_M:" + str(avg_HK_M / len(selected_QID)) + '\n')

        print("avg_push_time:" + str(avg_push_time / len(selected_QID)))
        print("avg_push_condu:" + str(avg_push_condu / len(selected_QID)))
        print("avg_push_size:" + str(avg_push_size / len(selected_QID)))
        print("avg_push_M:" + str(avg_push_M / len(selected_QID)) + '\n')

        print("avg_cr_time:" + str(avg_cr_time / len(selected_QID)))
        print("avg_cr_condu:" + str(avg_cr_condu / len(selected_QID)))
        print("avg_cr_size:" + str(avg_cr_size / len(selected_QID)))
        print("avg_cr_M:" + str(avg_cr_M / len(selected_QID)) + '\n')

        print("avg_dr_time:" + str(avg_dr_time / len(selected_QID)))
        print("avg_dr_condu:" + str(avg_dr_condu / len(selected_QID)))
        print("avg_dr_size:" + str(avg_dr_size / len(selected_QID)))
        print("avg_dr_M:" + str(avg_dr_M / len(selected_QID)) + '\n')