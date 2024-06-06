import r2pipe
import json
import datetime
import collections
import numpy as np
from scipy.sparse import csc_matrix


def r2_wait_result(r2, cmd):
    res = r2.cmd(cmd)
    start = datetime.datetime.now()
    while True:
        tmp = r2.cmd('e""')
        if len(tmp) == 0:
            break
        else:
            res += tmp
        if (datetime.datetime.now()-start).seconds >= 10:
            break
    return res


def r2_to_dict(r2, cmd):
    res = r2_wait_result(r2, cmd)
    r2j = json.loads(res)
    return r2j


def do_analysis(r2p):
    # analyze all (aa)
    # analyze function calls (aac)
    print('analyzing binary...')
    # r2p.cmd('aa;aac')
    r2_wait_result(r2p, 'aaaa')
    data = r2_to_dict(r2p, 'iIj')
    bits = data.get('bits', '')
    print('done')
    if bits != '32':
        return '64'
    return '32'


def get_functions(r2p):
    # analyze function list json
    # afllj contents:
    # | address        start address
    # | size           function size (realsize)
    # | nbbs           number of basic blocks
    # | edges          number of edges between basic blocks
    # | cc             cyclomatic complexity ( cc = edges - blocks + 2 * exit_blocks)
    # | cost           cyclomatic cost
    # | min bound      minimal address
    # | range          function size
    # | max bound      maximal address
    # | calls          number of caller functions
    # | locals         number of local variables
    # | args           number of function arguments
    # | xref           number of cross references
    # | frame          function stack size
    # | name           function name
    print('retrieving function data')
    data = r2_to_dict(r2p, 'afllj')
    return data


def get_adjacency_data(r2p, data):
    od = collections.OrderedDict()
    print(f'processing {len(data)} items')
    for func in data:
        offset = func.get('offset', 0)
        tmp = []
        # 'fcn_addr' is address of calling function
        refs = r2_to_dict(r2p, f'axtj {offset}')
        for ref in refs:
            tmp.append(ref.get('fcn_addr', 0))
        od[offset] = tmp
    return od


def get_adjacency_matrix(adjdict):
    adjlist = list(adjdict.keys())
    num_funcs = len(adjlist)
    adjmat = np.zeros(shape=(num_funcs, num_funcs))
    for i in range(0, num_funcs):
        key = adjlist[i]
        vals = adjdict[key]
        for v in vals:
            if v not in adjdict.keys():
                continue
            idx = adjlist.index(v)
            if idx == i:
                # ignore recursive functions
                continue
            # yes, for now we don't care about duplicates
            # print(f'adding reference: caller {v} to callee {key}')
            adjmat[i][idx] = 1
    return adjmat


def get_rank(G, alpha=.85, maxerr=1e-06, loops=100):
    n = G.shape[0]
    M = csc_matrix(G, dtype=np.float64)
    rsums = np.array(M.sum(1))[:, 0]
    ri, ci = M.nonzero()
    M.data /= rsums[ri]
    sink = rsums == 0

    ro = np.zeros(n)
    r = np.ones(n)
    last_sum = 9999999
    for _ in range(loops):
        if last_sum < (n*maxerr):
            break
        ro = r.copy()
        for i in range(0, n):
            Ii = np.array(M[:, i].todense())[:, 0]
            Si = sink / float(n)
            Ti = np.ones(n) / float(n)
            r[i] = ro.dot(Ii*alpha + Si*alpha + Ti*(1-alpha))
        csum = np.sum(np.abs(r-ro))
        if csum < last_sum:
            last_sum = csum
        else:
            break
    return r/sum(r)


def get_rank_2(G, alpha=0.85, maxerr=1e-6, max_iter=100):
    # updated version
    n = G.shape[0]
    M = csc_matrix(G, dtype=np.float64)

    rsums = np.array(M.sum(1))[:, 0]
    M.data /= rsums[M.nonzero()[0]]

    sink = (rsums == 0)

    # Initialize PageRank values
    r = np.ones(n) / n

    for _ in range(max_iter):
        r_old = r.copy()

        # Compute the new PageRank values
        sink_contrib = alpha * r_old[sink].sum() / n
        r = alpha * M.dot(r_old) + sink_contrib + (1 - alpha) / n

        # Check convergence
        if np.linalg.norm(r - r_old, ord=1) < maxerr:
            break

    return r / r.sum()


def main():
    r2p = r2pipe.open()
    # analyze the binary
    bits = do_analysis(r2p)
    # get functions
    data = get_functions(r2p)
    # get adjacency data
    adjdict = get_adjacency_data(r2p, data)
    # build adjacency matrix
    adjmat = get_adjacency_matrix(adjdict)
    np.set_printoptions(threshold=np.inf)
    # print(adjmat)
    r = get_rank_2(adjmat)
    result = {}
    for i in range(len(list(adjdict.keys()))):
        result[list(adjdict.keys())[i]] = r[i]
    # sort and print the results
    result = sorted(result.items(), key=lambda x: x[1])[::-1]
    for item in result:
        addr = '0x{0:0>16x}'.format(item[0])
        print(f'function: {addr} score: {item[1]}')
    return


if __name__ == '__main__':
    main()

# testing code
#     G =  np.array([
#     [0, 1, 1, 0],
#     [1, 0, 0, 1],
#     [0, 1, 0, 1],
#     [1, 0, 1, 0]
# ])
#     ranks = get_rank(G)
#     print(ranks)
#     ranks = get_rank_2(G)
#     print(ranks)
