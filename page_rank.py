import r2pipe
import json
import datetime
import networkx as nx
import pandas as pd


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
    r2_wait_result(r2p, 'aac')
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


def get_adjacency_dict(data):
    # parse the json and return a list of function refs
    xrefs = []
    # for each of our functions
    for func in data:
        # find all of the references within the function
        for ref in func.get('codexrefs', []):
            # make sure we have a call reference
            if ref.get('type', '') == 'CALL':
                # for each of our functions
                for func2 in data:
                    # make sure the call reference is within the bounds of the destination function
                    if ref['addr'] >= func2['minbound'] and ref['addr'] <= func2['maxbound']:
                        # add it to our results if we don't already have it
                        d = {'to': func['offset'], 'from': func2['offset']}
                        if d not in xrefs:
                            xrefs.append(d)
    return xrefs


def get_page_rank(funcs):
    df = pd.DataFrame(funcs, columns=['to', 'from'])
    G = nx.from_pandas_edgelist(df, 'to', 'from')
    rank = nx.pagerank(G)
    sorted_rank = [k for k in sorted(rank.items(), reverse=True, key=lambda item:item[1])]
    return sorted_rank


def main():
    r2p = r2pipe.open()
    # analyze the binary
    bits = do_analysis(r2p)
    # get functions
    data = get_functions(r2p)
    # build adjacency dict
    xrefs = get_adjacency_dict(data)
    # get sorted page rank scores
    pr = get_page_rank(xrefs)
    # print the results
    for item in pr:
        if bits == '64':
            addr = '0x{0:0>16x}'.format(item[0])
        else:
            addr = '0x{0:0>08x}'.format(item[0])
        print(f'address: {addr}, score: {item[1]}')


if __name__ == '__main__':
    main()
