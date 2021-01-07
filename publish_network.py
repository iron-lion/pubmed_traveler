#  ncbi api page: https://www.ncbi.nlm.nih.gov/research/bionlp/APIs/
#  I used biopython - Entrez lib (pip install biopython).
#
# toy version (v0.1)
# Youngjun Park (youngjun.park.bio@gmail.com)

import sys
import networkx as nx
import numpy as np
import json
import time
from itertools import chain
from Bio import Entrez
import matplotlib.pyplot as plt
from collections import Counter

import my_config

# if your api_key is not valid, connection could be slower.
# you can get your api_key @ ncbi api pages. see above comment.
Entrez.email = my_config.MY_EMAIL
Entrez.api_key = my_config.MY_NCBI_API_KEY

# valid api_key = 10 request / 1 sec, if not, 3 request / 1 sec.
# sleep time is just for safety
if (len(Entrez.api_key) > 1):
    SLEEP_TIME = 0.11
else:
    SLEEP_TIME = 0.34
#


def get_abstract(pmid):
    handle = Entrez.efetch(db='pubmed', id=pmid, retmode='text', rettype='medline')
    return handle.read()


#
# progenitor : older paper
def get_ref_progenitor(pmid):
    time.sleep(SLEEP_TIME)
    link_list = []
    links = Entrez.elink(dbfrom="pubmed", id=pmid, linkname="pubmed_pubmed_refs")
    record = Entrez.read(links)
    links.close()
    rec = record[0]["LinkSetDb"]

    try:
        records = record[0][u'LinkSetDb'][0][u'Link']

        for link in records:
            link_list.append(link[u'Id'])

        return link_list
    except:
        # print("-d")
        return []


#
# descendant : newer paper
def get_ref_descendant(pmid):
    time.sleep(SLEEP_TIME)
    link_list = []
    links = Entrez.elink(dbfrom="pubmed", id=pmid, linkname="pubmed_pubmed_citedin")
    record = Entrez.read(links)
    links.close()
    rec = record[0]["LinkSetDb"]

    try:
        records = record[0][u'LinkSetDb'][0][u'Link']

        for link in records:
            link_list.append(link[u'Id'])

        return link_list
    except:
        # print("-p")
        return []


#
# simplest way to identify REVIEW PAPER
def find_review(given_text):
    if((given_text.lower().find("review") != -1) or (given_text.lower().find("assessment") != -1)):
        return True
    else:
        return False


#
# do filter work.
# publication network was denser network than i expected.
def scan_neighbor_common_breadth_first(g, l, func_p, hit_list):
    # 1. neighbor scan
    next_list = []
    for (parent, el) in l:
        neighbors = func_p(el)
        next_list += [(el, n) for n in neighbors]

    next_list = list_to_common_list(next_list)

    # 2. add neighbor to graph
    for (p, n) in next_list:
        tmp_nodes = g.nodes
        if (n not in tmp_nodes):
            is_it_rev = find_review(get_abstract(n))
            if (is_it_rev is True):
                hit_list.append(n)
                # next_list.remove((p, n))
                # g.add_edge(p, n)
                # return next_list ### test you
            g.add_node(n, review=is_it_rev)
            g.add_edge(p, n)

    return next_list


#
# without filter work.
def scan_neighbor_breadth_first(g, l, func_p, hit_list):
    # 1. neighbor scan
    next_list = []
    for (parent, el) in l:
        neighbors = func_p(el)
        next_list += [(el, n) for n in neighbors]

    # 2. add neighbor to graph
    for (p, n) in next_list:
        tmp_nodes = g.nodes
        if (n not in tmp_nodes):
            is_it_rev = find_review(get_abstract(n))
            if (is_it_rev is True):
                hit_list.append(n)
                next_list.remove((p, n))
                # g.add_edge(p, n)
                # return next_list ### test you
            else:
                g.add_node(n, review=is_it_rev)
            g.add_edge(p, n)

    return next_list


#
# publication filter
# scan common one, first. Maybe we need different scan policy.
def list_to_common_list(q_list):
    next_dict = Counter([lx for _, lx in q_list])
    next_dict = {k : v for k, v in dict(next_dict).items() if v > my_config.COMMONALITY }

    threshold = my_config.COMMONALITY
    while (len(next_dict) > my_config.SEARCH_LIMIT):
        print ("WARNING: too many papers found, reconfigure threshold COMMONALITY: ", threshold)
        threshold += 1
        next_dict = {k : v for k, v in dict(next_dict).items() if v > threshold }

    if (len(next_dict) == 0):
        next_dict = {k : v for k, v in dict(next_dict).items() if v > threshold-1 }
        if (len(next_dict) == 0):
            print ("WARNING: leaf node")
            return []

    selected = next_dict.keys()

    r_list = [(a, b) for a, b in q_list if b in selected]
    return r_list[:my_config.SEARCH_LIMIT]


#
def init_building_network(do_g, up_g, paper_list):
    # First neighbor : _ :
    p1_list = []
    d1_list = []
    for p in paper_list:
        up_g.add_node(p, review=find_review(get_abstract(p)))
        do_g.add_node(p, review=find_review(get_abstract(p)))

        p1_list += [(p, n) for n in get_ref_progenitor(p)]
        d1_list += [(p, n) for n in get_ref_descendant(p)]

        for pe in p1_list:
            up_g.add_edge(pe[0], pe[1])
        for pe in d1_list:
            do_g.add_edge(pe[0], pe[1])

    p1_list = list_to_common_list(p1_list)
    d1_list = list_to_common_list(d1_list)

    cn_p = nx.number_connected_components(up_g)
    cn_d = nx.number_connected_components(do_g)
    print("initial network's progeneitor cluster # : ", cn_p)
    print("initial network's descendant cluster # : ", cn_d)

    return d1_list, p1_list, cn_d, cn_p


#
def expand_one_network(p_list, p_g, func_p, review_list):
    # Second +  generation search
    # next_list = scan_neighbor_breadth_first(p_g, p_list, func_p, review_list)
    next_list = scan_neighbor_common_breadth_first(p_g, p_list, func_p, review_list)
    cn = nx.number_connected_components(p_g)
    return next_list, cn


#
# result graph drawing part
def has_review_hub(gg, paper_list):
    g = nx.Graph()
    for p1 in paper_list:
        if p1 not in g.nodes:
            g.add_node(p1, review=gg.nodes[p1]['review'])
        for p2 in paper_list:
            if p1 is p2:
                continue

            paths = nx.all_simple_paths(gg, p1, p2)
            for p in paths:
                found = False
                ps = p[1:-1]
                for node in ps:
                    if (my_config.REVIEW_CHECK is True):
                        try:
                            if (gg.nodes[node]['review'] is True):
                                found = True
                        except:
                            gnode = gg.nodes[node]
                            gnode['review'] = find_review(get_abstract(node))
                            if (gnode['review'] is True):
                                # print("-", gnode)
                                found = True
                    else:  # Pass! result
                        found = True
                if (found is True):
                    for i in range(len(ps)-1):
                        g.add_edge(ps[i], ps[i+1])
    ncc = nx.number_connected_components(g)
    # print("result graph : ", ncc)
    return g, ncc


#
# not well implemented.
def network_save(gg, ppath):
    # nx.write_gpickle(gg, ppath)
    cyjs = nx.cytoscape_data(gg)
    with open(ppath, 'w') as fds:
        json.dump(cyjs, fds)


def network_load(ppath):
    # gg = nx.read_gpickle(ppath)

    # cyjs = np.load(ppath, allow_pickle=True).item()
    with open(ppath, 'r') as fds:
        cyjs = json.load(fds)
    return nx.cytoscape_graph(cyjs)


def init_network(dec_path, prog_path):
    g_d = nx.Graph()
    g_p = nx.Graph()
    """
    if os.path.exists(dec_path):
        g_d = network_load(dec_path)
        print( "loaded graph for descendent" )
    else:
        g_d = nx.Graph()

    if os.path.exists(prog_path):
        g_p = network_load(prog_path)
        print( "loaded graph for progenitor" )
    else:
        g_p = nx.Graph()
    """
    return g_d, g_p


#
# draw simple result network image on pdf
def print_two_result(dg, pg):
    options = {
        "font_size": 10,
        "node_size": 200,
        "node_color": "white",
        "edgecolors": "black",
        "linewidths": 3
    }

    plt.subplot(2, 1, 1)
    plt.margins(0.2)
    nx.draw_networkx(pg, **options)
    plt.axis("off")
    plt.ylabel('progenitor')

    node_labels = nx.get_node_attributes(dg, 'review')
    plt.subplot(2, 1, 2)
    plt.margins(0.2)
    nx.draw_networkx(dg, **options)
    plt.axis("off")
    plt.ylabel('descendant')

    plt.savefig(my_config.RESULT_PDF_PATH)


# MAIN ---- ---- ----
# ---- MAIN ---- ----
# ---- ---- MAIN ----
# ---- ---- ---- MAIN

if __name__ == '__main__':
    paper_list = sys.argv[1:]

    if (len(paper_list) == 0):
        paper_list = ["27883053", "29606308", "22258609", "31510697"]
    print("Query papers : ", paper_list)

    review_list = []
    g_d, g_p = init_network(my_config.GRAPH_D_PATH, my_config.GRAPH_P_PATH)

    # start
    first_dec_list, first_prog_list, n_cluster_d, n_cluster_p = init_building_network(g_d, g_p, paper_list)

    # progenitor
    prog_list = first_prog_list
    result_p_g, found_p_ncc = has_review_hub(g_p, paper_list)
    print("initial progenitor result components #:", found_p_ncc)
    ct = 1
    while ( (found_p_ncc > 1) and (len(prog_list) > 0) and (ct < my_config.LIMIT_ROUND) ):
        print("search round: ", ct, " - query_pmid #: ", len(prog_list),
              " - current component #: ", n_cluster_p,
              " - result component #: ", found_p_ncc)
        prog_list, n_cluster_p = expand_one_network(prog_list, g_p, get_ref_progenitor, review_list)
        result_p_g, found_p_ncc = has_review_hub(g_p, paper_list)
        ct += 1

    # descendant
    dec_list = first_dec_list
    result_d_g, found_d_ncc = has_review_hub(g_d, paper_list)
    print("initial descendant result components #: ", found_d_ncc)
    ct = 1
    while ( (found_d_ncc > 1) and (len(dec_list) > 0) and (ct < my_config.LIMIT_ROUND) ):
        print("search round: ", ct, " - query_pmid #: ", len(dec_list),
              " - current component #: ", n_cluster_d,
              " - result component #: ", found_d_ncc)
        dec_list, n_cluster_d = expand_one_network(dec_list, g_d, get_ref_descendant, review_list)
        result_d_g, found_d_ncc = has_review_hub(g_d, paper_list)
        ct += 1

    # end, save
    print_two_result(result_d_g, result_p_g)
    with open(my_config.LOG_PATH, 'a+') as log_fd:
        print(','.join(paper_list), " --> ",
              my_config.RESULT_PDF_PATH, " | ",
              my_config.GRAPH_D_PATH, ", ",
              my_config.GRAPH_P_PATH, file=log_fd)

    """
    network_save(g_d, "./tmp.d.cyjs")
    network_save(g_p, "./tmp.p.cyjs")
    """
    network_save(result_p_g, my_config.GRAPH_P_PATH)
    network_save(result_d_g, my_config.GRAPH_D_PATH)

    print(">______________")
    print("I found review keywords : ", len(review_list))
