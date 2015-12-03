# -*- coding: utf-8 -*-
#!/bin/python
import time

class Vertex_E(object):
    # 企业编号+'E','E',信用值，信用等级
    def __init__(self, idx, prop, creditscore, rating):
        self.idx = idx
        self.prop = prop
        self.creditscore = creditscore
        self.rating = rating

class Vertex_P(object):
    # 人编号+'P','P',身份证号
    def __init__(self, idx, prop, cerno):
        self.idx = idx
        self.prop = prop
        self.cerno = cerno

class Link(object):
    def __init__(self, link_weight, link_property):
        self.link_weight = link_weight
        self.link_property = link_property

class Graph(object):
    def __init__(self, v_list):
        self.creditscores = {}
        self.link_list = {}
        for vertex in v_list:
            if is_E(vertex):
                self.creditscores[vertex] = v_list[vertex].creditscore
            self.link_list[vertex] = {}
            link = Link(1.0, 'self')
            self.add_link(vertex, vertex, link)
    def add_link(self, v1, v2, link):
        if v2 in self.link_list[v1]:
            # 取权重最大的
            if self.link_list[v1][v2].link_weight < link.link_weight:
                self.link_list[v1][v2] = link
        else:
            self.link_list[v1][v2] = link

def is_P(v):
    return v[-1] == 'P'
def is_E(v):
    return v[-1] == 'E'

def load_vertex_E(vertex_fn, prop, v_list):
    vertex_file = open(vertex_fn, 'r')
    vertex_file.readline()
    blank_node_cnt = 0
    for line in vertex_file.readlines():
        vertex = line.strip().split(',')
        vertex[0] = vertex[0].strip('"')
        if vertex[9] == '':
            v = Vertex_E(vertex[0]+prop, prop, -1, 'E')
            blank_node_cnt = blank_node_cnt + 1
        else:
            v = Vertex_E(vertex[0]+prop, prop, float(vertex[9]), vertex[10])
        v_list[vertex[0]+prop] = v
    print "Blank E:" + str(blank_node_cnt)
    return v_list

def load_vertex_P(vertex_fn, prop, v_list):
    vertex_file = open(vertex_fn, 'r')
    vertex_file.readline()
    blank_node_cnt = 0
    for line in vertex_file.readlines():
        vertex = line.strip().split(',')
        vertex[0] = vertex[0].strip('"')
        if not vertex[1]:
            blank_node_cnt = blank_node_cnt + 1
        v = Vertex_P(vertex[0]+prop, prop, vertex[1])
        v_list[vertex[0]+prop] = v
    print "Blank P:" + str(blank_node_cnt)
    return v_list

def load_link(link_fn, graph):
    link_file = open(link_fn, 'r')
    link_file.readline()
    for line in link_file.readlines():
        link = line.strip().split(',')
        v1 = link[0].strip('"')+link[1].strip('"')
        v2 = link[2].strip('"')+link[3].strip('"')
        link_weight = float(link[4].strip('"'))
        link_prop_pos = link[5].strip('"')
        if link_prop_pos == 'FATHER':
            link_prop_neg = 'SON'
        elif link_prop_pos == 'SON':
            link_prop_neg = 'FATHER'
        elif link_prop_pos == 'OTHER':
            link_prop_neg = 'OTHER'
        else:
            print 'wrong link property!!!'
        if link_weight != 0:
            link_pos = Link(link_weight, link_prop_pos)
            link_neg = Link(link_weight, link_prop_neg)
            graph.add_link(v1, v2, link_pos)
            graph.add_link(v2, v1, link_neg)
    return graph
 
def pagerank(graph, v_list, damping=0.85, epsilon=1.0e-8):
    outlink_map = {}
    inlink_counts = {}
    
    # 初始化
    def new_node(node):
        if node not in outlink_map: 
            outlink_map[node] = set()
        if node not in inlink_counts: 
            inlink_counts[node] = 0
    
    # 计算加权出度入度
    for node in v_list:
        new_node(node)
        for nb_node in graph.link_list[node].keys():
            new_node(nb_node)
            if node == nb_node:
                continue

            if nb_node not in outlink_map[node]:
                if graph.link_list[node][nb_node].link_property == "OTHER" or graph.link_list[node][nb_node].link_property == "FATHER":
                    outlink_map[node].add(nb_node)
                    inlink_counts[nb_node] += graph.link_list[node][nb_node].link_weight
    print "degree compute finish!"
    all_nodes = set(outlink_map.keys())

    # 初始化权威值向量PageRank为1/n
    initial_value = 1.0 / len(all_nodes)
    print initial_value
    PageRanks = {}
    for node in outlink_map.keys():
        PageRanks[node] = initial_value
    print "init PageRanks finish!"

    # 迭代
    new_PageRanks = {}
    delta = 1.0
    n_iterations = 0
    
    print "begin to iteration!"
    virtual_node_rank = sum(PageRanks[node] * (1.0 - damping)/len(all_nodes) for node in all_nodes)
    print virtual_node_rank
    while delta > epsilon:
        new_PageRanks = {}
        sum_pagerank = 0
        for node, outlinks in outlink_map.items():
            new_PageRanks[node] = virtual_node_rank + sum(PageRanks[outlink]* (damping * graph.link_list[node][outlink].link_weight / inlink_counts[outlink]) for outlink in outlinks)
            sum_pagerank += new_PageRanks[node]
        print "=========== " + str(n_iterations) + "th iteration ==========="
        print "sum_pagerank:"
        print sum_pagerank
        delta = sum(abs(new_PageRanks[node] - PageRanks[node]) for node in new_PageRanks.keys())
        print "delta:"
        print delta
        PageRanks, new_PageRanks = new_PageRanks, PageRanks
        n_iterations += 1
    
    return PageRanks, n_iterations       

def groupscore_pagerank(PageRanks,graph,v_list):
    groupscores = {}
    for node in v_list:
        # 只计算E
        if is_P(node):
            continue
        # 分子
        denominator = 0
        # 分母
        numerator = 0
        # 二级节点集合
        # nb_nb_node: [nb_node,denominator,numerator]
        nb_nb_node_dict = {}
        # 一级邻边节点
        for nb_node in graph.link_list[node].keys():
            # 只有E有独立信用分，且分非负
            if (is_E(nb_node) and graph.creditscores[nb_node]>=0):
                denominator += graph.creditscores[nb_node] * graph.link_list[node][nb_node].link_weight * PageRanks[nb_node]
                numerator += graph.link_list[node][nb_node].link_weight * PageRanks[nb_node]
            # 二级邻边节点
            for nb_nb_node in graph.link_list[nb_node].keys():
                # 只有E有独立信用分
                if is_P(nb_nb_node):
                    continue
                if graph.creditscores[nb_nb_node]<0:
                    continue
                # 刨去node及一级邻边节点
                if (nb_nb_node != node and nb_nb_node not in graph.link_list[node].keys()):
                    temp_weight = graph.link_list[node][nb_node].link_weight * PageRanks[nb_node] * graph.link_list[nb_node][nb_nb_node].link_weight * PageRanks[nb_nb_node]
                    if nb_nb_node in nb_nb_node_dict.keys():
                        if temp_weight > nb_nb_node_dict[nb_nb_node][2]:
                            nb_nb_node_dict[nb_nb_node] = [nb_node, graph.creditscores[nb_nb_node] * temp_weight, temp_weight]
                    else:
                        nb_nb_node_dict[nb_nb_node] = [nb_node, graph.creditscores[nb_nb_node] * temp_weight, temp_weight]
        for v in nb_nb_node_dict.values():
            denominator += v[1]
            numerator += v[2]
        if denominator == 0:
            groupscores[node] = ""
            continue
        groupscores[node] = denominator/numerator
    return groupscores

if __name__ == '__main__':
    v_list = {}
    v_list = load_vertex_E("./data/EINFOALL_ANON.csv", 'E', v_list)
    print "total E: " + str(len(v_list))
    v_list = load_vertex_P("./data/PINFOALL_ANON.csv", 'P', v_list)
    print "total node: " + str(len(v_list))

    graph = Graph(v_list)
    graph = load_link("./data/LINK_ANON.csv", graph)
    
    print "load finish!"

    PageRanks, n_iterations = pagerank(graph, v_list, 0.8, 1.0e-8)

    # pagerank结果输出
    print "=========== PageRanks ==========="
    fs = open("pagerank.txt", 'w')
    fs.close()
    fs = open("pagerank.txt", 'a')
    fs.write(str(PageRanks))
    fs.close()
    # print PageRanks
    print "n_iterations:"
    print n_iterations

    # 计算“群体分”
    groupscores = groupscore_pagerank(PageRanks,graph,v_list)
    max_neg = 0
    thre = -100
    alert_count = 0
    fs = open("groupscores_pagerank.txt", 'w')
    fs.close()
    for node in groupscores:
        if groupscores[node]=="":
            print node + "\t" + str(groupscores[node]) + "\t" + str(graph.creditscores[node])
            fs = open("groupscores_pagerank.txt", 'a')
            fs.write( node + "\t" + str(groupscores[node]) + "\t" + str(graph.creditscores[node]) + "\n")
            fs.close()
        else:
            print node + "\t" + str(groupscores[node]) + "\t" + str(graph.creditscores[node]) + "\t" + str(groupscores[node]-graph.creditscores[node])
            fs = open("groupscores_pagerank.txt", 'a')
            fs.write( node + "\t" + str(groupscores[node]) + "\t" + str(graph.creditscores[node]) + "\t" + str(groupscores[node]-graph.creditscores[node]) + "\n")
            fs.close()
            if groupscores[node]-graph.creditscores[node]<max_neg:
                max_neg = groupscores[node]-graph.creditscores[node]
            if groupscores[node]-graph.creditscores[node]<thre:
                alert_count += 1
    print "=========== finish! ==========="
    print str(alert_count) + " companies " + str(thre) + " below independent score"


    print ""
    print ""
    print ""

    print "max_neg:\t"+ str(max_neg)
    # 查看
    max_neg_node = "anon_S3478E"
    print "=========== check:\t" + max_neg_node + " ==========="
    print max_neg_node + "\t" + str(groupscores[max_neg_node]) + "\t" + str(graph.creditscores[max_neg_node]) + "\t" + str(groupscores[max_neg_node]-graph.creditscores[max_neg_node])
    # 二级节点集合
    # nb_nb_node: [nb_node,denominator,numerator]
    nb_nb_node_dict = {}
    # 一级邻边节点
    for nb_node in graph.link_list[max_neg_node].keys():
        # 只有E有独立信用分，且分非负
        if (is_E(nb_node) and graph.creditscores[nb_node]>=0):
            weight = graph.link_list[max_neg_node][nb_node].link_weight * PageRanks[nb_node]
            print "level 1: " + nb_node + "\tweight:" + str(weight) + "\tcreditscore:" + str(graph.creditscores[nb_node])
        # 二级邻边节点
        for nb_nb_node in graph.link_list[nb_node].keys():
            # 只有E有独立信用分
            if is_P(nb_nb_node):
                continue
            if graph.creditscores[nb_nb_node]<0:
                continue
            # 刨去node及一级邻边节点
            if (nb_nb_node != max_neg_node and (nb_nb_node not in graph.link_list[max_neg_node].keys())):
                temp_weight = graph.link_list[max_neg_node][nb_node].link_weight * PageRanks[nb_node] * graph.link_list[nb_node][nb_nb_node].link_weight * PageRanks[nb_nb_node]
                if nb_nb_node in nb_nb_node_dict.keys():
                    if temp_weight > nb_nb_node_dict[nb_nb_node][2]:
                        nb_nb_node_dict[nb_nb_node] = [nb_node, graph.creditscores[nb_nb_node] * temp_weight, temp_weight]
                else:
                    nb_nb_node_dict[nb_nb_node] = [nb_node, graph.creditscores[nb_nb_node] * temp_weight, temp_weight]
    for k,v in nb_nb_node_dict.items():
        print "level 2: " + k + "\tweight:" + str(v[2]) + "\tcreditscore:" + str(graph.creditscores[k])

