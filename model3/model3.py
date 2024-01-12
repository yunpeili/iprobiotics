import os
from tqdm import tqdm
from collections import Counter
import json
import statistics

dna_list = ['A', 'G', 'T', 'C']

def all_possible_seq(k):
    
    def generate_seq(k, chars, current_string, results):
        if k == 0:
            results.append(current_string)
            return
        for char in chars:
            generate_seq(k - 1, chars, current_string + char, results)

    results = []
    generate_seq(k, dna_list, '', results)
    return results


ambiguity_symbol_dict = {
    'A': ['A'], 'T': ['T'], 'C': ['C'], 'G': ['G'],
    'W': ['A', 'T'], 'S': ['C', 'G'], 'M': ['A', 'C'], 
    'K': ['G', 'T'], 'R': ['A', 'G'], 'Y': ['C', 'T'], 
    'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'], 
    'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']
}

def nucleic_acid_expand(dna_seq):
    results = ['']
    for base in dna_seq:
        results = [ seq + new_seq 
                  for seq in results 
                    for new_seq in ambiguity_symbol_dict[base]]

    return results


import itertools
# from line_profiler import profile

expand_dict = {}

def nucleic_acid_expand_product(dna_seq):
    if dna_seq in expand_dict:
        return expand_dict[dna_seq]
    else:
        bases = [ambiguity_symbol_dict[base] for base in dna_seq]
        combinations = [''.join(comb) for comb in itertools.product(*bases)]
        expand_dict[dna_seq] = combinations
    return combinations

# with open('all_possible_seq.txt', 'r') as f:
#     lines = f.readlines()
#     all_possible_list = []
#     for line in lines:
#         all_possible_list.append(line)



# Test 
# print(nucleic_acid_count('AWCBN'))
# print(nucleic_acid_expand('AWCBN'))

#创建kmer序列与对应tf的字典
# @profile
# def get_dict(seq):
#     curr_dict = {}
#     for k in range(2, 9):
#         dna_segs = all_possible_seq(k)                #k-mer序列
#         new_dict = {seg:0 for seg in dna_segs}        #创建一个初始字典，key为k-mer序列，value为0
#         curr_dict.update(new_dict)                         #在k值更新之后更新字典
#         all_occur = len(seq) - k + 1                  #k-mer序列的数量
#         for x in range(all_occur):                    #计算每个k-mer序列的出现次数
#             key = seq[x: x+k]                         
#             # curr_dict[key] += 1
#             # curr_dict.update({ key: curr_dict[key]+1 })
#             # curr_dict.update({ key1: curr_dict[key]+0.5, key2: curr_dict[key]+0.5 })
#             # curr_dict.update({ lambda N: keyN: curr_dict[key]+1/N })
#             # ignore seqences which have more than 1 'N'
#             # if key.count('N') > 3:
#             #     continue
#             expand_result = nucleic_acid_expand_product(dna_seq=key)
#             for key in expand_result:
#                 curr_dict[key] += 1/len(expand_result)
        
#         # curr_dict = {k:v for k,v in curr_dict.items()}#字典储存在curr_dict 
        
#     # for k in list(curr_dict.keys()):
#     #     if curr_dict[k] == float(0):
#     #         del curr_dict[k]
            
#     l = sum(curr_dict.values())

#     for key in curr_dict.keys():                          
#         curr_dict[key] = curr_dict[key] / l           #计算tf
#     return curr_dict

# def get_dict(seq):
#     curr_dict = {}
#     k = 6
#     dna_segs = all_possible_seq(k)                #k-mer序列
#     new_dict = {seg:0 for seg in dna_segs}        #创建一个初始字典，key为k-mer序列，value为0
#     curr_dict.update(new_dict)                         #在k值更新之后更新字典
#     all_occur = len(seq) - k + 1                  #k-mer序列的数量

#     for x in range(all_occur):                    #计算每个k-mer序列的出现次数
#         key = seq[x: x+k]     
#         expand_result = nucleic_acid_expand_product(dna_seq=key)
#         for key in expand_result:
#             curr_dict[key] += 1/len(expand_result)

#     # for k in list(curr_dict.keys()):
#     #     if curr_dict[k] == float(0):
#     #         del curr_dict[k]

#     l = sum(curr_dict.values())
#     for key in curr_dict.keys():              
#         curr_dict[key] = curr_dict[key] / l           #计算tf
#     return curr_dict



###10.31 version
###first step of IFS, to generate our kmer:tf dictionary, prepared for abstracting the fine features
# def get_dict(seq):
#     init_dict = {}
#     k = 9
#     dna_segs = all_possible_seq(k)
#     curr_dict = {seg:0 for seg in dna_segs}
#     init_dict.update(curr_dict)
#     for x in range(len(seq)-k+1):
#         key = seq[x: x+k]
#         if key not in init_dict.keys():
#             continue
#         else:
#             init_dict[key] += 1

#     l = sum(init_dict.values())
#     for key in init_dict.keys():
#         init_dict[key] /= l
    
#     values = list(init_dict.values())
#     mean = statistics.mean(values)
#     std = statistics.stdev(values)
#     for key in init_dict:
#         value = init_dict[key]
#         normalized_value = (value-mean)/std
#         init_dict[key] = normalized_value

#     return init_dict

###11.2 version
###second step of IFS, to abastract the core features
def get_dict(seq):
    with open('/home/yunpei/probiotics/model3/v2_sorted_fine_feature5-9.json', 'r') as f:
        fine_feature_list = json.load(f)
        fine_feature_list = fine_feature_list[:85]
    init_dict = {seg:0 for seg in fine_feature_list}
    for k in range(5, 10):
        for x in range(len(seq)-k+1):
            key = seq[x: x+k]
            if key not in init_dict.keys():
                continue
            else:
                init_dict[key] += 1
    l = sum(init_dict.values())
    for k in init_dict.keys():
        init_dict[k] /= l
    
    values = list(init_dict.values())
    mean = statistics.mean(values)
    std = statistics.stdev(values)
    for key in init_dict:
        value = init_dict[key]
        normalized_value = (value-mean)/std
        init_dict[key] = normalized_value

    return init_dict


def feat2str(int_float_dict):
    feat_str = ''
    for index, value in int_float_dict.items():
        feat_str += str(index) + ':' + str(value) + ' '  
    #feat_str = str(label) + '\t' + feat_str 
    return feat_str

folder_path_0 = '/home/yunpei/probiotics/dataset/longest_pro_lacto'
label_0 = 0
dirs_0 = [path for path in sorted(os.listdir(folder_path_0))]
seq_list_0 = []
filename_list_0 =[]
for input_file in dirs_0:
# for input_file in dirs:
    with open(os.path.join(folder_path_0, input_file)) as seqs:
        dna_seq = seqs.read()
        seq_list_0.append(dna_seq)
        filename_list_0.append(input_file)

folder_path_1 = '/home/yunpei/probiotics/dataset/longest_nonpro_lacto'
label_1 = 1
dirs_1 = [path for path in sorted(os.listdir(folder_path_1))]
seq_list_1 = []
filename_list_1 =[]
for input_file in dirs_1:
# for input_file in dirs:
    with open(os.path.join(folder_path_1, input_file)) as seqs:
        dna_seq = seqs.read()
        seq_list_1.append(dna_seq)
        filename_list_1.append(input_file)

# folder_path_2 = '/home/yunpei/probiotics/dataset/longest_other_probiotics'
# label_2 = 2
# dirs_2 = [path for path in sorted(os.listdir(folder_path_2))]
# seq_list_2 = []
# filename_list_2 =[]
# for input_file in dirs_2:
# # for input_file in dirs:
#     with open(os.path.join(folder_path_2, input_file)) as seqs:
#         dna_seq = seqs.read()
#         seq_list_2.append(dna_seq)
#         filename_list_2.append(input_file)

import multiprocessing as mp
from multiprocessing import Pool, Queue, Manager, Process
# num_process = mp.cpu_count()//2
num_process = 32
print(f'using {num_process} cores')


def get_feat_str_list(seq_list, feat_str_list, queue, filename_list):
    for seq, filename in zip(seq_list, filename_list):
        feat_dict = get_dict(seq)
    # #给kmer序列标号并创建字典
    # index_ATCG_map = {ATCG: int(i) for i, ATCG in enumerate(feat_dict.keys())}
    # #创建index和tf的字典
    # feat_dict = {index_ATCG_map[ATCG]: value for i, (ATCG, value) in enumerate(feat_dict.items())}
        feat_dict = {i: value for i, (ATCG, value) in enumerate(feat_dict.items())}
        feat_str = feat2str(feat_dict)
        feat_str_list.append(feat_str)
        # print(filename)
        queue.put(1)

def tqdm_listener(total:int, queue:Queue):
    pbar = tqdm(total=total)
    while True:
        if not queue.empty():
            k = queue.get()
            if k == 1:
                pbar.update(1)
            else:
                break
    pbar.close()

def process_feat(seq_list, filename_list):

    manager = Manager()
    # all_feat_str_list = manager.list()
    feat_str_list = manager.list()
    tqdm_queue = manager.Queue()
    # iter_count = manager.Value('i', 0)
    len(seq_list)

    batch_size = len(seq_list) // num_process
    batches = [seq_list[i:i+batch_size] for i in range(0, len(seq_list), batch_size)]
    filename_batches = [filename_list[i:i+batch_size] for i in range(0, len(seq_list), batch_size)]

    tqdm_process = Process(target=tqdm_listener, args=(len(seq_list), tqdm_queue))
    tqdm_process.start()

    # start all process
    processes = []

    for batch, fb in zip(batches,filename_batches):
        p = Process(target=get_feat_str_list, args=(batch, feat_str_list, tqdm_queue, fb))
        processes.append(p)
        p.start()

    # wait all finish
    for p in processes:
        p.join()

    tqdm_queue.put(-1)
    tqdm_process.join()

    # with Pool(processes=num_process) as pool:
    #     results = []
    #     for dna_seq in seq_list:
    #         result = pool.apply_async(get_feat_str, args=(dna_seq, shared_expand_dict))
    #         results.append(result)

    #     with tqdm(total=len(seq_list)) as pbar:
    #         for result in results:
    #             feat_str_list.append(result.get())
    #             pbar.update(1)

            # for result in pool.imap(process_seq, seq_list):
            
    # pool.close()
    # pool.join()

    return list(feat_str_list)


pro_lacto_list = process_feat(seq_list_0, filename_list_0)
nonpro_lacto_list = process_feat(seq_list_1, filename_list_1)

with open('./model3/test.tsv', 'w') as f:
    for feat_str in pro_lacto_list:
        f.write(str(label_0) + '\t' + feat_str + '\n')
    for feat_str in nonpro_lacto_list:
        f.write(str(label_1) + '\t' + feat_str + '\n')