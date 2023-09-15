import os
from tqdm import tqdm
from collections import Counter

#生成所有k-mer序列
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
    'A': {'A': 1}, 'T': {'T': 1}, 'C': {'C': 1}, 'G': {'G': 1},
    'W': {'A': 0.5, 'T': 0.5}, 'S': {'C': 0.5, 'G': 0.5},
    'M': {'A': 0.5, 'C': 0.5}, 'K': {'G': 0.5, 'T': 0.5},
    'R': {'A': 0.5, 'G': 0.5}, 'Y': {'C': 0.5, 'T': 0.5},
    'B': {'C': 1/3, 'G': 1/3, 'T': 1/3}, 'D': {'A': 1/3, 'G': 1/3, 'T': 1/3},
    'H': {'A': 1/3, 'C': 1/3, 'T': 1/3}, 'V': {'A': 1/3, 'C': 1/3, 'G': 1/3},
    'N': {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
}
# e.g. 'W':['A', 'T'] --> 'A':0.5, 'T':0.5
#      'D':['A', 'G', 'T'] --> 'A':1/3, 'G':1/3, 'T':1/3
#      'A':['A'] --> 'A':1/1

print(ambiguity_symbol_dict)

def nucleic_acid_count(dna_seq):
    result = {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0}

    for char in dna_seq:
        possib_chars = ambiguity_symbol_dict[char]

        for pc, possibility in possib_chars.items():
            result[pc] += possibility
        
    return result

# Test 
# print(nucleic_acid_count('AWCBN'))

#创建kmer序列与对应tf的字典
def get_dict(seq):
    curr_dict = {}
    for k in range(2, 5):
        dna_segs = all_possible_seq(k)                #k-mer序列
        new_dict = {seg:0 for seg in dna_segs}        #创建一个初始字典，key为k-mer序列，value为0
        curr_dict += new_dict                         #在k值更新之后更新字典
        all_occur = len(seq) - k + 1                  #k-mer序列的数量
        for x in range(all_occur):                    #计算每个k-mer序列的出现次数
            key = seq[x: x+k]                         
            # curr_dict[key] += 1
            # curr_dict.update({ key: curr_dict[key]+1 })
            # curr_dict.update({ key1: curr_dict[key]+0.5, key2: curr_dict[key]+0.5 })
            # curr_dict.update({ lambda N: keyN: curr_dict[key]+1/N })
            count_result = nucleic_acid_count(dna_seq=key)
            for char, value in count_result:
                pass
        
        # curr_dict = {k:v for k,v in curr_dict.items()}#字典储存在curr_dict 

    for key in curr_dict.keys():                      
        l = len(seq) - len(key) + 1
        curr_dict[key] = curr_dict[key] / l           #计算tf
    return(curr_dict)

def feat2str(int_float_dict):
    feat_str = ''
    for index, value in int_float_dict.items():
        feat_str += str(index) + ':' + str(value) + ' '  
    #feat_str = str(label) + '\t' + feat_str 
    return feat_str

folder_path = '/home/nusri/下载/dna_feat/dataset/merged_nonprobiotics'
label = 0
dirs = [path for path in os.listdir(folder_path) if path.endswith('txt')]
seq_list = []
filename_list =[]
for input_file in dirs[:3]:
# for input_file in dirs:
    with open(os.path.join(folder_path, input_file)) as seqs:
        dna_seq = seqs.read()
        seq_list.append(dna_seq)
        filename_list.append(input_file)


feat_str_list =[]
for seq, filename in tqdm(zip(seq_list,filename_list),total=len(seq_list)):
    feat_str = ''
    try:
        feat_dict = get_dict(seq)
    except Exception as e:
        print('ERROR: ', filename)
        raise('ERROR')
    #给kmer序列标号并创建字典
    index_ATCG_map = {ATCG: int(i) for i, ATCG in enumerate(feat_dict.keys())}
    #创建index和tf的字典
    feat_dict = {index_ATCG_map[ATCG]: value for i, (ATCG, value) in enumerate(feat_dict.items())}
    feat_str += feat2str(feat_dict)
    feat_str_list.append(feat_str)

with open('result_nonprobio.tsv', 'w') as f:
    for feat_str in feat_str_list:
        f.write(str(label) + '\t' + feat_str + '\n')