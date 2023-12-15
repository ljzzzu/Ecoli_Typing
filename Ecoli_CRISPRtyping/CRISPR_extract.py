import re,os
import pandas as pd

#提取每个fasta文档中的CRISPR序列
def CRISPR_extract(dic, silence=False):
    pattern1 = '\w[G,A,T]\w[T,G,A][T,G][T,A][A,C,T][T,C][C,A,T]{3}[C,G,T]\w[C,A,T][T,A,G][C,G,A][G,A]\w{2}[C,G,T][A,G][G,T][G,A,T][G,A][A,T]\w[C,A,T]{3}'
    pattern2 = '[G,T][T,A,C]\w{2}[A,T][C,T][C,A,T][C,A][T,C][C,G,A]\w{2}[C,G,T][C,G,T][A,T,C][G,A,T]\w[C,G,A][G,A,T]{3}[A,G][T,G,A][A,T][A,C][A,C,T]\w[C,A,T]\w'
    CRISPR_sequence_list = []
    CRISPR_list = []
    CRISPR_dic = {}
    flag = 0
    for x in dic:
        CRISPR_dic[x] = []
        line = dic[x]
        for match in re.finditer(pattern1, line):
            repeat_s = match.start()
            repeat_e = match.end()
            list2 = [repeat_s, repeat_e]
            CRISPR_dic[x].append(list2)
        for match_1 in re.finditer(pattern2, line):
            repeat_s = match_1.start()
            repeat_e = match_1.end()
            list2 = [repeat_s, repeat_e]
            CRISPR_dic[x].append(list2)
        if len(CRISPR_dic[x]) <= 1:
            del CRISPR_dic[x]
        else:
            for key, value in CRISPR_dic.items():
                CRISPR_dic = {}
                i = 0
                CRISPR_sequence = []
                if flag != x.split('>')[1]:
                    flag = x.split('>')[1]
                    CRISPR_sequence_list.append(CRISPR_list)
                    if not silence:
                        print("CRISPR from file" + str(int(flag) - 1) + " is extracted")
                    CRISPR_list = []
                CRISPR_sequence.append(key)
                while i < len(value):
                    repeat_s1 = value[i][0]
                    repeat_e1 = value[i][1]
                    repeat_sequence = line[repeat_s1:repeat_e1]
                    a = len(repeat_sequence)
                    CRISPR_sequence.append([repeat_s1, repeat_e1])
                    CRISPR_sequence.append(repeat_sequence)
                    if i != len(value) - 1:
                        repeat_s2 = value[i + 1][0]
                        if repeat_s2 - repeat_s1 > 200:
                            CRISPR_sequence.append("")
                            CRISPR_sequence.append("")
                            CRISPR_sequence.append("CRISPR")
                            CRISPR_sequence.append("")
                            CRISPR_sequence.append("")
                            CRISPR_sequence.append("")
                        else:
                            spacer_s1 = value[i][1]
                            spacer_e1 = value[i + 1][0]
                            spacer_sequence = line[spacer_s1:spacer_e1]
                            CRISPR_sequence.append(spacer_sequence)
                            b = len(spacer_sequence)
                            CRISPR_sequence.append([a, b])
                    i = i + 1
            CRISPR_list.append(CRISPR_sequence)
    return CRISPR_sequence_list


#将CRISPR_sequence_list 转换为符合csv文档的格式
def transform_for_csv(CRISPR_sequence_list1):
    CRISPR_sequence_final = []
    for CRISPR_fasta_each in CRISPR_sequence_list1:
        for CRISPR_each in CRISPR_fasta_each:
            list3 = ["","","","",""]
            CRISPR_sequence_final.append(list3)
            CRISPR_each.append("")
            CRISPR_each.append("")
            i = -1
            for CRISPR_each_element in CRISPR_each:
                i = i + 1
                if CRISPR_each_element:
                    if i == 0:
                        list1 = []
                        list1.append(CRISPR_each_element)
                        list1.append("")
                        list1.append("")
                        list1.append("")
                        list1.append("")
                        CRISPR_sequence_final.append(list1)
                    elif i % 4 == 1:
                        list2 = []
                        list2.append("")
                        list2.append(CRISPR_each_element)
                    elif i % 4 == 2:
                        list2.append(CRISPR_each_element)
                    elif i % 4 == 3:
                        list2.append(CRISPR_each_element)
                    elif i % 4 == 0:
                        list2.append(CRISPR_each_element)
                        CRISPR_sequence_final.append(list2)
                else:
                    if i % 4 == 2:
                        list2.append("")
                    elif i % 4 == 3:
                        list2.append("")
                    if i % 4 == 0:
                        list2.append("")
                        CRISPR_sequence_final.append(list2)
    return CRISPR_sequence_final
# 过滤掉CRISPR中单个repeat的情况
def select_seq(CRISPR_sequence_final):
    CRISPR_dic = {}
    i = -1
    for line in CRISPR_sequence_final:
        i = i + 1
        if i >=1:
            if line[0] != "":
                if line[0][0] == ">":
                    x = line[0].strip()
                    CRISPR_dic[x] = []
            else:
                if line[1] == "CRISPR":
                    if CRISPR_sequence_final[i + 1][3] != "":
                        CRISPR_dic[x].append(CRISPR_sequence_final[i])
                elif line[1] == "":
                    continue
                else:
                    if line[3] != "":
                        if len(line[3]) <= 40:
                            CRISPR_dic[x].append(line)
#将长度大于85bp且小于95bp的序列截断
                        elif len(line[3]) >= 85 and len(line[3]) <= 95:
                            spacer_1 = line[3][:32]
                            repeat = line[3][32:61]
                            spacer_2 = line[3][61:]
                            line[3] = spacer_1
                            line[4] = [29,32]
                            CRISPR_dic[x].append(line)
                            new_line = []
                            new_line.append("")
                            new_line.append("new_line")
                            new_line.append(repeat)
                            new_line.append(spacer_2)
                            spacer_length = len(spacer_2)
                            new_line.append([29, spacer_length])
                            CRISPR_dic[x].append(new_line)
#将长度大于96bp且小于160bp的序列截断
                        elif len(line[3]) > 96 and len(line[3]) <= 160:
                            spacer_1 = line[3][:32]
                            repeat_1 = line[3][32:61]
                            spacer_2 = line[3][61:93]
                            repeat_2 = line[3][93:122]
                            spacer_3 = line[3][122:]
                            line[3] = spacer_1
                            line[4] = [29,32]
                            CRISPR_dic[x].append(line)
                            new_line = []
                            new_line.append("")
                            new_line.append("new_line")
                            new_line.append(repeat_1)
                            new_line.append(spacer_2)
                            spacer_length = len(spacer_2)
                            new_line.append([29, spacer_length])
                            CRISPR_dic[x].append(new_line)
                            new_line = []
                            new_line.append("")
                            new_line.append("new_line")
                            new_line.append(repeat_2)
                            new_line.append(spacer_3)
                            spacer_length = len(spacer_3)
                            new_line.append([29, spacer_length])
                            CRISPR_dic[x].append(new_line)
                    elif CRISPR_sequence_final[i - 1][3] != "":
                        CRISPR_dic[x].append(CRISPR_sequence_final[i])
    return CRISPR_dic
#
def seq_to_csv(CRISPR_dic):
    col=['name','start_end','repeats','spacers','length','6']
    CRISPR_sequence_final = []
    for key,value in CRISPR_dic.items():
        if CRISPR_dic[key] != []:
            list1 = []
            list1.append(key)
            list1.append("")
            list1.append("")
            list1.append("")
            list1.append("")
            CRISPR_sequence_final.append(list1)
            for i in value:
                CRISPR_sequence_final.append(i)
    data2 = pd.DataFrame(data=CRISPR_sequence_final, columns=None, index=None)
    return data2



















