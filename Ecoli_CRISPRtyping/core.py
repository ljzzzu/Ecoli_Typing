"""Ecoli_CRISPRtyping

Usage:
    CRISPRtyping [--type=TYPE --path=PATH]

"""
import docopt
import sys
import os
# 添加自定义模块所在的路径
module_path = os.path.abspath(os.path.join(os.path.dirname(__file__)))
sys.path.append(module_path)
import file_transform
import CRISPR_extract
import pandas as pd
import CRISPR_type

def getpath(s=CRISPR_type.get_separator()):
    CRISPR_path=os.path.abspath(os.path.join(os.path.dirname(__file__)))+s+'spacer_database'
    order_path=os.path.abspath(os.path.join(os.path.dirname(__file__)))+s+'order_database'
    inpath = os.path.abspath(os.path.join(os.path.dirname(__file__))) + s+'test'
    return [CRISPR_path, order_path, inpath]

def write_updata(path,upd):
    ss=''
    for i in upd:
        for ii in i[0]:
            ss=ss+ii+','
        ss=ss+'\n'
    with open(path, 'a') as file:
        file.write(ss)

def write_updata_CCT(path,upd):
    ss=''
    for i in upd:
        ss=ss+i[0][0]+','+i[0][1]+'\n'
    with open(path, 'a') as file:
        file.write(ss)

def write_updata_spacer(path,upd):
    ss=''
    for i in upd:
        ss=ss+i+','+upd[i]+'\n'
    with open(path,'a') as file:
        file.write(ss)

def CCT_typing(inpath,n_max=500,CRISPR_path=getpath()[0],order_path=getpath()[1],output=True,s=CRISPR_type.get_separator()):
    dic1 = file_transform.getFileName(inpath)
    list_CCT=[]
    for ik in range(len(dic1) // n_max + 1):
        dict_fasta = file_transform.line_merge_block(inpath, dic1, n_start=1 + ik * n_max, n_max=n_max)
        if dict_fasta:
            CRISPR_sequence_list = CRISPR_extract.CRISPR_extract(dict_fasta, silence=False)
            CRISPR_seq_final = CRISPR_extract.transform_for_csv(CRISPR_sequence_list)
            CRISPR_dict = CRISPR_extract.select_seq(CRISPR_seq_final)
            data = CRISPR_extract.seq_to_csv(CRISPR_dict)
            data_filter = CRISPR_type.repeat_filter(data)
            data1 = pd.read_csv(CRISPR_path + s+'CRISPR1T.csv', header=None)
            data2 = pd.read_csv(CRISPR_path + s+'CRISPR1F.csv', header=None)
            data3 = pd.read_csv(CRISPR_path +s+ 'CRISPR2T.csv', header=None)
            data4 = pd.read_csv(CRISPR_path +s+ 'CRISPR2F.csv', header=None)
            data5 = pd.read_csv(order_path + s+'CRISPR1.csv', header=None)
            data6 = pd.read_csv(order_path + s+'CRISPR2.csv', header=None)
            data7 = pd.read_csv(order_path + s+'CCT.csv', header=None)
            CRISPR1_dic = file_transform.transform_for_dict(data1, data2)
            CRISPR2_dic = file_transform.transform_for_dict(data3, data4)
            order_dic1 = file_transform.transform_order_for_dict(data5)
            order_dic2 = file_transform.transform_order_for_dict(data6)
            order_final = file_transform.transform_final_for_dict(data7)
            data_all = CRISPR_type.CRISPR_filter_form(data_filter, CRISPR1_dic, CRISPR2_dic, dic1)
            data_CRISPR = CRISPR_type.confirm_CRISPR(data_all)
            data_order1 = CRISPR_type.CRISPR_order(data_CRISPR[0], order_dic1, 5)
            data_order2 = CRISPR_type.CRISPR_order(data_CRISPR[1], order_dic2, 6)
            data_CCT = CRISPR_type.confirm_CCT(data_order1, data_order2, order_final)
            list_CCT.append(data_CCT)
            if output:
                directory = "output"
                if not os.path.exists(directory):
                    os.mkdir(directory)
                try:
                    s_crispr1 = ''
                    for index, row in data_order1.iterrows():
                        if row[1]:
                            s_crispr1 = s_crispr1 + row[0] + ','
                            for i in row[1]:
                                s_crispr1 = s_crispr1 + str(i) + ','
                            s_crispr1 = s_crispr1 + '\n'
                    with open(inpath+s+'output'+s+'spacer1.csv', 'a') as file:
                        file.write(s_crispr1)
                    s_crispr2 = ''
                    for index, row in data_order2.iterrows():
                        if row[1]:
                            s_crispr2 = s_crispr2 + row[0] + ','
                            for i in row[1]:
                                s_crispr2 = s_crispr2 + str(i) + ','
                            s_crispr2 = s_crispr2 + '\n'
                    with open(inpath+s+'output'+s+'spacer2.csv', 'a') as file:
                        file.write(s_crispr2)
                except:
                    print('cannot save result')
    data_CCT=list_CCT[0]
    for i in range(1,len(list_CCT)):
        data_CCT=pd.concat([data_CCT,list_CCT[i]])
    if output:
        try:
            data_CCT.to_csv(inpath + s+'output'+s+'final_CCT.csv', index=False)
        except:
            print('cannot save result')
    return data_CCT

def CCT_typing_remove(inpath,n_max=500,CRISPR_path=getpath()[0],order_path=getpath()[1],output=True,s=CRISPR_type.get_separator()):
    dic1 = file_transform.getFileName(inpath)
    list_CCT=[]
    e_list = []
    for ik in range(len(dic1) // n_max + 1):
        dict_fasta = file_transform.line_merge_block(inpath, dic1, n_start=1 + ik * n_max, n_max=n_max)
        if dict_fasta:
            CRISPR_sequence_list = CRISPR_extract.CRISPR_extract(dict_fasta, silence=False)
            CRISPR_seq_final = CRISPR_extract.transform_for_csv(CRISPR_sequence_list)
            CRISPR_dict = CRISPR_extract.select_seq(CRISPR_seq_final)
            data = CRISPR_extract.seq_to_csv(CRISPR_dict)
            data_filter = CRISPR_type.repeat_filter(data)
            data1 = pd.read_csv(CRISPR_path + s+'CRISPR1T.csv', header=None)
            data2 = pd.read_csv(CRISPR_path + s+'CRISPR1F.csv', header=None)
            data3 = pd.read_csv(CRISPR_path +s+ 'CRISPR2T.csv', header=None)
            data4 = pd.read_csv(CRISPR_path +s+ 'CRISPR2F.csv', header=None)
            data5 = pd.read_csv(order_path + s+'CRISPR1.csv', header=None)
            data6 = pd.read_csv(order_path + s+'CRISPR2.csv', header=None)
            data7 = pd.read_csv(order_path + s+'CCT.csv', header=None)
            CRISPR1_dic = file_transform.transform_for_dict(data1, data2)
            CRISPR2_dic = file_transform.transform_for_dict(data3, data4)
            order_dic1 = file_transform.transform_order_for_dict(data5)
            order_dic2 = file_transform.transform_order_for_dict(data6)
            order_final = file_transform.transform_final_for_dict(data7)
            data_all = CRISPR_type.CRISPR_filter_form(data_filter, CRISPR1_dic, CRISPR2_dic, dic1)
            data_CRISPR = CRISPR_type.confirm_CRISPR(data_all)
            e_list.append(CRISPR_type.remove_low_quality_data(data_CRISPR[0]))
            e_list.append(CRISPR_type.remove_low_quality_data(data_CRISPR[1]))
            data_order1 = CRISPR_type.CRISPR_order(data_CRISPR[0], order_dic1, 5)
            data_order2 = CRISPR_type.CRISPR_order(data_CRISPR[1], order_dic2, 6)
            data_CCT = CRISPR_type.confirm_CCT(data_order1, data_order2, order_final)
            list_CCT.append(data_CCT)
            if output:
                directory = "output"
                if not os.path.exists(directory):
                    os.mkdir(directory)
                try:
                    s_crispr1 = ''
                    for index, row in data_order1.iterrows():
                        if row[1]:
                            s_crispr1 = s_crispr1 + row[0] + ','
                            for i in row[1]:
                                s_crispr1 = s_crispr1 + str(i) + ','
                            s_crispr1 = s_crispr1 + '\n'
                    with open(inpath+s+'output'+s+'spacer1.csv', 'a') as file:
                        file.write(s_crispr1)
                    s_crispr2 = ''
                    for index, row in data_order2.iterrows():
                        if row[1]:
                            s_crispr2 = s_crispr2 + row[0] + ','
                            for i in row[1]:
                                s_crispr2 = s_crispr2 + str(i) + ','
                            s_crispr2 = s_crispr2 + '\n'
                    with open(inpath+s+'output'+s+'spacer2.csv', 'a') as file:
                        file.write(s_crispr2)
                except:
                    print('cannot save result')
    data_CCT=list_CCT[0]
    for i in range(1,len(list_CCT)):
        data_CCT=pd.concat([data_CCT,list_CCT[i]])
    ss=''
    for i in e_list:
        for j in i:
            ss=ss+j+'\n'
    if output:
        try:
            data_CCT.to_csv(inpath + s+'output'+s+'final_CCT.csv', index=False)
            with open(inpath + s + 'output' + s + 'lowq.txt', 'w') as file:
                file.write(ss)
        except:
            print('cannot save result')
    return [data_CCT,e_list]

def CCT_typing_update(inpath,n_max=500,CRISPR_path=getpath()[0],order_path=getpath()[1],output=True,new_spacer=True,new_CCT=True,s=CRISPR_type.get_separator()):
    dic1 = file_transform.getFileName(inpath)
    list_CCT=[]
    for ik in range(len(dic1) // n_max + 1):
        dict_fasta = file_transform.line_merge_block(inpath, dic1, n_start=1 + ik * n_max, n_max=n_max)
        if dict_fasta:
            CRISPR_sequence_list = CRISPR_extract.CRISPR_extract(dict_fasta, silence=False)
            CRISPR_seq_final = CRISPR_extract.transform_for_csv(CRISPR_sequence_list)
            CRISPR_dict = CRISPR_extract.select_seq(CRISPR_seq_final)
            data = CRISPR_extract.seq_to_csv(CRISPR_dict)
            data_filter = CRISPR_type.repeat_filter(data)
            data1 = pd.read_csv(CRISPR_path + s+'CRISPR1T.csv', header=None)
            data2 = pd.read_csv(CRISPR_path + s+'CRISPR1F.csv', header=None)
            data3 = pd.read_csv(CRISPR_path +s+ 'CRISPR2T.csv', header=None)
            data4 = pd.read_csv(CRISPR_path +s+ 'CRISPR2F.csv', header=None)
            data5 = pd.read_csv(order_path + s+'CRISPR1.csv', header=None)
            data6 = pd.read_csv(order_path + s+'CRISPR2.csv', header=None)
            data7 = pd.read_csv(order_path + s+'CCT.csv', header=None)
            CRISPR1_dic = file_transform.transform_for_dict(data1, data2)
            CRISPR2_dic = file_transform.transform_for_dict(data3, data4)
            order_dic1 = file_transform.transform_order_for_dict(data5)
            order_dic2 = file_transform.transform_order_for_dict(data6)
            order_final = file_transform.transform_final_for_dict(data7)
            filter_all = CRISPR_type.CRISPR_filter_form(data_filter, CRISPR1_dic, CRISPR2_dic, dic1)
            CRISPR1_n = file_transform.CRISPR_identifier(data1, 1)
            CRISPR2_n = file_transform.CRISPR_identifier(data3, 1)
            data_CRISPR = CRISPR_type.confirm_CRISPR_updata(filter_all, CRISPR1_n, CRISPR2_n)
            if new_spacer:
                if data_CRISPR[2] != {}:
                    t1 = {}
                    for i in data_CRISPR[2]:
                        t1[file_transform.complementary(i)] = data_CRISPR[2][i]
                    write_updata_spacer(getpath()[0] + s + 'CRISPR1T.csv', data_CRISPR[2])
                    write_updata_spacer(getpath()[0] + s + 'CRISPR1F.csv', t1)
                if data_CRISPR[3] != {}:
                    t2 = {}
                    for i in data_CRISPR[3]:
                        t2[file_transform.complementary(i)] = data_CRISPR[3][i]
                    write_updata_spacer(getpath()[1] + s + 'CRISPR2T.csv', data_CRISPR[3])
                    write_updata_spacer(getpath()[0] + s + 'CRISPR2F.csv', t2)
            data_order1 = CRISPR_type.CRISPR_order_updata(data_CRISPR[0], order_dic1, 5)
            data_order2 = CRISPR_type.CRISPR_order_updata(data_CRISPR[1], order_dic2, 6)
            data_CCT = CRISPR_type.confirm_CCT_updata(data_order1[0], data_order2[0], order_final)
            if new_CCT:
                try:
                    write_updata(getpath()[1] + s + 'CRISPR1.csv', data_order1[1])
                    write_updata(getpath()[1] + s + 'CRISPR2.csv', data_order2[1])
                    write_updata_CCT(getpath()[1] + s + 'CCT.csv', data_CCT[1])
                except:
                    print('can not update')
            list_CCT.append(data_CCT)
            if output:
                directory = "output"
                if not os.path.exists(directory):
                    os.mkdir(directory)
                try:
                    s_crispr1 = ''
                    for index, row in data_order1[0].iterrows():
                        if row[1]:
                            s_crispr1 = s_crispr1 + row[0] + ','
                            for i in row[1]:
                                s_crispr1 = s_crispr1 + str(i) + ','
                            s_crispr1 = s_crispr1 + '\n'
                    with open(inpath+s+'output'+s+'spacer1.csv', 'a') as file:
                        file.write(s_crispr1)
                    s_crispr2 = ''
                    for index, row in data_order2[0].iterrows():
                        if row[1]:
                            s_crispr2 = s_crispr2 + row[0] + ','
                            for i in row[1]:
                                s_crispr2 = s_crispr2 + str(i) + ','
                            s_crispr2 = s_crispr2 + '\n'
                    with open(inpath+s+'output'+s+'spacer2.csv', 'a') as file:
                        file.write(s_crispr2)
                except:
                    print('cannot save result')
    data_CCT=list_CCT[0][0]
    for i in range(1,len(list_CCT)):
        data_CCT=pd.concat([data_CCT,list_CCT[i][0]])
    if output:
        try:
            data_CCT.to_csv(inpath + s+'output'+s+'final_CCT.csv', index=False)
        except:
            print('cannot save result')
    return data_CCT

def main(s=CRISPR_type.get_separator()):
    args = docopt.docopt(__doc__)
    type_method = args['--type'] or 'normal'
    in_path=args['--path']
    if in_path:
        nowpath=in_path
    else:
        nowpath=os.path.dirname(os.path.abspath(__file__))
    print(nowpath)
    if type_method=='normal':
        data_CCT=CCT_typing(inpath=nowpath)
    elif type_method == 'test':
        data_CCT = CCT_typing(inpath=getpath()[2])
    elif type_method == 'updata':
        data_CCT = CCT_typing_update(inpath=nowpath)
    elif type_method == 'high':
        data_CCT =CCT_typing_remove(nowpath)
