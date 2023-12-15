# -*- coding:utf-8 -*-
import os
import pandas as pd

#调整输入文件格式并用字典储存
def line_merge(in_path,silence=True):
    os.chdir(in_path)
    flag=0
    dic={}
    for item in os.listdir():
        if item.endswith(".fna") or item.endswith(".fasta"):
            flag=flag+1
            with open(item,'r') as f:
                for line in f.readlines():
                    if line.startswith(">"):
                        key ='>'+str(flag)+line.strip()
                        dic[key] = ''
                    else:
                        dic[key]=dic[key]+line.strip()
            if not silence:
                print(item+" finish")
    dic['>' + str(flag + 1) + '>'] = 'GAGTTCCCCGCGTCAGCGGGGATAAACCGGAGTTCCCCGCGTCAGCGGGGATAAACCG'
    return dic

#获取fna与fasta结尾的文件名并用列表储存
def getFileName(in_path):
    os.chdir(in_path)
    flag = 0
    dic1 = {}
    for item in os.listdir():
        if item.endswith(".fna") or item.endswith(".fasta"):
            flag = flag + 1
            dic1[str(flag)]=item
    return dic1

#分块调整输入文件格式并用列表储存
def line_merge_block(in_path,dic1,n_start=1,n_max=500,silence=True):
    os.chdir(in_path)
    dic = {}
    for flag in range(n_start, n_start + n_max):
        if str(flag) in dic1:
            with open(dic1[str(flag)], 'r') as f:
                for line in f.readlines():
                    if line.startswith(">"):
                        key = '>' + str(flag) + line.strip()
                        dic[key] = ''
                    else:
                        dic[key] = dic[key] + line.strip()
            if not silence:
                print(dic1[str(flag)] + " finish")
    dic['>' + str(n_start + n_max + 1) + '>'] = 'GAGTTCCCCGCGTCAGCGGGGATAAACCGGAGTTCCCCGCGTCAGCGGGGATAAACCG'
    return dic

#将excel文件整理成字典的形式
def transform_for_dict(excel,excell):
    excel_dic={}
    excel1=excel.dropna()
    excell1=excell.dropna()
    for index, row in excel1.iterrows():
        excel_dic[row[0]]=row[1]
    for index, row in excell1.iterrows():
        excel_dic[row[0]]=row[1]
    return excel_dic

#将排列方式转换为字典
def transform_order_for_dict(excel):
    order={}
    for index,row in excel.iterrows():
        lt=[]
        for i in range(len(row)):
            if pd.isna(row[i]):
                continue
            else:
                lt.append(str(int(row[i])))
        order[str(lt)]=int(index+1)
    return order

#将CCT转换为字典
def transform_final_for_dict(excel):
    order={}
    for index,row in excel.iterrows():
        lt=[]
        for i in range(len(row)):
            if pd.isna(row[i]):
                lt.append('')
            else:
                lt.append(str(int(row[i])))
        order[str(lt)]=int(index+1)
    return order

#获取间隔序列编号的最大值
def CRISPR_identifier(df,col):
    df1=df.dropna()
    df1[2]=df1[col].str.split('_').str[-1]
    df1[2]=df1[2].astype(int)
    identifier=df1.sort_values(by=2, ascending=False).iloc[0,2]
    return identifier
#碱基互补配对
def complementary(s):
    s1=''
    for i in s:
        if i=='C':
            s1=s1+ 'G'
        elif i=='G':
            s1=s1+ 'C'
        elif i=='A':
            s1=s1+ 'T'
        elif i=='T':
            s1=s1+ 'A'
    return s1[::-1]

#独热编码
def one_hot(filename):
    dict_all={}
    with open(filename, 'r') as file:
        lines = file.readlines()
    for line in lines:
        t=line[:-1].split(',')
        dict_all[t[0]]=[]
        for i in range(1,len(t)):
            if t[i]!='':
                dict_all[t[0]].append(t[i])
    sp_all=[]
    for i in dict_all.values():
        sp_all=sp_all+i
    sp_all=list(set(sp_all))
    l_all=[]
    for j in dict_all.keys():
        lt=[j]
        for i in sp_all:
            if i in dict_all[j]:
                lt.append(1)
            else:
                lt.append(0)
        l_all.append(lt)
    col=['ID']+sp_all
    result= pd.DataFrame(l_all, columns=col)
    return result

#删除低质量序列
def delete_ID(filename,e_list):
    s=''
    with open(filename, 'r') as file:
        lines = file.readlines()
    for line in lines:
        t = line[:-1].split(',')[0]
        if t not in e_list:
            s=s+line
    with open(filename, 'w') as file:
        file.write(s)