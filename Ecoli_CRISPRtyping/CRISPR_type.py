import pandas as pd
import os
import numpy as np
import re
import platform
#获取不同系统路径的分隔符
def get_separator():
    s='/'
    if platform.system()=='Windows':
        s='\\'
    return s
#过滤CRISPR提取的结果
def repeat_filter(df,output=False,filename='final_result_CRISPR.csv',s=get_separator()):
    data_filter=df.copy()
    drop_index=[]
    for index,row in data_filter.iterrows():
        if row[2]!='':
            if row[4]!='':
                if row[4] not in [[29,32],[29,33]]:
                    drop_index.append(index)
    data_filter=data_filter.drop(drop_index)
    drop_index=[]
    dic_filter={}
    dic_n={}
    for index,row in data_filter.iterrows():
        if row[2]=='':
            if np.array(list(dic_n.values())).sum()==len(dic_n):
                for x in dic_filter:
                    drop_index=drop_index+list(dic_filter[x])
            dic_filter={}
            dic_n={}
        else:
            if row[2] not in dic_filter:
                dic_filter[row[2]]=[]
                dic_n[row[2]]=0
            dic_filter[row[2]].append(index)
            dic_n[row[2]]=dic_n[row[2]]+1
    if np.array(list(dic_n.values())).sum()==len(dic_n):
        for x in dic_filter:
            drop_index=drop_index+list(dic_filter[x])
    data_filter=data_filter.drop(drop_index)
    if output:
        directory = "output"
        if not os.path.exists(directory):
            os.mkdir(directory)
        try:
            data_filter.to_csv(s+'output'+s+filename,index=False)
        except:
            print('cannot save result')
    return data_filter

#整理过滤后CRISPR的文件形式
def CRISPR_filter_form(data2,CRISPR1_dic,CRISPR2_dic,dic1):
    datat=data2.copy()
    datat[5]=datat[3]
    datat[6]=datat[3]
    datat[5]=datat[5].map(CRISPR1_dic)
    datat[6]=datat[6].map(CRISPR2_dic)
    data_CRISPR=datat.drop(datat[datat.iloc[:, 3] == ""].index).reset_index()
    datat[7]=datat[0].str.split('>').str[1]
    datat[7]=datat[7].map(dic1)
    data_name=datat.dropna(subset=[datat.columns[7]])
    data_name = data_name.drop_duplicates(subset=[data_name.columns[7]]).reset_index()
    data_all=pd.concat([data_name,data_CRISPR])
    data_all[8]=data_all['index']
    data_all=data_all.sort_values(by='index', ascending=True)
    data_all=data_all.drop_duplicates(subset=[data_all.columns[9]])
    data_all=data_all[['index',7,1,2,3,4,5,6]]
    return data_all

#将CRISPR1和CRISPR2分开，并将新出现的间隔序列编号
def confirm_CRISPR_updata(data_all,CRISPR1_n,CRISPR2_n):
    data_all=data_all.set_index('index')
    data_both=data_all.dropna(subset=[data_all.columns[5],data_all.columns[6]])
    data_CRISPR1=data_all.dropna(subset=[data_all.columns[5]])
    data_CRISPR2=data_all.dropna(subset=[data_all.columns[6]])
    list_none=list(set(list(data_CRISPR1.index)+list(data_CRISPR2.index)))
    data_none=data_all[data_all[3]!='']
    data_none=data_none.drop(list_none)
    data_name=data_all.dropna(subset=[data_all.columns[0]])
    dict_type={}
    list_none=list(data_none.index)
    list_both=list(data_both.index)
    list_CRISPR1=list(data_CRISPR1.index)
    list_CRISPR2=list(data_CRISPR2.index)
    for i in list_CRISPR1:
        dict_type[i]=1
    for i in list_CRISPR2:
        dict_type[i]=2
    for i in list_both:
        dict_type[i]=1.5
    for i in list_none:
        dict_type[i]=1.5
    for i in list_both+list_none:
        if i-1 not in dict_type and i+1 not in dict_type:
            dict_type[i]=1.5
        elif i-1 in dict_type and i+1 in dict_type:
            if dict_type[i-1]==dict_type[i+1]:
                dict_type[i]=dict_type[i-1]
            else:
                CRISPR_sum=0
                CRISPR_n=0
                flag=i-1
                while flag in dict_type:
                    CRISPR_n=CRISPR_n+1
                    CRISPR_sum=CRISPR_sum+dict_type[flag]
                    flag=flag-1
                flag=i+1
                while flag in dict_type:
                    CRISPR_n=CRISPR_n+1
                    CRISPR_sum=CRISPR_sum+dict_type[flag]
                    flag=flag+1
                if CRISPR_sum/CRISPR_n<1.5:
                    dict_type[i]=1
                else:
                    dict_type[i]=2
        else:
            CRISPR_sum=0
            CRISPR_n=0
            flag=i-1
            while flag in dict_type:
                CRISPR_n=CRISPR_n+1
                CRISPR_sum=CRISPR_sum+dict_type[flag]
                flag=flag-1
            flag=i+1
            while flag in dict_type:
                CRISPR_n=CRISPR_n+1
                CRISPR_sum=CRISPR_sum+dict_type[flag]
                flag=flag+1
            if CRISPR_sum/CRISPR_n<1.5:
                dict_type[i]=1
            else:
                dict_type[i]=2
    for i in list_both:
        if dict_type[i]==1:
            data_CRISPR2=data_CRISPR2.drop(i)
        elif dict_type[i]==2:
            data_CRISPR1=data_CRISPR1.drop(i)
        else:
            data_CRISPR1=data_CRISPR1.drop(i)
            data_CRISPR2=data_CRISPR2.drop(i)
    t1={}
    t2={}
    for i in list_none:
        if dict_type[i]==1:
            if data_none.loc[i,3] not in t1:
                t1[data_none.loc[i,3]]='>CRISPR1_'+str(CRISPR1_n)
                CRISPR1_n=CRISPR1_n+1
            data_none.loc[i,5]=t1[data_none.loc[i,3]]
            data_CRISPR1 = pd.concat([data_CRISPR1, data_none.loc[i].to_frame().T], ignore_index=True)
        elif dict_type[i]==2:
            if data_none.loc[i,3] not in t2:
                t2[data_none.loc[i,3]]='>CRISPR2_'+str(CRISPR2_n)
                CRISPR2_n=CRISPR2_n+1
            data_none.loc[i,6]=t2[data_none.loc[i,3]]
            data_CRISPR2 = pd.concat([data_CRISPR2, data_none.loc[i].to_frame().T], ignore_index=True)
    data_CRISPR1=pd.concat([data_name.reset_index(),data_CRISPR1.reset_index()]).sort_values(by='index', ascending=True)
    data_CRISPR2=pd.concat([data_name.reset_index(),data_CRISPR2.reset_index()]).sort_values(by='index', ascending=True)
    data_CRISPR1=data_CRISPR1.drop(6,axis=1)
    data_CRISPR2=data_CRISPR2.drop(5,axis=1)
    return [data_CRISPR1,data_CRISPR2,t1,t2]



#整理提取CRISPR后的文件形式
def CRISPR_extract_form(data2,CRISPR1_dic,CRISPR2_dic,dic1):
    datat=data2.copy()
    datat[5]=datat[3]
    datat[6]=datat[3]
    datat[5]=datat[5].map(CRISPR1_dic)
    datat[6]=datat[6].map(CRISPR2_dic)
    data_CRISPR1=datat.dropna(subset=[datat.columns[5]]).reset_index()
    data_CRISPR2=datat.dropna(subset=[datat.columns[6]]).reset_index()
    datat[7]=datat[0].str.split('>').str[1]
    datat[7]=datat[7].map(dic1)
    data_name=datat.dropna(subset=[datat.columns[7]])
    data_name = data_name.drop_duplicates(subset=[data_name.columns[7]]).reset_index()
    data_all=pd.concat([data_name,data_CRISPR1,data_CRISPR2])
    data_all[8]=data_all['index']
    data_all=data_all.sort_values(by='index', ascending=True)
    data_all=data_all.drop_duplicates(subset=[data_all.columns[9]])
    data_all=data_all[['index',7,1,2,3,4,5,6]]
    return data_all

#将CRISPR1和CRISPR2分开
def confirm_CRISPR(data_all):
    data_all=data_all.set_index('index')
    data_both=data_all.dropna(subset=[data_all.columns[5],data_all.columns[6]])
    data_CRISPR1=data_all.dropna(subset=[data_all.columns[5]])
    data_CRISPR2=data_all.dropna(subset=[data_all.columns[6]])
    data_name=data_all.dropna(subset=[data_all.columns[0]])
    dict_type={}
    list_both=list(data_both.index)
    list_CRISPR1=list(data_CRISPR1.index)
    list_CRISPR2=list(data_CRISPR2.index)
    for i in list_CRISPR1:
        dict_type[i]=1
    for i in list_CRISPR2:
        dict_type[i]=2
    for i in list_both:
        dict_type[i]=1.5
    for i in list_both:
        if i-1 not in dict_type and i+1 not in dict_type:
            dict_type[i]=1.5
        elif i-1 in dict_type and i+1 in dict_type:
            if dict_type[i-1]==dict_type[i+1]:
                dict_type[i]=dict_type[i-1]
            else:
                CRISPR_sum=0
                CRISPR_n=0
                flag=i-1
                while flag in dict_type:
                    CRISPR_n=CRISPR_n+1
                    CRISPR_sum=CRISPR_sum+dict_type[flag]
                    flag=flag-1
                flag=i+1
                while flag in dict_type:
                    CRISPR_n=CRISPR_n+1
                    CRISPR_sum=CRISPR_sum+dict_type[flag]
                    flag=flag+1
                if CRISPR_sum/CRISPR_n<1.5:
                    dict_type[i]=1
                else:
                    dict_type[i]=2
        else:
            CRISPR_sum=0
            CRISPR_n=0
            flag=i-1
            while flag in dict_type:
                CRISPR_n=CRISPR_n+1
                CRISPR_sum=CRISPR_sum+dict_type[flag]
                flag=flag-1
            flag=i+1
            while flag in dict_type:
                CRISPR_n=CRISPR_n+1
                CRISPR_sum=CRISPR_sum+dict_type[flag]
                flag=flag+1
            if CRISPR_sum/CRISPR_n<1.5:
                dict_type[i]=1
            else:
                dict_type[i]=2
    for i in list_both:
        if dict_type[i]==1:
            data_CRISPR2=data_CRISPR2.drop(i)
        elif dict_type[i]==2:
            data_CRISPR1=data_CRISPR1.drop(i)
        else:
            data_CRISPR1=data_CRISPR1.drop(i)
            data_CRISPR2=data_CRISPR2.drop(i)
    data_CRISPR1=pd.concat([data_name.reset_index(),data_CRISPR1.reset_index()]).sort_values(by='index', ascending=True)
    data_CRISPR2=pd.concat([data_name.reset_index(),data_CRISPR2.reset_index()]).sort_values(by='index', ascending=True)
    data_CRISPR1=data_CRISPR1.drop(6,axis=1)
    data_CRISPR2=data_CRISPR2.drop(5,axis=1)
    return [data_CRISPR1,data_CRISPR2]

#获取CRISPR系统的排列方式，不存在的赋予新的编号
def CRISPR_order_updata(df_CRISPR, order_dic, col, output=False, filename='order.csv',s=get_separator()):
    l_CRISPR1 = []
    lt = ['', [], 0]
    flag = 0
    n = len(order_dic)
    pattern1 = '\w[G,A,T]\w[T,G,A][T,G][T,A][A,C,T][T,C][C,A,T]{3}[C,G,T]\w[C,A,T][T,A,G][C,G,A][G,A]\w{2}[C,G,T][A,G][G,T][G,A,T][G,A][A,T]\w[C,A,T]{3}'
    for index, row in df_CRISPR.iterrows():
        if not pd.isna(row[7]):
            l_CRISPR1.append(lt)
            flag = 0
            lt = [row[7], [], 0]
        else:
            lt[1].append(row[col].split('_')[-1])
            matches = re.search(pattern1, row[2])
            if matches:
                lt[2] = 1
    l_CRISPR1.append(lt)
    for i in l_CRISPR1:
        if i[1] == []:
            l_CRISPR1.remove(i)
    for i in range(len(l_CRISPR1)):
        if l_CRISPR1[i][2] == 1:
            l_CRISPR1[i][1] = l_CRISPR1[i][1][::-1]
    data_order = pd.DataFrame(l_CRISPR1, columns=['name', 'order', 's'])
    data_order['type'] = data_order['order'].astype(str)
    data_order['type'] = data_order['type'].map(order_dic)
    updata = []
    for index, row in data_order.iterrows():
        if pd.isna(row[-1]):
            if row[1] != []:
                n = n + 1
                updata.append([row[1], n])
                data_order.iloc[index, -1] = n

    if output:
        directory = "output"
        if not os.path.exists(directory):
            os.mkdir(directory)
        try:
            data_order.to_csv(s+'output' +s+ filename, index=False)
        except:
            print('cannot save result')
    return [data_order, updata]

#获取CRISPR系统的排列方式
def CRISPR_order(df_CRISPR,order_dic,col,output=False,filename='order.csv',s=get_separator()):
    l_CRISPR1=[]
    lt=['',[],0]
    flag=0
    pattern1 = '\w[G,A,T]\w[T,G,A][T,G][T,A][A,C,T][T,C][C,A,T]{3}[C,G,T]\w[C,A,T][T,A,G][C,G,A][G,A]\w{2}[C,G,T][A,G][G,T][G,A,T][G,A][A,T]\w[C,A,T]{3}'
    for index,row in df_CRISPR.iterrows():
        if not pd.isna(row[7]):
            l_CRISPR1.append(lt)
            flag=0
            lt=[row[7],[],0]
        else:
            lt[1].append(row[col].split('_')[-1])
            matches = re.search(pattern1,row[2])
            if matches:
                lt[2]=1
    l_CRISPR1.append(lt)
    for i in l_CRISPR1:
        if i[1]==[]:
            l_CRISPR1.remove(i)
    for i in range(len(l_CRISPR1)):
        if l_CRISPR1[i][2]==1:
            l_CRISPR1[i][1]=l_CRISPR1[i][1][::-1]
    data_order=pd.DataFrame(l_CRISPR1,columns=['name','order','s'])
    data_order['type']=data_order['order'].astype(str)
    data_order['type']=data_order['type'].map(order_dic)
    if output:
        directory = "output"
        if not os.path.exists(directory):
            os.mkdir(directory)
        try:
            data_order.to_csv(s+'output'+s+filename,index=False)
        except:
            print('cannot save result')
    return data_order
#获取CCT分型
def confirm_CCT(order1,order2,order_final,output=False,filename='CCT.csv',s=get_separator()):
    order1 = order1.drop('s', axis=1)
    order2 = order2.drop('s', axis=1)
    order1.rename(columns={'order': 'order1','type':'type1'}, inplace=True)
    order2.rename(columns={'order': 'order2','type':'type2'}, inplace=True)
    data_order=pd.merge(order1,order2,on='name',how='outer')
    data_order=data_order.fillna(0)
    list_CCT=[]
    for index,row in data_order.iterrows():
        list_CCT.append([row[0], [str(int(row[2])), str(int(row[4]))]])
    CCT=pd.DataFrame(list_CCT,columns=['name','type'])
    CCT['CCT']=CCT['type'].astype(str)
    CCT['CCT']=CCT['CCT'].map(order_final)
    if output:
        directory = "output"
        if not os.path.exists(directory):
            os.mkdir(directory)
        try:
            CCT.to_csv(s+'output'+s+filename,index=False)
        except:
            print('cannot save result')
    return CCT

#获取CCT分型，并将新出现的赋予新的编号
def confirm_CCT_updata(order1,order2,order_final,output=False,filename='CCT.csv',s=get_separator()):
    n=len(order_final)
    order1=order1.drop('s',axis=1)
    order2=order2.drop('s',axis=1)
    order1.rename(columns={'order': 'order1','type':'type1'}, inplace=True)
    order2.rename(columns={'order': 'order2','type':'type2'}, inplace=True)
    data_order=pd.merge(order1,order2,on='name',how='outer')
    data_order=data_order.fillna(0)
    list_CCT=[]
    for index,row in data_order.iterrows():
        list_CCT.append([row[0],[str(int(row[2])),str(int(row[4]))]])
    CCT=pd.DataFrame(list_CCT,columns=['name','type'])
    CCT['CCT']=CCT['type'].astype(str)
    CCT['CCT']=CCT['CCT'].map(order_final)
    updata=[]
    for index,row in CCT.iterrows():
        if pd.isna(row[-1]):
            n=n+1
            updata.append([row[1],n])
            CCT.iloc[index,-1]=n
    if output:
        directory = "output"
        if not os.path.exists(directory):
            os.mkdir(directory)
        try:
            CCT.to_csv(s+'output'+s+filename,index=False)
        except:
            print('cannot save result')
    return [CCT,updata]

#返回CRISPR序列截断的菌株
def remove_low_quality_data(df):
    e_list=[]
    t_list=[]
    flag=0
    for index,row in df.iterrows():
        if pd.isna(row[7]):
            t_list.append(row['index'])
            if len(t_list)>=2:
                if t_list[-1]-t_list[-2]>2:
                    flag=1
        else:
            if flag:
                e_list.append(name)
            name=row[7]
            t_list=[]
            flag = 0
    return e_list